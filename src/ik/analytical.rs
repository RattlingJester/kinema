use nalgebra::{Isometry3, Matrix3, RealField, SVector, Vector3};
use simba::scalar::SubsetOf;

use crate::{Error, ik::IkSolver, kinematics::Chain};

// 6-DOF ARM
const DOF: usize = 6;
// Root node + 6 joints + tool joint
const JOINTS: usize = 8;

/// Closed-form analytical IK solver for 6-DOF manipulators with spherical wrist
/// Assumptions:
///   - Joint 1: base rotation about Z
///   - Joint 2: shoulder rotation about Y
///   - Joint 3: elbow rotation about Y
///   - Joints 4/5/6: spherical wrist — axes intersect at wrist center
///
/// `d` — link offsets [d1, d2, d3, d4, d5, d6] (d6 = tool length along approach axis)
/// `a` — link lengths  [a1, a2, a3, a4, a5, a6]
/// `alpha` — twist angles
#[derive(Debug, Default)]
#[cfg_attr(feature = "defmt", derive(defmt::Format))]
pub struct AnalyticalIK<T: RealField + SubsetOf<f64> + Copy> {
	pub d1: T,
	pub a1: T,
	pub a2: T,
	pub a3: T,
	pub d6: T,
}

impl<T: RealField + SubsetOf<f64> + Copy> IkSolver<DOF, JOINTS, T> for AnalyticalIK<T> {
	fn solve(&self, chain: &mut Chain<DOF, JOINTS, T>, target: Isometry3<T>) -> Result<(), Error> {
		let orig_pos = chain.joint_positions();

		let solution = self.solve_closest(&target, chain);

		match solution {
			Some(pose) => {
				chain.set_joint_positions(pose)?;
				chain.update_transforms();
				Ok(())
			}
			None => {
				chain.set_joint_positions(orig_pos)?;
				chain.update_transforms();
				Err(Error::IkNotConverged {
					tries:    1,
					pos_diff: nalgebra::try_convert(
						(target.translation.vector - chain.end_transform().translation.vector)
							.fixed_rows::<3>(0)
							.into_owned(),
					)
					.unwrap(),
					rot_diff: nalgebra::try_convert(
						chain
							.end_transform()
							.rotation
							.rotation_to(&target.rotation)
							.scaled_axis(),
					)
					.unwrap(),
				})
			}
		}
	}
}

#[derive(Debug, Clone)]
#[cfg_attr(feature = "defmt", derive(defmt::Format))]
pub struct IkSolution<T: RealField + SubsetOf<f64> + Copy> {
	pub pose:     SVector<T, 6>,
	/// false if any joint exceeds its limits
	pub feasible: bool,
}

impl<T: RealField + SubsetOf<f64> + Copy> Default for IkSolution<T> {
	fn default() -> Self {
		Self {
			pose:     SVector::zeros(),
			feasible: false,
		}
	}
}

impl<T: RealField + SubsetOf<f64> + Copy> AnalyticalIK<T> {
	pub fn new(d1: T, a1: T, a2: T, a3: T, d6: T) -> Self {
		Self { d1, a1, a2, a3, d6 }
	}

	/// Example solution selection strategy. Selects solution closest to current pose
	pub fn solve_closest(
		&self,
		target: &Isometry3<T>,
		chain: &mut Chain<DOF, JOINTS, T>,
	) -> Option<SVector<T, DOF>> {
		let current = chain.joint_positions();
		let (solutions, count) = self.solve(target, chain);

		#[cfg(feature = "debug")]
		for (i, s) in solutions[..count].iter().enumerate() {
			eprintln!(
				"solution {i}: feasible={}, joints={:.4?}",
				s.feasible,
				s.pose.as_slice()
			);
			for (j, _, node) in chain.iter_movable() {
				eprintln!(
					"  joint {j}: value={:.4}, limits=[{:.4}, {:.4}]",
					nalgebra::try_convert::<T, f32>(s.pose[j]).unwrap(),
					nalgebra::try_convert::<T, f32>(node.joint.limits.min).unwrap(),
					nalgebra::try_convert::<T, f32>(node.joint.limits.max).unwrap(),
				);
			}
		}

		solutions[..count]
			.iter()
			.filter(|s| s.feasible)
			.min_by(|a, b| {
				let da = (a.pose - current).norm();
				let db = (b.pose - current).norm();
				da.partial_cmp(&db).unwrap()
			})
			.map(|s| s.pose)
	}

	/// Solve for all feasible configurations. Returns up to 8 solutions.
	/// Chain must have update_transforms() called before this.
	pub fn solve(
		&self,
		target: &Isometry3<T>,
		chain: &Chain<DOF, JOINTS, T>,
	) -> ([IkSolution<T>; 8], usize) {
		let mut solutions: [IkSolution<T>; 8] = core::array::from_fn(|_| IkSolution::default());

		let p_tcp = target.translation.vector;

		let r_target: Matrix3<T> = target.rotation.to_rotation_matrix().into_inner();

		let approach = Vector3::new(r_target[(0, 2)], r_target[(1, 2)], r_target[(2, 2)]);

		let p_wc = p_tcp - approach * self.d6;

		let two = T::one() + T::one();

		let rxy = (p_wc.x * p_wc.x + p_wc.y * p_wc.y).sqrt();

		let r_eff = rxy - self.a1;
		let h_eff = p_wc.z - self.d1;

		let d_sw2 = r_eff * r_eff + h_eff * h_eff;
		// let d_sw = d_sw2.sqrt();

		let cos_theta3 =
			(d_sw2 - self.a2 * self.a2 - self.a3 * self.a3) / (two * self.a2 * self.a3);

		if cos_theta3.abs() > T::one() {
			return (solutions, 0);
		}

		let sin_theta3_pos = (T::one() - cos_theta3 * cos_theta3).sqrt();
		let sin_theta3_neg = -sin_theta3_pos;

		let theta1_front = T::atan2(p_wc.y, p_wc.x);
		let theta1_back = atan2_flip(p_wc.y, p_wc.x);

		let mut sol_idx = 0usize;

		for &theta1 in &[theta1_front, theta1_back] {
			let r_eff_local = if theta1 == theta1_front {
				rxy - self.a1
			} else {
				-(rxy + self.a1)
			};

			for &sin3 in &[sin_theta3_pos, sin_theta3_neg] {
				let theta3 = T::atan2(sin3, cos_theta3);

				let alpha = T::atan2(h_eff, r_eff_local);
				let beta = T::atan2(self.a3 * sin3, self.a2 + self.a3 * cos_theta3);
				let theta2 = alpha - beta;

				let r03 = rotation_0_3(theta1, theta2, theta3);
				let r36 = r03.transpose() * r_target;

				for &wrist_flip in &[false, true] {
					let theta5_pos_or_neg = T::atan2(
						(r36[(0, 2)] * r36[(0, 2)] + r36[(1, 2)] * r36[(1, 2)]).sqrt(),
						r36[(2, 2)],
					);
					let theta5 = if wrist_flip {
						-theta5_pos_or_neg
					} else {
						theta5_pos_or_neg
					};

					let (theta4, theta6) = if theta5.abs() < T::from_f64(1e-6).unwrap() {
						let t4 = T::atan2(r36[(1, 0)], r36[(0, 0)]);
						(t4, T::zero())
					} else if (theta5 - T::pi()).abs() < T::from_f64(1e-6).unwrap() {
						let t4 = T::atan2(-r36[(1, 0)], -r36[(0, 0)]);
						(t4, T::zero())
					} else {
						let sign = if wrist_flip { -T::one() } else { T::one() };
						let t4 = T::atan2(sign * r36[(1, 2)], sign * r36[(0, 2)]);
						let t6 = T::atan2(sign * r36[(2, 1)], -sign * r36[(2, 0)]);
						(t4, t6)
					};

					let pose = SVector::<T, 6>::from([
						wrap_pi(theta1),
						wrap_pi(theta2),
						wrap_pi(theta3),
						wrap_pi(theta4),
						wrap_pi(theta5),
						wrap_pi(theta6),
					]);

					let feasible = check_limits(pose, chain);

					solutions[sol_idx] = IkSolution { pose, feasible };
					sol_idx += 1;
				}
			}
		}

		(solutions, sol_idx)
	}
}

fn rotation_0_3<T: RealField + Copy>(t1: T, t2: T, t3: T) -> Matrix3<T> {
	let rz = rz(t1);
	let ry2 = ry(t2);
	let ry3 = ry(t3);
	rz * ry2 * ry3
}

#[inline]
fn rz<T: RealField + Copy>(a: T) -> Matrix3<T> {
	let (s, c) = (a.sin(), a.cos());
	Matrix3::new(
		c,
		-s,
		T::zero(),
		s,
		c,
		T::zero(),
		T::zero(),
		T::zero(),
		T::one(),
	)
}

#[inline]
fn ry<T: RealField + Copy>(a: T) -> Matrix3<T> {
	let (s, c) = (a.sin(), a.cos());
	Matrix3::new(
		c,
		T::zero(),
		s,
		T::zero(),
		T::one(),
		T::zero(),
		-s,
		T::zero(),
		c,
	)
}

fn check_limits<T, const DOF: usize, const JOINTS: usize>(
	pose: SVector<T, DOF>,
	chain: &Chain<DOF, JOINTS, T>,
) -> bool
where
	T: RealField + SubsetOf<f64> + Copy,
{
	for (idx, _, node) in chain.iter_movable() {
		let v = pose[idx];
		if v < node.joint.limits.min || v > node.joint.limits.max {
			return false;
		}
	}
	true
}

#[inline]
fn wrap_pi<T: RealField + Copy>(a: T) -> T {
	let pi = T::pi();
	let two_pi = pi.clone() + pi.clone();
	let mut v = a;
	while v >= pi {
		v = v - two_pi.clone();
	}
	while v < -pi.clone() {
		v = v + two_pi.clone();
	}
	v
}

#[inline]
fn atan2_flip<T: RealField + Copy>(y: T, x: T) -> T {
	T::atan2(-y, -x)
}
