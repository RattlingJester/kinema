use nalgebra::{Isometry3, Matrix3, RealField, SVector};
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
					nalgebra::try_convert::<T, f64>(s.pose[j]).unwrap(),
					nalgebra::try_convert::<T, f64>(node.joint.limits.min).unwrap(),
					nalgebra::try_convert::<T, f64>(node.joint.limits.max).unwrap(),
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
		let mut count = 0;

		let two: T = nalgebra::convert(2.0);
		let r = target.rotation.to_rotation_matrix();
		let p = target.translation.vector;

		let z6 = r.matrix().column(2).into_owned();
		let wrist_center = p - z6 * self.d6;
		let wx = wrist_center[0];
		let wy = wrist_center[1];
		let wz = wrist_center[2];

		let theta1_a = wy.atan2(wx);
		let theta1_b = theta1_a + T::pi();

		for &theta1 in &[theta1_a, theta1_b] {
			let r_xy = wx * theta1.cos() + wy * theta1.sin() - self.a1;

			let z2 = wz - self.d1;

			let d_sw_sq = r_xy * r_xy + z2 * z2;

			let cos_theta3 =
				(d_sw_sq - self.a2 * self.a2 - self.a3 * self.a3) / (two * self.a2 * self.a3);

			if cos_theta3 < nalgebra::convert(-1.0) || cos_theta3 > T::one() {
				continue;
			}

			let sin_theta3_pos = (T::one() - cos_theta3 * cos_theta3).sqrt();

			for &sin_theta3 in &[sin_theta3_pos, -sin_theta3_pos] {
				let theta3 = sin_theta3.atan2(cos_theta3);

				let k1 = self.a2 + self.a3 * sin_theta3;
				let k2 = self.a3 * cos_theta3;
				let theta2 = r_xy.atan2(z2) - k1.atan2(k2);

				let r03 = self.calculate_r03(theta1, theta2, theta3);
				let r36 = r03.transpose() * r.matrix();

				let wrist_options = extract_wrist_yzy(&r36);

				for (t4, t5, t6) in wrist_options {
					if count < 8 {
						let mut pose = SVector::<T, 6>::zeros();
						pose[0] = theta1;
						pose[1] = theta2;
						pose[2] = theta3;
						pose[3] = t4;
						pose[4] = t5;
						pose[5] = t6;

						solutions[count] = IkSolution {
							feasible: self.check_limits(chain, &pose),
							pose,
						};
						count += 1;
					}
				}
			}
		}
		(solutions, count)
	}

	fn calculate_r03(&self, t1: T, t2: T, t3: T) -> Matrix3<T> {
		let (c1, s1) = (t1.cos(), t1.sin());
		let (c2, s2) = (t2.cos(), t2.sin());
		let (c3, s3) = (t3.cos(), t3.sin());

		let r01 = Matrix3::new(
			c1,
			-s1,
			T::zero(),
			s1,
			c1,
			T::zero(),
			T::zero(),
			T::zero(),
			T::one(),
		);

		let r12_fixed = Matrix3::new(
			T::one(),
			T::zero(),
			T::zero(),
			T::zero(),
			T::zero(),
			-T::one(),
			T::zero(),
			T::one(),
			T::zero(),
		);
		let r12_joint = Matrix3::new(
			c2,
			-s2,
			T::zero(),
			s2,
			c2,
			T::zero(),
			T::zero(),
			T::zero(),
			T::one(),
		);
		let r12 = r12_fixed * r12_joint;

		let r23 = Matrix3::new(
			c3,
			-s3,
			T::zero(),
			s3,
			c3,
			T::zero(),
			T::zero(),
			T::zero(),
			T::one(),
		);

		r01 * r12 * r23
	}

	fn check_limits<const DOF: usize, const JOINTS: usize>(
		&self,
		chain: &Chain<DOF, JOINTS, T>,
		joints: &SVector<T, DOF>,
	) -> bool {
		chain.iter_movable().all(|(i, _, node)| {
			let lim = &node.joint.limits;
			joints[i] >= lim.min && joints[i] <= lim.max
		})
	}
}

fn extract_wrist_yzy<T: RealField + Copy>(r: &Matrix3<T>) -> [(T, T, T); 2] {
	let sy = (r[(0, 1)] * r[(0, 1)] + r[(2, 1)] * r[(2, 1)]).sqrt();

	if sy < T::default_epsilon() {
		let t4 = T::zero();
		let t5 = if r[(1, 1)] > T::zero() {
			T::zero()
		} else {
			T::pi()
		};
		let t6 = r[(2, 0)].atan2(r[(0, 0)]);
		[(t4, t5, t6), (t4, t5, t6)]
	} else {
		let t4_a = r[(2, 1)].atan2(r[(0, 1)]);
		let t5_a = sy.atan2(r[(1, 1)]);
		let t6_a = r[(1, 2)].atan2(-r[(1, 0)]);

		let t4_b = (-r[(2, 1)]).atan2(-r[(0, 1)]);
		let t5_b = (-sy).atan2(r[(1, 1)]);
		let t6_b = (-r[(1, 2)]).atan2(r[(1, 0)]);

		[(t4_a, t5_a, t6_a), (t4_b, t5_b, t6_b)]
	}
}
