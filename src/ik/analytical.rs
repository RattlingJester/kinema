use nalgebra::{Isometry3, Matrix3, RealField, Rotation3, SVector, Vector3};
use simba::scalar::SubsetOf;

use crate::{Error, ik::IkSolver, kinematics::Chain};

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
	/// DH parameter: link offsets
	pub d:     [T; 6],
	/// DH parameter: link lengths
	pub a:     [T; 6],
	/// DH parameter: twist angles
	pub alpha: [T; 6],
}

impl<T: RealField + SubsetOf<f64> + Copy> IkSolver<6, 7, T> for AnalyticalIK<T> {
	fn solve(&self, chain: &mut Chain<6, 7, T>, target: Isometry3<T>) -> Result<(), Error> {
		let orig_pos = chain.joint_positions();

		let solution = self.solve_closest(&target, &*chain);

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
	pub fn new(d: [T; 6], a: [T; 6], alpha: [T; 6]) -> Self {
		Self { d, a, alpha }
	}

	pub fn from_chain(chain: &Chain<6, 7, T>) -> Self {
		let mut d = [T::zero(); 6];
		let mut a = [T::zero(); 6];
		let mut alpha = [T::zero(); 6];

		for (i, _, node) in chain.iter_movable() {
			let origin = &node.joint.origin;

			// Undo the joint origin rotation to get translation in DH frame
			let t_dh = origin.rotation.inverse() * origin.translation.vector;

			a[i] = t_dh[0];
			d[i] = t_dh[2];

			let r = origin.rotation.to_rotation_matrix();
			alpha[i] = (-r[(1, 2)]).atan2(r[(2, 2)]);
		}

		eprintln!(
			"d:     {:.4?}",
			d.map(|v| nalgebra::try_convert::<T, f64>(v).unwrap())
		);
		eprintln!(
			"a:     {:.4?}",
			a.map(|v| nalgebra::try_convert::<T, f64>(v).unwrap())
		);
		eprintln!(
			"alpha: {:.4?}",
			alpha.map(|v| nalgebra::try_convert::<T, f64>(v).unwrap())
		);

		Self { d, a, alpha }
	}

	// pub fn from_chain<const DOF: usize, const JOINTS: usize>(
	// 	chain: &Chain<DOF, JOINTS, T>,
	// ) -> Self {
	// 	let mut d = [T::zero(); 6];
	// 	let mut a = [T::zero(); 6];
	// 	let mut alpha = [T::zero(); 6];

	// 	for (i, _, node) in chain.iter_movable() {
	// 		let origin = &node.joint.origin;

	// 		// a_i: translation along X of the joint origin
	// 		a[i] = origin.translation.vector[0];
	// 		// d_i: translation along Z of the joint origin
	// 		d[i] = origin.translation.vector[2];
	// 		// alpha_i: rotation about X (twist angle)
	// 		let euler = origin.rotation.euler_angles();
	// 		alpha[i] = euler.0; // roll = rotation about X
	// 	}

	// 	#[cfg(feature = "debug")]
	// 	{
	// 		eprintln!("d:     {:?}", d.map(|v| v));
	// 		eprintln!("a:     {:?}", a.map(|v| v));
	// 		eprintln!("alpha: {:?}", alpha.map(|v| v));
	// 		eprintln!("wrist center from start: ...");
	// 	}

	// 	Self { d, a, alpha }
	// }

	/// Example solution selection strategy. Selects solution closest to current pose
	pub fn solve_closest(
		&self,
		target: &Isometry3<T>,
		chain: &Chain<6, 7, T>,
	) -> Option<SVector<T, 6>> {
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
		chain: &Chain<6, 7, T>,
	) -> ([IkSolution<T>; 8], usize) {
		let mut solutions: [IkSolution<T>; 8] = core::array::from_fn(|_| IkSolution::default());
		let mut count = 0;

		let two: T = nalgebra::convert(2.0);

		let r = target.rotation.to_rotation_matrix();
		let p = target.translation.vector;

		let approach = r.matrix().column(2).into_owned();
		let wrist_center = p - approach * self.d[5];

		let wx = wrist_center[0];
		let wy = wrist_center[1];
		let wz = wrist_center[2];

		#[cfg(feature = "debug")]
		eprintln!("wrist_center: [{wx:.4}, {wy:.4}, {wz:.4}]");

		let theta1_a = wy.atan2(wx);
		let theta1_b = theta1_a + T::pi();

		for &theta1 in &[theta1_a, theta1_b] {
			let r_xy = (wx * wx + wy * wy).sqrt() - self.a[0];
			let z2 = wz - self.d[0];

			let a2 = self.a[1];
			let a3 = self.a[2];

			let d_sw_sq = r_xy * r_xy + z2 * z2;

			let cos_theta3 = (d_sw_sq - a2 * a2 - a3 * a3) / (two * a2 * a3);

			#[cfg(feature = "debug")]
			eprintln!(
				"theta1={theta1:.4}, r_xy={r_xy:.4}, z2={z2:.4}, d_sw={:.4}, cos_theta3={cos_theta3:.4}",
				d_sw_sq.sqrt()
			);

			if cos_theta3 < nalgebra::convert(-1.0) || cos_theta3 > T::one() {
				continue;
			}

			let sin_theta3_pos = (T::one() - cos_theta3 * cos_theta3).sqrt();

			for &sin_theta3 in &[sin_theta3_pos, -sin_theta3_pos] {
				let theta3 = sin_theta3.atan2(cos_theta3);

				let k1 = a2 + a3 * cos_theta3;
				let k2 = a3 * sin_theta3;
				let theta2 = r_xy.atan2(z2) - k2.atan2(k1);

				let r03 = self.r03(theta1, theta2, theta3);
				let r36 = r03.transpose() * r.matrix();

				let (theta4, theta5, theta6) = euler_zyz(&r36);

				for &t5 in &[theta5, -theta5] {
					let (t4, t6) = if t5.abs() < nalgebra::convert(1e-6_f64) {
						(T::zero(), (-r36[(0, 1)]).atan2(r36[(0, 0)]))
					} else if t5 == theta5 {
						(theta4, theta6)
					} else {
						(theta4 + T::pi(), theta6 + T::pi())
					};

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

	fn r03(&self, t1: T, t2: T, t3: T) -> Matrix3<T> {
		let r1 = Rotation3::from_axis_angle(&Vector3::z_axis(), t1);
		let r2 = Rotation3::from_axis_angle(&Vector3::y_axis(), t2);
		let r3 = Rotation3::from_axis_angle(&Vector3::y_axis(), t3);
		(r1 * r2 * r3).into_inner()
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

/// Extract ZYZ Euler angles from a 3×3 rotation matrix.
/// R = Rz(a) * Ry(b) * Rz(c)
fn euler_zyz<T: RealField + Copy>(r: &Matrix3<T>) -> (T, T, T) {
	let sy = (r[(0, 2)] * r[(0, 2)] + r[(1, 2)] * r[(1, 2)]).sqrt();

	// let singular = sy < nalgebra::convert(1e-6_f64);
	let singular = sy < T::default_epsilon();

	if singular {
		// Gimbal lock: b ≈ 0, only a+c is determined
		let a = T::zero();
		let b = T::zero();
		let c = (-r[(0, 1)]).atan2(r[(0, 0)]);
		(a, b, c)
	} else {
		let a = r[(1, 2)].atan2(r[(0, 2)]);
		let b = sy.atan2(r[(2, 2)]);
		let c = r[(2, 1)].atan2(-r[(2, 0)]);
		(a, b, c)
	}
}
