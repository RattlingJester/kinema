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
	pub fn new(d: [T; 6], a: [T; 6], alpha: [T; 6]) -> Self {
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
		chain: &mut Chain<6, 7, T>,
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
		chain: &mut Chain<6, 7, T>,
	) -> ([IkSolution<T>; 8], usize) {
		let mut solutions: [IkSolution<T>; 8] = core::array::from_fn(|_| IkSolution::default());
		let mut count = 0;

		let r = target.rotation.to_rotation_matrix();
		let p = target.translation.vector;

		let approach = r.matrix().column(2).into_owned();
		let wrist_center = p - approach * self.d[5];

		let theta1_a = wrist_center[1].atan2(wrist_center[0]);
		let theta1_b = theta1_a + T::pi();

		for &theta1 in &[theta1_a, theta1_b] {
			if let Some((theta2, theta3)) =
				self.solve_arm_profile_with_chain(wrist_center, theta1, chain)
			{
				let r03 = self.r03(theta1, theta2, theta3);
				let r36 = r03.transpose() * r.matrix();

				let wrist_options = Self::extract_euler_zyz_pairs(&r36);

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

	fn solve_arm_profile_with_chain(
		&self,
		target_wc: Vector3<T>,
		theta1: T,
		chain: &mut Chain<6, 7, T>,
	) -> Option<(T, T)> {
		let mut t2 = T::zero();
		let mut t3 = T::zero();

		let mut current_pose = SVector::<T, 6>::zeros();
		current_pose[0] = theta1;

		let wrist_node_idx = 4;

		for _ in 0..15 {
			current_pose[1] = t2;
			current_pose[2] = t3;

			chain.set_joint_positions_clamped(current_pose);
			chain.update_transforms();

			let active_joints = chain.joint_positions();
			let clamped_t2 = active_joints[1];
			let clamped_t3 = active_joints[2];

			let current_wc = chain.nodes[wrist_node_idx]
				.world_transform
				.translation
				.vector;

			let error = target_wc - current_wc;
			if error.norm_squared() < nalgebra::convert(1e-12) {
				return Some((clamped_t2, clamped_t3));
			}

			let full_jacobian = chain.jacobian();
			let j_sub = full_jacobian.fixed_view::<3, 2>(0, 1);

			let j_t = j_sub.transpose();
			let j_jt = j_t * j_sub;

			if let Some(inv_j_jt) = j_jt.try_inverse() {
				let delta_theta = inv_j_jt * j_t * error;
				t2 = clamped_t2 + delta_theta[0];
				t3 = clamped_t3 + delta_theta[1];
			} else {
				return None; // Singularity encountered
			}
		}

		current_pose[1] = t2;
		current_pose[2] = t3;
		chain.set_joint_positions_clamped(current_pose);
		chain.update_transforms();

		let final_wc = chain.nodes[wrist_node_idx]
			.world_transform
			.translation
			.vector;
		let final_joints = chain.joint_positions();

		if (target_wc - final_wc).norm() < nalgebra::convert(1e-4) {
			Some((final_joints[1], final_joints[2]))
		} else {
			None
		}
	}

	fn forward_arm_profile(&self, t1: T, t2: T, t3: T) -> (T, T) {
		let r03 = self.r03(t1, t2, t3);
		let wx = r03[(0, 3)] + r03[(0, 2)] * self.d[3];
		let wy = r03[(1, 3)] + r03[(1, 2)] * self.d[3];
		let wz = r03[(2, 3)] + r03[(2, 2)] * self.d[3];

		((wx * wx + wy * wy).sqrt(), wz)
	}

	fn r03(&self, t1: T, t2: T, t3: T) -> Matrix3<T> {
		let r1 = Rotation3::from_axis_angle(&Vector3::z_axis(), t1);
		let r2 = Rotation3::from_axis_angle(&Vector3::y_axis(), t2);
		let r3 = Rotation3::from_axis_angle(&Vector3::y_axis(), t3);
		(r1 * r2 * r3).into_inner()
	}

	fn extract_euler_zyz_pairs(r: &Matrix3<T>) -> [(T, T, T); 2] {
		let zero = T::zero();
		let pi = T::pi();

		if r[(2, 2)].abs() > nalgebra::convert(0.99999) {
			let t5 = if r[(2, 2)] > zero { zero } else { pi };
			let t4 = zero;
			let t6 = (-r[(0, 1)]).atan2(r[(0, 0)]);
			return [(t4, t5, t6), (t4, t5, t6)];
		}

		let t5_a = r[(2, 2)].acos();
		let t5_b = -t5_a;

		let t4_a = r[(1, 2)].atan2(r[(0, 2)]);
		let t6_a = r[(2, 1)].atan2(-r[(2, 0)]);

		let t4_b = t4_a + pi;
		let t6_b = t6_a + pi;

		[(t4_a, t5_a, t6_a), (t4_b, t5_b, t6_b)]
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
