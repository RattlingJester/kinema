use nalgebra::{Isometry3, RealField, SMatrix, SVector, Vector3, Vector6};
use simba::scalar::SubsetOf;

use crate::{Error, kinematics::Chain, node::NodeIDx};

// #[derive(Debug)]
#[cfg_attr(feature = "defmt", derive(defmt::Format))]
pub struct JacobianIK<const JOINTS: usize, T: RealField> {
	pub allowable_error_dist:  T,
	pub allowable_error_angle: T,
	pub jacobian_mult:         T,
	pub max_try:               usize,
	#[allow(clippy::type_complexity)]
	pub nullpace_fn:           Option<&'static (dyn Fn(&[T]) -> [T; JOINTS] + Send + Sync)>,
}

impl<const J: usize, T: RealField + SubsetOf<f64> + Copy> Default for JacobianIK<J, T> {
	fn default() -> Self {
		Self::new(
			nalgebra::convert(0.001),
			nalgebra::convert(0.1),
			nalgebra::convert(1.0),
			1000,
		)
	}
}

impl<const JOINTS: usize, T: RealField + SubsetOf<f64> + Copy> JacobianIK<JOINTS, T> {
	pub fn new(
		allowable_error_dist: T,
		allowable_error_angle: T,
		jacobian_mult: T,
		max_tries: usize,
	) -> Self {
		Self {
			allowable_error_dist,
			allowable_error_angle,
			jacobian_mult,
			max_try: max_tries,
			nullpace_fn: None,
		}
	}

	fn iteration<const DOF: usize>(
		&self,
		chain: &mut Chain<DOF, JOINTS, T>,
		target: Isometry3<T>,
		op_space: [bool; 6],
		ignored_joints: &[usize],
	) -> Result<SVector<T, DOF>, Error> {
		let required_dof = op_space.iter().filter(|x| **x).count();
		let orig_positions = chain.joint_positions();
		let available_dof = DOF - ignored_joints.len();

		let t_n = chain.end_transform();
		let err_full = calc_pose_diff_with_constraints::<DOF, T>(&target, &t_n, op_space);

		let mut err = SVector::<T, 6>::zeros();
		for (i, &use_i) in op_space.iter().enumerate() {
			if use_i {
				err[i] = err_full[i];
			}
		}

		let jacobi_full = chain.jacobian();
		let mut jacobi = SMatrix::<T, 6, DOF>::zeros();

		for src_r in 0..6 {
			if !op_space[src_r] {
				continue;
			}
			for src_c in 0..DOF {
				if ignored_joints.contains(&src_c) {
					continue;
				}
				jacobi[(src_r, src_c)] = jacobi_full[(src_r, src_c)];
			}
		}

		let mut d_q = SVector::<T, DOF>::zeros();

		if available_dof > required_dof {
			// Redundant system: Solve via Damped Least Squares (DLS) Pseudo-Inverse
			// Formula: J_pinv = J^T * (J * J^T + lambda^2 * I)^-1
			let eps: T = nalgebra::convert(0.0001);
			let lambda_sq = eps * eps;

			let mut jj_t = jacobi * jacobi.transpose();

			for i in 0..6 {
				if op_space[i] {
					jj_t[(i, i)] += lambda_sq;
				} else {
					jj_t[(i, i)] = T::one();
				}
			}

			let jj_t_inv = jj_t.try_inverse().ok_or(Error::MathError)?;
			let jacobi_pinv = jacobi.transpose() * jj_t_inv;

			match self.nullpace_fn {
				Some(ref f) => {
					let subtask_full = f(orig_positions.as_slice());
					let mut subtask = SVector::<T, DOF>::zeros();
					for src_c in 0..DOF {
						if !ignored_joints.contains(&src_c) {
							subtask[src_c] = subtask_full[src_c];
						}
					}

					let identity = SMatrix::<T, DOF, DOF>::identity();
					d_q = jacobi_pinv * err + (identity - jacobi_pinv * jacobi) * subtask;
				}
				None => {
					d_q = jacobi_pinv * err;
				}
			}
		} else {
			// Square or Over-constrained system: Solve via Normal Equations
			// Formula: (J^T * J)^-1 * J^T * err
			// j_t_j: DOF x DOF matrix
			let mut j_t_j = jacobi.transpose() * jacobi;

			for i in 0..DOF {
				if ignored_joints.contains(&i) {
					j_t_j[(i, i)] = T::one();
				}
			}

			let j_t_j_inv = j_t_j.try_inverse().ok_or(Error::MathError)?;
			d_q = j_t_j_inv * jacobi.transpose() * err;
		};

		let mut positions_vec = SVector::<T, DOF>::zeros();
		for i in 0..DOF {
			if ignored_joints.contains(&i) {
				positions_vec[i] = orig_positions[i];
			} else {
				positions_vec[i] = orig_positions[i] + self.jacobian_mult * d_q[i];
			}
		}

		chain.set_joint_positions_clamped(positions_vec)?;
		chain.update_transforms();

		Ok(calc_pose_diff_with_constraints(
			&target,
			&chain.end_transform(),
			op_space,
		))
	}

	pub fn solve<const DOF: usize>(
		&self,
		chain: &mut Chain<DOF, JOINTS, T>,
		target: Isometry3<T>,
		constraints: &Constraints<JOINTS>,
	) -> Result<(), Error> {
		let op_space = constraints.operational_space();
		let orig_positions = chain.joint_positions();

		let mut last_target_distance = None;

		let mut ignored_joints = [0_usize; JOINTS];
		for (count, joint) in constraints.ignored_joints.iter().flatten().enumerate() {
			ignored_joints[count] = *joint;
		}

		for _ in 0..self.max_try {
			let target_diff = self.iteration(chain, target, op_space, &ignored_joints)?;
			let (len_diff, rot_diff) = target_diff_to_len_rot_diff(&target_diff, op_space);
			if len_diff.norm() < self.allowable_error_dist
				&& rot_diff.norm() < self.allowable_error_angle
			{
				let non_checked_positions = chain.joint_positions();
				chain.set_joint_positions_clamped(non_checked_positions)?;
				chain.update_transforms();
				return Ok(());
			}
			last_target_distance = Some((len_diff, rot_diff));
		}

		chain.set_joint_positions(orig_positions)?;
		chain.update_transforms();

		Err(Error::IkNotConverged {
			tries:    self.max_try,
			pos_diff: nalgebra::try_convert(last_target_distance.as_ref().unwrap().0).unwrap(),
			rot_diff: nalgebra::try_convert(last_target_distance.as_ref().unwrap().1).unwrap(),
		})
	}
}

fn target_diff_to_len_rot_diff<const DOF: usize, T: RealField>(
	target_diff: &SVector<T, DOF>,
	operational_space: [bool; 6],
) -> (Vector3<T>, Vector3<T>) {
	let mut len_diff = Vector3::zeros();
	let mut index = 0;
	for i in 0..3 {
		if operational_space[i] {
			len_diff[i] = target_diff[index].clone();
			index += 1;
		}
	}
	let mut rot_diff = Vector3::zeros();
	for i in 0..3 {
		if operational_space[i + 3] {
			rot_diff[i] = target_diff[index].clone();
			index += 1;
		}
	}
	(len_diff, rot_diff)
}

fn calc_pose_diff<T: RealField>(a: &Isometry3<T>, b: &Isometry3<T>) -> Vector6<T> {
	let p_diff = a.translation.vector.clone() - b.translation.vector.clone();
	let w_diff = b.rotation.rotation_to(&a.rotation).scaled_axis();
	Vector6::new(
		p_diff[0].clone(),
		p_diff[1].clone(),
		p_diff[2].clone(),
		w_diff[0].clone(),
		w_diff[1].clone(),
		w_diff[2].clone(),
	)
}

fn calc_pose_diff_with_constraints<const DOF: usize, T: RealField>(
	a: &Isometry3<T>,
	b: &Isometry3<T>,
	operational_space: [bool; 6],
) -> SVector<T, DOF> {
	let full_diff = calc_pose_diff(a, b);
	let mut diff = SVector::from_element(T::zero());
	let mut index = 0;
	for (i, use_i) in operational_space.iter().enumerate() {
		if *use_i {
			diff[index] = full_diff[i].clone();
			index += 1;
		}
	}
	diff
}

/// Utility function to create nullspace function using reference joint positions.
///
/// H(q) = 1/2(q-q^)T W (q-q^)
/// dH(q) / dq = W (q-q^)
///
/// Taken from k crate
pub fn create_reference_positions_nullspace_function<const DOF: usize, T: RealField + Copy>(
	reference_positions: SVector<T, DOF>,
	weight_vector: SVector<T, DOF>,
) -> impl Fn(&[T]) -> [T; DOF] {
	move |positions| {
		let mut derivative_vec = [T::zero(); DOF];
		for i in 0..DOF {
			derivative_vec[i] = weight_vector[i] * (positions[i] - reference_positions[i]);
		}
		derivative_vec
	}
}

#[derive(Debug)]
pub struct Constraints<const J: usize> {
	pub position_x:     bool,
	pub position_y:     bool,
	pub position_z:     bool,
	pub rotation_x:     bool,
	pub rotation_y:     bool,
	pub rotation_z:     bool,
	pub ignored_joints: [Option<NodeIDx>; J],
}

impl<const J: usize> Default for Constraints<J> {
	fn default() -> Self {
		Constraints {
			position_x:     true,
			position_y:     true,
			position_z:     true,
			rotation_x:     true,
			rotation_y:     true,
			rotation_z:     true,
			ignored_joints: [None; J],
		}
	}
}

impl<const J: usize> Constraints<J> {
	fn operational_space(&self) -> [bool; 6] {
		let mut arr = [true; 6];
		arr[0] = self.position_x;
		arr[1] = self.position_y;
		arr[2] = self.position_z;
		arr[3] = self.rotation_x;
		arr[4] = self.rotation_y;
		arr[5] = self.rotation_z;
		arr
	}
}
