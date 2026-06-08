use nalgebra::{Isometry3, RealField, SMatrix, SVector, Vector3, Vector6};
use simba::scalar::SubsetOf;

use crate::{Error, ik::IkSolver, kinematics::Chain};

/// Numerical inverse Jacobian IK solver. Still WIP, took the solution from k crate by openrr
#[cfg_attr(feature = "defmt", derive(defmt::Format))]
pub struct JacobianIK<const JOINTS: usize, T: RealField> {
	pub allowable_error_dist:  T,
	pub allowable_error_angle: T,
	pub jacobian_mult:         T,
	pub max_try:               usize,
	#[allow(clippy::type_complexity)]
	pub nullpace_fn:           Option<&'static (dyn Fn(&[T]) -> [T; JOINTS] + Send + Sync)>,
}

impl<const D: usize, const J: usize, T: RealField + SubsetOf<f64> + Copy> IkSolver<D, J, T>
	for JacobianIK<J, T>
{
	fn solve(&self, chain: &mut Chain<D, J, T>, target: Isometry3<T>) -> Result<(), Error> {
		let orig_positions = chain.joint_positions();
		let mut last_target_distance = None;

		for _ in 0..self.max_try {
			let target_diff = self.iteration(chain, target)?;
			let (len_diff, rot_diff) = target_diff_to_len_rot_diff(&target_diff);
			if len_diff.norm() < self.allowable_error_dist
				&& rot_diff.norm() < self.allowable_error_angle
			{
				let non_checked_positions = chain.joint_positions();
				chain.set_joint_positions_clamped(non_checked_positions);
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
	) -> Result<Vector6<T>, Error> {
		let orig_positions = chain.joint_positions();

		let err = calc_pose_diff(&target, &chain.end_transform());

		let jacobi = chain.jacobian(); // Expected shape: smatrix<T, 6, DOF>

		let mut d_q = SVector::zeros();

		if DOF > 6 {
			// Redundant system: J_pinv = J^T * (J * J^T + lambda^2 * I)^-1
			// (J * J^T) is always a fixed 6x6 matrix regardless of your joint count!
			let eps: T = nalgebra::convert(0.01); // Damping factor to handle 3D singularities
			let lambda_sq = eps * eps;

			let mut jj_t = jacobi * jacobi.transpose(); // 6x6 Matrix
			for i in 0..6 {
				jj_t[(i, i)] += lambda_sq;
			}

			let jj_t_inv = jj_t.try_inverse().ok_or(Error::MathError)?;
			let jacobi_pinv = jacobi.transpose() * jj_t_inv; // DOF x 6 Matrix

			match &self.nullpace_fn {
				Some(f) => {
					let subtask = SVector::from_row_slice(&f(orig_positions.as_slice()));
					let identity = SMatrix::identity();
					let projector = identity - (jacobi_pinv * jacobi);
					d_q = (jacobi_pinv * err) + (projector * subtask);
				}
				None => {
					d_q = jacobi_pinv * err;
				}
			}
		} else {
			// Square or Under-determined system (DOF <= 6): J_pinv = (J^T * J + lambda^2 * I)^-1 * J^T
			// (J^T * J) is a fixed DOF x DOF matrix
			let eps: T = nalgebra::convert(0.01);
			let lambda_sq = eps * eps;

			let mut j_t_j = jacobi.transpose() * jacobi; // DOF x DOF Matrix
			for i in 0..DOF {
				j_t_j[(i, i)] += lambda_sq;
			}

			let j_t_j_inv = j_t_j.try_inverse().ok_or(Error::MathError)?;
			d_q = j_t_j_inv * jacobi.transpose() * err;
		}

		// 4. Update configuration space
		let mut positions_vec = SVector::zeros();
		for i in 0..DOF {
			positions_vec[i] = orig_positions[i] + self.jacobian_mult * d_q[i];
		}

		chain.set_joint_positions_clamped(positions_vec);
		chain.update_transforms();

		// 5. Return target differences for convergence checks
		Ok(calc_pose_diff(&target, &chain.end_transform()))
	}
}

fn target_diff_to_len_rot_diff<T: RealField>(target_diff: &Vector6<T>) -> (Vector3<T>, Vector3<T>) {
	let len_diff = Vector3::new(
		target_diff[0].clone(),
		target_diff[1].clone(),
		target_diff[2].clone(),
	);

	let rot_diff = Vector3::new(
		target_diff[3].clone(),
		target_diff[4].clone(),
		target_diff[5].clone(),
	);

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
