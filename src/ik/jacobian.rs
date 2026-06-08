use nalgebra::{Isometry3, RealField, SMatrix, SVector, Vector3, Vector6};
use simba::scalar::SubsetOf;

use crate::{Error, ik::IkSolver, kinematics::Chain};

/// Numerical inverse Jacobian IK solver.
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

			let mut target_diff_6 = Vector6::<T>::zeros();
			let c_max = if J < 6 { J } else { 6 };
			for i in 0..c_max {
				target_diff_6[i] = target_diff[i].clone();
			}

			let (len_diff, rot_diff) = target_diff_to_len_rot_diff(&target_diff_6);

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

	fn iteration<const DOF: usize, const TASK_SPACE: usize>(
		&self,
		chain: &mut Chain<DOF, TASK_SPACE, T>,
		target: Isometry3<T>,
	) -> Result<SVector<T, TASK_SPACE>, Error> {
		let orig_positions = chain.joint_positions();
		let err = calc_pose_diff(&target, &chain.end_transform());
		let jacobi = chain.jacobian(); // Shape matches active parameters: SMatrix<T, TASK_SPACE, JOINTS>

		let mut jacobi_6x6 = SMatrix::<T, 6, 6>::zeros();
		let active_cols = if JOINTS < 6 { JOINTS } else { 6 };
		let active_rows = if TASK_SPACE < 6 { TASK_SPACE } else { 6 };

		for r in 0..active_rows {
			for c in 0..active_cols {
				jacobi_6x6[(r, c)] = jacobi[(r, c)].clone();
			}
		}

		let svd = jacobi_6x6.svd(true, true);
		let u = svd.u.ok_or(Error::MathError)?;
		let v_t = svd.v_t.ok_or(Error::MathError)?;
		let singular_values = svd.singular_values;

		let eps_machine = T::default_epsilon();
		let eps_factor: T = nalgebra::convert(100.0);
		let tolerance = eps_machine * eps_factor;
		let lambda_max: T = nalgebra::convert(0.15);

		let mut s_pinv = SMatrix::<T, 6, 6>::zeros();
		for i in 0..6 {
			let sigma = singular_values[i].clone();
			let lambda_sq = if sigma < tolerance {
				let ratio = sigma.clone() / tolerance.clone();
				let lambda = lambda_max.clone() * (T::one() - ratio);
				lambda.clone() * lambda
			} else {
				T::zero()
			};

			let denominator = (sigma.clone() * sigma.clone()) + lambda_sq;
			if denominator > T::zero() {
				s_pinv[(i, i)] = sigma / denominator;
			}
		}

		let jacobi_pinv_6x6 = v_t.transpose() * s_pinv * u.transpose();
		let d_q_6 = &jacobi_pinv_6x6 * &err;

		let mut d_q = SVector::<T, JOINTS>::zeros();
		for i in 0..active_cols {
			d_q[i] = d_q_6[i].clone();
		}

		if DOF > 6 {
			if let Some(f) = &self.nullpace_fn {
				let subtask = SVector::<T, DOF>::from_row_slice(&f(orig_positions.as_slice()));

				let mut j_pinv_j = SMatrix::<T, JOINTS, JOINTS>::zeros();
				for r in 0..active_cols {
					for c in 0..JOINTS {
						let mut sum = T::zero();
						for k in 0..active_rows {
							sum += jacobi_pinv_6x6[(r, k)].clone() * jacobi[(k, c)].clone();
						}
						j_pinv_j[(r, c)] = sum;
					}
				}

				for r in 0..DOF {
					let mut projector_row_sum = T::zero();
					for c in 0..JOINTS {
						let mut proj_elem = -j_pinv_j[(r, c)].clone();
						if r == c {
							proj_elem += T::one();
						}
						projector_row_sum += proj_elem * subtask[c].clone();
					}
					d_q[r] += projector_row_sum;
				}
			}
		}

		let max_joint_step: T = nalgebra::convert(0.05);
		let mut positions_vec = SVector::<T, DOF>::zeros();
		for i in 0..DOF {
			let mut delta = self.jacobian_mult.clone() * d_q[i].clone();
			if delta > max_joint_step {
				delta = max_joint_step.clone();
			} else if delta < -max_joint_step {
				delta = -max_joint_step.clone();
			}
			positions_vec[i] = orig_positions[i].clone() + delta;
		}

		chain.set_joint_positions_clamped(positions_vec);
		chain.update_transforms();

		let full_diff = calc_pose_diff(&target, &chain.end_transform());
		let mut out_diff = SVector::<T, TASK_SPACE>::zeros();
		for i in 0..active_rows {
			out_diff[i] = full_diff[i].clone();
		}
		Ok(out_diff)
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
