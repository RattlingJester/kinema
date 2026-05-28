use nalgebra::{Isometry3, RealField};
use simba::scalar::SubsetOf;

use crate::{kinematics::Chain, node::NodeIDx};

#[derive(Debug)]
pub struct Constraints<const J: usize> {
	pub position_x:    bool,
	pub position_y:    bool,
	pub position_z:    bool,
	pub rotation_x:    bool,
	pub rotation_y:    bool,
	pub rotation_z:    bool,
	pub inored_joints: [NodeIDx; J],
}

#[derive(Debug)]
pub struct JacobianIK<const J: usize, T: RealField, F: Fn(&[T]) -> [T; J] + Send + Sync> {
	pub allowable_error_dist:  T,
	pub allowable_error_angle: T,
	pub jacobian_mult:         T,
	pub max_tries:             usize,
	pub nullpace_fn:           Option<F>,
}

impl<const J: usize, T: RealField, F: Fn(&[T]) -> [T; J] + Send + Sync> Default
	for JacobianIK<J, T, F>
{
	fn default() -> Self {
		Self {
			allowable_error_dist:  nalgebra::convert(0.001),
			allowable_error_angle: nalgebra::convert(0.1),
			jacobian_mult:         nalgebra::convert(1.0),
			max_tries:             1000,
			nullpace_fn:           None,
		}
	}
}

impl<const J: usize, T: RealField + SubsetOf<f64>, F: Fn(&[T]) -> [T; J] + Send + Sync>
	JacobianIK<J, T, F>
{
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
			max_tries,
			nullpace_fn: None,
		}
	}

	pub fn solve<const DOF: usize>(
		&self,
		arm: Chain<DOF, J, T>,
		target: Isometry3<T>,
		constraints: &Constraints<J>,
	) -> Result<(), crate::Error> {
		let orig_pos = arm.joints_positions();

		todo!()
	}

	fn solve_with_constraints<const DOF: usize>(
		&self,
		arm: Chain<DOF, J, T>,
		target: Isometry3<T>,
		constraints: &Constraints<J>,
	) -> Result<(), crate::Error> {
		todo!()
	}
}
