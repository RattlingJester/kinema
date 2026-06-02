#![allow(unused)]

use nalgebra::{Isometry3, RealField, SVector};
use simba::scalar::SubsetOf;

use crate::{Error, kinematics::Chain, node::NodeIDx};

#[derive(Debug)]
pub struct JacobianIK<const JOINTS: usize, T: RealField, F: Fn(&[T]) -> [T; JOINTS] + Send + Sync> {
	pub allowable_error_dist:  T,
	pub allowable_error_angle: T,
	pub jacobian_mult:         T,
	pub max_try:               usize,
	pub nullpace_fn:           Option<F>,
}

impl<const J: usize, T: RealField + SubsetOf<f64>, F: Fn(&[T]) -> [T; J] + Send + Sync> Default
	for JacobianIK<J, T, F>
{
	fn default() -> Self {
		Self::new(
			nalgebra::convert(0.001),
			nalgebra::convert(0.1),
			nalgebra::convert(1.0),
			1000,
		)
	}
}

impl<const JOINTS: usize, T: RealField + SubsetOf<f64>, F: Fn(&[T]) -> [T; JOINTS] + Send + Sync>
	JacobianIK<JOINTS, T, F>
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
			max_try: max_tries,
			nullpace_fn: None,
		}
	}

	fn solve<const DOF: usize>(
		&mut self,
		chain: &Chain<DOF, JOINTS, T>,
		constraints: Constraints<JOINTS>,
	) -> Result<SVector<T, DOF>, Error> {
		let op_space = constraints.operational_space();
		let orig_positions = chain.joint_positions();

		let required_dof = op_space.iter().filter(|x| **x).count();
		let available_dof = DOF - constraints.inored_joints.len();
		if available_dof < required_dof {
			return Err(Error::SizeMismatch {
				provided: available_dof,
				expected: required_dof,
			});
		}

		let mut last_target_distance = None;

		for _ in 0..self.max_try {}

		todo!()
	}
}

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
