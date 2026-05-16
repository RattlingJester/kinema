use nalgebra::{Isometry3, RealField, Translation3, Unit, UnitQuaternion, Vector3};

#[derive(Debug, Default)]
pub struct JointLimit<T: RealField> {
	/// Radians
	pub min:      T,
	/// Radians
	pub max:      T,
	/// Radians/s
	pub velocity: T,
	/// N*m
	pub effort:   T,
}

#[derive(Debug, PartialEq)]
pub enum JointType<T: RealField> {
	Fixed,
	Revolute { axis: Unit<Vector3<T>> },
	Prismatic { axis: Unit<Vector3<T>> },
}

#[derive(Debug)]
pub struct Joint<T: RealField> {
	pub name:       &'static str,
	pub joint_type: JointType<T>,
	pub pos:        T,
	pub limits:     JointLimit<T>,
	pub origin:     Isometry3<T>,
}

impl<T: RealField> Joint<T> {
	pub fn local_transform(&self) -> Isometry3<T> {
		match &self.joint_type {
			JointType::Fixed => Isometry3::identity(),
			JointType::Revolute { axis } => Isometry3::from_parts(
				Translation3::new(T::zero(), T::zero(), T::zero()),
				UnitQuaternion::from_axis_angle(axis, self.pos.clone()),
			),
			JointType::Prismatic { axis } => Isometry3::from_parts(
				Translation3::from(axis.clone().into_inner() * self.pos.clone()),
				UnitQuaternion::identity(),
			),
		}
	}
}
