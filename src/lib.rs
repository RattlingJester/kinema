#![no_std]

pub mod ik;
pub mod joint;
pub mod kinematics;
pub mod node;
pub mod trajectory;

pub use nalgebra::{Isometry3, SMatrix, Translation3, Unit, UnitQuaternion, Vector3};

#[cfg_attr(feature = "defmt", derive(defmt::Format))]
#[derive(Debug)]
pub enum Error {
	SizeMismatch { provided: usize, expected: usize },
}
