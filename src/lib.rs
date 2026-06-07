#![cfg_attr(not(feature = "std"), no_std)]

pub mod ik;
pub mod joint;
pub mod kinematics;
pub mod node;
pub mod trajectory;
pub mod visual;

pub use nalgebra::{Isometry3, SMatrix, SVector, Translation3, Unit, UnitQuaternion, Vector3};

#[cfg(not(feature = "std"))]
pub(crate) const MAX_NAME_LEN: usize = 32;

#[cfg_attr(feature = "defmt", derive(defmt::Format))]
#[derive(Debug, thiserror::Error)]
pub enum Error {
	// #[error("Mismatch")]
	// SizeMismatch { provided: usize, expected: usize },
	#[error("Node count mismatch, got: {got}")]
	NodeCountMismatch { got: usize },
	#[error("Movable count mismatch, got: {got}")]
	MovableCountMismatch { got: usize },
	#[error("Unknown parent link")]
	UnknownParentLink,
	#[error("Unsupported joint type")]
	UnsupportedJointType,
	#[error("Jacobian error")]
	MathError,
	#[error(
		"Inverse kinematics not converged.\nTries: {tries},\npos_diff: {pos_diff},\nrot_diff: {rot_diff}"
	)]
	IkNotConverged {
		tries:    usize,
		pos_diff: Vector3<f64>,
		rot_diff: Vector3<f64>,
	},
	#[cfg(feature = "urdf")]
	#[error("URDF parse error: {0}")]
	UrdfError(#[from] urdf_rs::UrdfError),
}
