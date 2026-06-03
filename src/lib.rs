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
	#[error("")]
	SizeMismatch { provided: usize, expected: usize },
	#[error("")]
	NodeCountMismatch { got: usize },
	#[error("")]
	MovableCountMismatch { got: usize },
	#[error("")]
	UnknownParentLink,
	#[error("")]
	UnsupportedJointType,
	#[error("")]
	MathError,
	#[error("")]
	IkNotConverged {
		tries:    usize,
		pos_diff: Vector3<f64>,
		rot_diff: Vector3<f64>,
	},
	#[cfg(feature = "urdf")]
	#[error("")]
	UrdfError(#[from] urdf_rs::UrdfError),
}
