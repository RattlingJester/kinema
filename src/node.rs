use nalgebra::{Isometry3, RealField};

use crate::joint::Joint;

pub type NodeIDx = usize;

pub struct Node<T: RealField> {
	pub parent:          Option<NodeIDx>,
	pub joint:           Joint<T>,
	pub world_transform: Isometry3<T>,
}
