use nalgebra::{Isometry3, RealField};

use crate::joint::Joint;

pub type NodeIDx = usize;

#[cfg_attr(feature = "defmt", derive(defmt::Format))]
#[derive(Debug)]
pub struct Node<T: RealField> {
	pub parent:          Option<NodeIDx>,
	pub joint:           Joint<T>,
	pub world_transform: Isometry3<T>,
}
