use nalgebra::{Isometry3, RealField};
use simba::scalar::SubsetOf;

#[cfg(feature = "visuals")]
use crate::visual::Visual;

use crate::joint::Joint;

pub type NodeIDx = usize;

#[cfg_attr(feature = "defmt", derive(defmt::Format))]
#[derive(Debug)]
pub struct Node<T: RealField + SubsetOf<f64>> {
	pub parent:          Option<NodeIDx>,
	pub joint:           Joint<T>,
	pub world_transform: Isometry3<T>,
	#[cfg(feature = "visuals")]
	pub link:            Option<Link<T>>,
}

#[cfg(feature = "visuals")]
#[cfg_attr(feature = "defmt", derive(defmt::Format))]
#[derive(Debug)]
pub struct Link<T: RealField + simba::scalar::SubsetOf<f64>> {
	pub name:    String,
	pub visuals: Vec<Visual<T>>,
}
