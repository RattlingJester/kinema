use nalgebra::{Isometry3, RealField, SMatrix, SVector};

use crate::{
	Error,
	joint::JointType,
	node::{Node, NodeIDx},
};

/// DOF   = number of movable joints
/// JOINTS = DOF + 1 (root node counts too)
#[cfg_attr(feature = "defmt", derive(defmt::Format))]
#[derive(Debug)]
pub struct Chain<const DOF: usize, const JOINTS: usize, T: RealField> {
	nodes:         [Node<T>; JOINTS],
	movable_nodes: [NodeIDx; DOF],
}

impl<const DOF: usize, const JOINTS: usize, T: RealField> Chain<DOF, JOINTS, T> {
	pub const fn new(nodes: [Node<T>; JOINTS], movable_nodes: [NodeIDx; DOF]) -> Self {
		Self {
			nodes,
			movable_nodes,
		}
	}

	pub fn end_transform(&self) -> Isometry3<T> {
		self.nodes[DOF].world_transform.clone()
	}

	pub fn joints_positions(&self) -> SVector<T, DOF> {
		SVector::from_fn(|i, _| self.nodes[self.movable_nodes[i]].joint.pos.clone())
	}

	pub fn set_joints_positions(&mut self, pos: SVector<T, DOF>) -> Result<(), Error> {
		if pos.len() > DOF {
			return Err(Error::SizeMismatch {
				provided: pos.len(),
				expected: DOF,
			});
		}

		for (i, &idx) in self.movable_nodes.iter().enumerate() {
			self.nodes[idx].joint.pos = pos[i].clone();
		}

		Ok(())
	}

	/// Recompute world transforms bottom-up.
	/// Relies on nodes being stored in topological (parent-before-child) order.
	pub fn update_transforms(&mut self) {
		for i in 0..JOINTS {
			let parent_world = self.nodes[i]
				.parent
				.map(|p| self.nodes[p].world_transform.clone())
				.unwrap_or_else(Isometry3::identity);

			// Split borrow: read joint fields before mutating world_transform
			let local = self.nodes[i].joint.origin.clone() * self.nodes[i].joint.local_transform();

			self.nodes[i].world_transform = parent_world * local;
		}
	}

	/// Geometric Jacobian (6 × DOF): [linear; angular]
	/// Call update_transforms() before this.
	pub fn jacobian(&self) -> SMatrix<T, 6, DOF> {
		let p_n = self.end_transform().translation.vector.clone();

		SMatrix::from_fn(|row, col| {
			let idx = self.movable_nodes[col];
			let t_i = &self.nodes[idx].world_transform;

			match &self.nodes[idx].joint.joint_type {
				JointType::Revolute { axis } => {
					let a_i = t_i.rotation.clone() * axis.clone();
					let dp_i = a_i.cross(&(p_n.clone() - t_i.translation.vector.clone()));
					// rows 0-2: linear,  rows 3-5: angular
					if row < 3 {
						dp_i[row].clone()
					} else {
						a_i[row - 3].clone()
					}
				}
				JointType::Prismatic { axis } => {
					let a_i = t_i.rotation.clone() * axis.clone();
					if row < 3 { a_i[row].clone() } else { T::zero() }
				}
				JointType::Fixed => panic!("fixed joint in movable_nodes — bug in Chain::new()"),
			}
		})
	}

	pub fn iter(&self) -> impl Iterator<Item = (NodeIDx, &Node<T>)> {
		self.nodes.iter().enumerate()
	}

	pub fn iter_movable(&self) -> impl Iterator<Item = (usize, NodeIDx, &Node<T>)> {
		self.movable_nodes
			.iter()
			.enumerate()
			.map(|(dof_idx, &id)| (dof_idx, id, &self.nodes[id]))
	}
}
