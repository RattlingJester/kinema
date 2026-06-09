use nalgebra::{Isometry3, RealField, SMatrix, SVector};
#[cfg(feature = "urdf")]
use nalgebra::{Translation3, Unit, UnitQuaternion, Vector3};
use simba::scalar::SubsetOf;

use crate::{
	Error,
	joint::{JointLimit, JointType},
	node::{Node, NodeIDx},
};
#[cfg(feature = "urdf")]
use crate::{joint::Joint, node::Link};

/// DOF   = number of movable joints
/// JOINTS = DOF + 1 (root node counts too)
#[cfg_attr(feature = "defmt", derive(defmt::Format))]
#[derive(Debug)]
pub struct Chain<const DOF: usize, const JOINTS: usize, T: RealField + SubsetOf<f64>> {
	pub nodes:     [Node<T>; JOINTS],
	movable_nodes: [NodeIDx; DOF],
}

impl<const DOF: usize, const JOINTS: usize, T: RealField + SubsetOf<f64>> Chain<DOF, JOINTS, T> {
	#[cfg(feature = "urdf")]
	pub fn from_urdf<P: AsRef<std::path::Path>>(path: &P) -> Result<Self, Error>
	where
		Chain<DOF, JOINTS, T>: TryFrom<urdf_rs::Robot, Error = Error>,
	{
		let robot = urdf_rs::utils::read_urdf_or_xacro(path)?;

		let chain = Chain::try_from(robot)?;

		Ok(chain)
	}

	pub const fn new(nodes: [Node<T>; JOINTS], movable_nodes: [NodeIDx; DOF]) -> Self {
		Self {
			nodes,
			movable_nodes,
		}
	}

	/// Return DOF + 1 world transform (should be tool joint)
	pub fn end_transform(&self) -> Isometry3<T> {
		self.nodes[DOF + 1].world_transform.clone()
	}

	pub fn joint_positions(&self) -> SVector<T, DOF> {
		SVector::from_fn(|i, _| self.nodes[self.movable_nodes[i]].joint.pos.clone())
	}

	pub fn joint_limits(&self) -> [JointLimit<T>; DOF] {
		let mut limits = core::array::from_fn(|_| JointLimit {
			min:      T::zero(),
			max:      T::zero(),
			velocity: T::zero(),
			effort:   T::zero(),
		});

		for (idx, _, node) in self.iter_movable() {
			limits[idx] = node.joint.limits.clone();
		}

		limits
	}

	pub fn set_joint_positions(&mut self, pos: SVector<T, DOF>) -> Result<(), Error> {
		for (i, &idx) in self.movable_nodes.iter().enumerate() {
			self.nodes[idx].joint.pos = pos[i].clone();
		}

		Ok(())
	}

	pub fn set_joint_positions_clamped(&mut self, pos: SVector<T, DOF>) {
		for (i, &idx) in self.movable_nodes.iter().enumerate() {
			let limits = &self.nodes[idx].joint.limits.clone();
			self.nodes[idx].joint.pos =
				nalgebra::clamp(pos[i].clone(), limits.min.clone(), limits.max.clone());
		}
	}

	/// Recompute world transforms bottom-up.
	/// Relies on nodes being stored in topological (parent-before-child) order.
	pub fn update_transforms(&mut self) {
		for i in 0..JOINTS {
			let parent_world = self.nodes[i]
				.parent
				.map(|p| self.nodes[p].world_transform.clone())
				.unwrap_or_else(Isometry3::identity);

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

#[cfg(feature = "urdf")]
impl<const DOF: usize, const JOINTS: usize, T> TryFrom<urdf_rs::Robot> for Chain<DOF, JOINTS, T>
where
	T: RealField + SubsetOf<f64> + Copy,
{
	type Error = Error;

	fn try_from(robot: urdf_rs::Robot) -> Result<Self, Self::Error> {
		use std::collections::HashMap;

		use crate::visual::Visual;

		if robot.links.len() != JOINTS {
			return Err(Error::NodeCountMismatch {
				got: robot.links.len(),
			});
		}

		let movable_count = robot
			.joints
			.iter()
			.filter(|j| !matches!(j.joint_type, urdf_rs::JointType::Fixed))
			.count();

		if movable_count != DOF {
			return Err(Error::MovableCountMismatch { got: movable_count });
		}

		let mut link_visuals: HashMap<String, Vec<Visual<T>>> = robot
			.links
			.iter()
			.map(|l| {
				let visuals = l.visual.iter().map(Visual::from_urdf).collect();
				(l.name.clone(), visuals)
			})
			.collect();

		let mut link_to_node: HashMap<String, usize> = HashMap::with_capacity(JOINTS);

		let mut nodes: Vec<Node<T>> = Vec::with_capacity(JOINTS);
		let mut movable: Vec<usize> = Vec::with_capacity(DOF);

		let root_link = robot
			.links
			.first()
			.ok_or(Error::NodeCountMismatch { got: 0 })?;

		let root_visuals = link_visuals.remove(&root_link.name).unwrap_or_default();

		nodes.push(Node {
			parent:          None,
			joint:           Joint {
				name:       format!("{}_root", root_link.name),
				joint_type: JointType::Fixed,
				pos:        T::zero(),
				limits:     JointLimit {
					min:      T::zero(),
					max:      T::zero(),
					effort:   T::zero(),
					velocity: T::zero(),
				},
				origin:     Isometry3::identity(),
			},
			link:            Some(Link {
				name:    root_link.name.clone(),
				visuals: root_visuals,
			}),
			world_transform: Isometry3::identity(),
		});

		link_to_node.insert(root_link.name.clone(), 0);

		for j in &robot.joints {
			let xyz = Vector3::new(
				nalgebra::convert(j.origin.xyz[0]),
				nalgebra::convert(j.origin.xyz[1]),
				nalgebra::convert(j.origin.xyz[2]),
			);
			let rotation = UnitQuaternion::from_euler_angles(
				nalgebra::convert(j.origin.rpy[0]),
				nalgebra::convert(j.origin.rpy[1]),
				nalgebra::convert(j.origin.rpy[2]),
			);
			let origin = Isometry3::from_parts(Translation3::from(xyz), rotation);

			let axis = Unit::new_normalize(Vector3::new(
				nalgebra::convert(j.axis.xyz[0]),
				nalgebra::convert(j.axis.xyz[1]),
				nalgebra::convert(j.axis.xyz[2]),
			));

			let joint_type = match &j.joint_type {
				urdf_rs::JointType::Fixed => JointType::Fixed,
				urdf_rs::JointType::Revolute => JointType::Revolute { axis },
				urdf_rs::JointType::Prismatic => JointType::Prismatic { axis },
				_ => return Err(Error::UnsupportedJointType),
			};

			let parent = link_to_node
				.get(&j.parent.link)
				.copied()
				.ok_or(Error::UnknownParentLink)?;

			let child_visuals = link_visuals.remove(&j.child.link).unwrap_or_default();

			let node_idx = nodes.len();
			let is_movable = !matches!(joint_type, JointType::Fixed);

			nodes.push(Node {
				parent:          Some(parent),
				joint:           Joint {
					name: j.name.clone(),
					joint_type,
					pos: T::zero(),
					limits: JointLimit {
						min:      nalgebra::convert(j.limit.lower),
						max:      nalgebra::convert(j.limit.upper),
						velocity: nalgebra::convert(j.limit.velocity),
						effort:   nalgebra::convert(j.limit.effort),
					},
					origin,
				},
				link:            Some(Link {
					name:    j.child.link.clone(),
					visuals: child_visuals,
				}),
				world_transform: Isometry3::identity(),
			});

			if is_movable {
				movable.push(node_idx);
			}

			link_to_node.insert(j.child.link.clone(), node_idx);
		}

		let nodes_arr: [Node<T>; JOINTS] = nodes
			.try_into()
			.unwrap_or_else(|_| unreachable!("length validated in step 1"));

		let movable_arr: [usize; DOF] = movable
			.try_into()
			.unwrap_or_else(|_| unreachable!("length validated in step 1"));

		Ok(Chain::new(nodes_arr, movable_arr))
	}
}
