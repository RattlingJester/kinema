use std::f32::consts::{FRAC_PI_2, TAU};

use kinema::{
	Isometry3, SVector, Translation3, Unit, UnitQuaternion, Vector3,
	joint::{Joint, JointLimit, JointType},
	kinematics::Chain,
	node::Node,
};

pub fn robot_chain() -> Chain<6, 8, f32> {
	let iso = |x, y, z, roll, pitch, yaw| {
		Isometry3::from_parts(
			Translation3::new(x, y, z),
			UnitQuaternion::from_euler_angles(roll, pitch, yaw),
		)
	};

	let revolute_z = || JointType::Revolute {
		axis: Unit::new_normalize(Vector3::z()),
	};

	let deg_lim = |lo: f32, hi: f32, effort: f32, velocity: f32| JointLimit {
		min: lo.to_radians(),
		max: hi.to_radians(),
		effort,
		velocity,
	};

	// 0. base_link — root, fixed
	let base = Node {
		parent:          None,
		joint:           Joint {
			name:       "base_link".try_into().unwrap(),
			joint_type: JointType::Fixed,
			pos:        0.0,
			limits:     JointLimit::default(),
			origin:     Isometry3::identity(),
		},
		world_transform: Isometry3::identity(),
	};

	// 1. joint_1 -> link_1  (parent: 0)
	//   origin xyz="0 0 0" rpy="0 0 0"
	let j1 = Node {
		parent:          Some(0),
		joint:           Joint {
			name:       "joint_1".try_into().unwrap(),
			joint_type: revolute_z(),
			pos:        0.0,
			limits:     deg_lim(-180.0, 180.0, 30.0, 31.41),
			origin:     Isometry3::identity(),
		},
		world_transform: Isometry3::identity(),
	};

	// 2. joint_2 -> link_2  (parent: 1)
	//   origin xyz="0.071 0 0.292" rpy="PI/2 0 0"
	let j2 = Node {
		parent:          Some(1),
		joint:           Joint {
			name:       "joint_2".try_into().unwrap(),
			joint_type: revolute_z(),
			pos:        0.0,
			limits:     deg_lim(-145.0, 90.0, 25.0, TAU),
			origin:     iso(0.071, 0.0, 0.292, FRAC_PI_2, 0.0, 0.0),
		},
		world_transform: Isometry3::identity(),
	};

	// 3. joint_3 -> link_3  (parent: 2)
	//   origin xyz="0 0.295 0" rpy="0 0 0"
	let j3 = Node {
		parent:          Some(2),
		joint:           Joint {
			name:       "joint_3".try_into().unwrap(),
			joint_type: revolute_z(),
			pos:        0.0,
			limits:     deg_lim(-65.0, 145.0, 1_000.0, TAU),
			origin:     iso(0.0, 0.295, 0.0, 0.0, 0.0, 0.0),
		},
		world_transform: Isometry3::identity(),
	};

	// 4. joint_4 -> link_4  (parent: 3)
	//   origin xyz="0.255 0 0" rpy="0 PI/2 0"
	let j4 = Node {
		parent:          Some(3),
		joint:           Joint {
			name:       "joint_4".try_into().unwrap(),
			joint_type: revolute_z(),
			pos:        0.0,
			limits:     deg_lim(-360.0, 360.0, 2.0, 11.0),
			origin:     iso(0.255, 0.0, 0.0, 0.0, FRAC_PI_2, 0.0),
		},
		world_transform: Isometry3::identity(),
	};

	// 5. joint_5 -> link_5  (parent: 4)
	//   origin xyz="0 0 0" rpy="0 -PI/2 0"
	let j5 = Node {
		parent:          Some(4),
		joint:           Joint {
			name:       "joint_5".try_into().unwrap(),
			joint_type: revolute_z(),
			pos:        0.0,
			limits:     deg_lim(-90.0, 90.0, 3.6, 43.6),
			origin:     iso(0.0, 0.0, 0.0, 0.0, -FRAC_PI_2, 0.0),
		},
		world_transform: Isometry3::identity(),
	};

	// 6. joint_6 -> link_6  (parent: 5)
	//   origin xyz="0 0 0" rpy="0 PI/2 0"
	let j6 = Node {
		parent:          Some(5),
		joint:           Joint {
			name:       "joint_6".try_into().unwrap(),
			joint_type: revolute_z(),
			pos:        0.0,
			limits:     deg_lim(-360.0, 360.0, 0.125, 15.7),
			origin:     iso(0.0, 0.0, 0.0, 0.0, FRAC_PI_2, 0.0),
		},
		world_transform: Isometry3::identity(),
	};

	// 7. tool_fixed -> tool  (parent: 6)
	//   origin xyz="0 0 0.059" rpy="0 0 0"  — fixed
	let tool = Node {
		parent:          Some(6),
		joint:           Joint {
			name:       "tool_fixed".try_into().unwrap(),
			joint_type: JointType::Fixed,
			pos:        0.0,
			limits:     JointLimit::default(),
			origin:     iso(0.0, 0.0, 0.059, 0.0, 0.0, 0.0),
		},
		world_transform: Isometry3::identity(),
	};

	Chain::new(
		[base, j1, j2, j3, j4, j5, j6, tool],
		[1, 2, 3, 4, 5, 6], // movable = nodes 1-6
	)
}

fn main() {
	let mut chain = robot_chain();

	chain
		.set_joint_positions(SVector::from([0.0, FRAC_PI_2, 0.0, 0.0, 0.0, 0.0]))
		.unwrap();
	chain.update_transforms();

	println!("TCP position = {}", chain.end_transform());
}
