use std::f32::consts::{FRAC_PI_2, TAU};

use kinema::{
	Isometry3, Translation3, Unit, UnitQuaternion, Vector3,
	joint::{Joint, JointLimit, JointType},
	kinematics::Chain,
	node::Node,
};

fn main() {
	let _chain = robot_chain();
}

pub fn robot_chain() -> Chain<6, 7, f32> {
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

	let zero_lim = || JointLimit {
		min:      0.0,
		max:      0.0,
		effort:   0.0,
		velocity: 0.0,
	};
	let identity = Isometry3::identity;

	// ── nodes must be in topological order ─────────────

	//  0 · base_link — root, fixed
	let base = Node {
		parent:          None,
		joint:           Joint {
			name:       "base_link",
			joint_type: JointType::Fixed,
			pos:        0.0,
			limits:     zero_lim(),
			origin:     identity(),
		},
		world_transform: identity(),
	};

	//  1 · joint_1 → link_1  (parent: 0)
	//    origin  xyz="0 0 0.212"   rpy="0 0 0"
	//    limits  -180..180 deg  30 Nm  31.41 rad/s
	let j1 = Node {
		parent:          Some(0),
		joint:           Joint {
			name:       "joint_1",
			joint_type: revolute_z(),
			pos:        0.0,
			limits:     deg_lim(-180.0, 180.0, 30.0, 31.41),
			origin:     iso(0.0, 0.0, 0.212, 0.0, 0.0, 0.0),
		},
		world_transform: identity(),
	};

	//  2 · joint_2 → link_2  (parent: 1)
	//    origin  xyz="0.071 0 0.08003"   rpy="π/2 0 0"
	//    limits  -145..90 deg  25 Nm  6.28 rad/s
	let j2 = Node {
		parent:          Some(1),
		joint:           Joint {
			name:       "joint_2",
			joint_type: revolute_z(),
			pos:        0.0,
			limits:     deg_lim(-145.0, 90.0, 25.0, TAU),
			origin:     iso(0.071, 0.0, 0.080_030, FRAC_PI_2, 0.0, 0.0),
		},
		world_transform: identity(),
	};

	//  3 · joint_3 → link_3  (parent: 2)
	//    origin  xyz="0 0.290 0"   rpy="0 0 0"
	//    limits  -65..145 deg  1000 Nm  6.28 rad/s
	let j3 = Node {
		parent:          Some(2),
		joint:           Joint {
			name:       "joint_3",
			joint_type: revolute_z(),
			pos:        0.0,
			limits:     deg_lim(-65.0, 145.0, 1_000.0, TAU),
			origin:     iso(0.0, 0.290, 0.0, 0.0, 0.0, 0.0),
		},
		world_transform: identity(),
	};

	//  4 · joint_4 → link_4  (parent: 3)
	//    origin  xyz="0 0 0"   rpy="0 π/2 0"
	//    limits  -360..360 deg  2 Nm  11 rad/s
	let j4 = Node {
		parent:          Some(3),
		joint:           Joint {
			name:       "joint_4",
			joint_type: revolute_z(),
			pos:        0.0,
			limits:     deg_lim(-360.0, 360.0, 2.0, 11.0),
			origin:     iso(0.0, 0.0, 0.0, 0.0, FRAC_PI_2, 0.0),
		},
		world_transform: identity(),
	};

	//  5 · joint_5 → link_5  (parent: 4)
	//    origin  xyz="0 0 0.250"   rpy="0 -π/2 0"
	//    limits  -90..90 deg  3.6 Nm  43.6 rad/s
	let j5 = Node {
		parent:          Some(4),
		joint:           Joint {
			name:       "joint_5",
			joint_type: revolute_z(),
			pos:        0.0,
			limits:     deg_lim(-90.0, 90.0, 3.6, 43.6),
			origin:     iso(0.0, 0.0, 0.250, 0.0, -FRAC_PI_2, 0.0),
		},
		world_transform: identity(),
	};

	//  6 · joint_6 → link_6  (parent: 5)
	//    origin  xyz="0.059 0 0"   rpy="0 π/2 0"
	//    limits  -360..360 deg  0.125 Nm  15.7 rad/s
	let j6 = Node {
		parent:          Some(5),
		joint:           Joint {
			name:       "joint_6",
			joint_type: revolute_z(),
			pos:        0.0,
			limits:     deg_lim(-360.0, 360.0, 0.125, 15.7),
			origin:     iso(0.059, 0.0, 0.0, 0.0, FRAC_PI_2, 0.0),
		},
		world_transform: identity(),
	};

	Chain::new(
		[base, j1, j2, j3, j4, j5, j6],
		core::array::from_fn(|i| i + 1), // movable = nodes 1-6
	)
}
