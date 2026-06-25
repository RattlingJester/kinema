# kinema

A `no_std` robot kinematics library for bare-metal embedded systems, built on
[nalgebra](https://nalgebra.org) and largely ispired by [k](https://crates.io/crates/k) crate from openrr collection.

Provides forward kinematics and Jacobian inverse kinematics for serial-chain robot arms.

## Building a chain

### NO-STD:

Building it like that is kinda tedious, maybe I'll improve the API later.

```rust
pub fn my_robot() -> Chain<2, 4, f32> {
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

   	// 0 · base_link — root, fixed
   	let base = Node {
  		parent:     None,
  		joint:      Joint {
     			name:       "base_link".try_into().unwrap(),
     			joint_type: JointType::Fixed,
     			pos:        0.0,
     			limits:     JointLimit::default(),
     			origin:     Isometry3::identity(),
  		},
  		world_transform: Isometry3::identity(),
   	};

   	// 1 · joint_1 → link_1  (parent: 0)
   	//   origin xyz="0 0 0" rpy="0 0 0"
   	let j1 = Node {
  		parent:     Some(0),
  		joint:      Joint {
     			name:       "joint_1".try_into().unwrap(),
     			joint_type: revolute_z(),
     			pos:        0.0,
     			limits:     deg_lim(-180.0, 180.0, 30.0, 31.41),
     			origin:     Isometry3::identity(),
  		},
  		world_transform: Isometry3::identity(),
   	};

   	// 2 · joint_2 → link_2  (parent: 1)
   	//   origin xyz="0.071 0 0.292" rpy="π/2 0 0"
   	let j2 = Node {
       	parent:     Some(1),
      	joint:      Joint {
      		name:       "joint_2".try_into().unwrap(),
      		joint_type: revolute_z(),
      		pos:        0.0,
      		limits:     deg_lim(-145.0, 90.0, 25.0, TAU),
      		origin:     iso(0.071, 0.0, 0.292, FRAC_PI_2, 0.0, 0.0),
  		},
  		world_transform: Isometry3::identity(),
   	};

    // ... repeat for remaining joints ...

    let tool = Node {
        parent:     Some(2),
        joint:      Joint {
           	name:       "tool_fixed".try_into().unwrap(),
           	joint_type: JointType::Fixed,
           	pos:        0.0,
           	limits:     JointLimit::default(),
           	origin:     iso(0.0, 0.0, 0.059, 0.0, 0.0, 0.0),
        },
        world_transform: Isometry3::identity(),
    };

    Chain::new(
        [base, j1, j2, tool],
        core::array::from_fn(|i| i + 1), // movable = nodes 1-3
    )
}
```

#### Node ordering requirement

Nodes **must** be added in topological order (every parent before its
children). `update_transforms` relies on this to compute world poses in a
single forward pass without recursion.

### Load from URDF:

URDF parsing requires standard library and enabled "urdf" feature
```rust
let mut chain = Chain::from_urdf("robot.urdf").unwrap();
```
