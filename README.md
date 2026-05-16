# kinema

A `no_std` robot kinematics library for bare-metal embedded systems, built on
[nalgebra](https://nalgebra.org) and largely ispired by [k](https://crates.io/crates/k) crate from openrr collection.

Provides forward kinematics and geometric Jacobian computation for serial-chain robot arms with a fully static memory layout. 

## TODO
 - Inverse kinematics

## Building a chain
Not very convenient for now, may improve later.

```rust
use core::f32::consts::FRAC_PI_2;
use nalgebra::{Isometry3, Translation3, Unit, UnitQuaternion, Vector3};
use kinema::{Chain, Joint, JointLimit, JointType, Node, NodeIDx};

pub fn my_robot() -> Chain<6, 7, f32> {
    let iso = |x, y, z, roll, pitch, yaw| {
        Isometry3::from_parts(
            Translation3::new(x, y, z),
            UnitQuaternion::from_euler_angles(roll, pitch, yaw),
        )
    };
    let z = || JointType::Revolute { axis: Unit::new_normalize(Vector3::z()) };
    let id = Isometry3::identity;

    let base = Node {
        parent: None,
        joint: Joint {
            name: "base", joint_type: JointType::Fixed, pos: 0.0,
            limits: JointLimit { min: 0.0, max: 0.0, effort: 0.0, velocity: 0.0 },
            origin: id(),
        },
        world_transform: id(),
    };

    let j1 = Node {
        parent: Some(NodeIDx(0)),
        joint: Joint {
            name: "joint_1", joint_type: z(), pos: 0.0,
            limits: JointLimit { min: -3.14, max: 3.14, effort: 30.0, velocity: 31.41 },
            origin: iso(0.0, 0.0, 0.212, 0.0, 0.0, 0.0),
        },
        world_transform: id(),
    };

    // ... repeat for remaining joints ...

    Chain::new(
        [base, j1, /* j2, j3, j4, j5, j6 */],
        core::array::from_fn(|i| NodeIDx(i + 1)),
    )
}
```

See [`examples/six_dof.rs`](examples/six_dof.rs) for a complete 6-DOF arm.

## Usage

```rust
static CHAIN: StaticCell<Chain<6, 7, f32>> = StaticCell::new();

let chain =  CHAIN.init(my_robot());

// Set joint angles (radians)
let q = &[0.1, -0.5, 0.3, 0.0, 0.8, -0.2];
chain.set_joint_positions(q).unwrap();

// Forward kinematics
chain.update_transforms();
let end_effector: Isometry3<f32> = chain.end_transform();
```

## Node ordering requirement

Nodes **must** be added in topological order (every parent before its
children). `update_transforms` relies on this to compute world poses in a
single forward pass without recursion.
