use std::{env, path::PathBuf};

use proc_macro::TokenStream;
use quote::quote;
use syn::{LitStr, parse_macro_input};

#[proc_macro]
pub fn load_urdf(input: TokenStream) -> TokenStream {
	let input_lit = parse_macro_input!(input as LitStr);

	let workspace_dir = env::var("CARGO_WORKSPACE_DIR")
		.map(PathBuf::from)
		.unwrap_or_else(|_| PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap()));

	let mut urdf_path = workspace_dir;
	urdf_path.push(input_lit.value());

	let robot = urdf_rs::read_file(&urdf_path).expect("Failed to parse URDF file");

	let mut node_tokens = Vec::new();
	let mut movable_indices = Vec::new();

	let make_isometry_tokens = |x: f32, y: f32, z: f32, r: f32, p: f32, y_angle: f32| {
		let q = nalgebra::UnitQuaternion::from_euler_angles(r, p, y_angle);
		let coords = q.as_vector();
		let qx = coords.x;
		let qy = coords.y;
		let qz = coords.z;
		let qw = coords.w;

		quote! {
			nalgebra::Isometry3 {
				translation: nalgebra::Translation3 {
					vector: nalgebra::Vector3::new(#x, #y, #z),
				},
				rotation: nalgebra::Unit::new_unchecked(
					nalgebra::Quaternion::new(#qw, #qx, #qy, #qz)
				),
			}
		}
	};

	let identity_iso = make_isometry_tokens(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	node_tokens.push(quote! {
		Node {
			parent: None,
			joint: Joint {
				joint_type: JointType::Fixed,
				pos: 0.0,
				limits: JointLimit { min: 0.0, max: 0.0, effort: 0.0, velocity: 0.0 },
				origin: #identity_iso,
			},
			world_transform: #identity_iso,
		}
	});

	for (i, joint) in robot.joints.iter().enumerate() {
		let current_node_idx = i + 1;
		let parent_token = quote! { Some(#i) };

		let joint_type_token = match joint.joint_type {
			urdf_rs::JointType::Fixed => quote! { JointType::Fixed },
			urdf_rs::JointType::Revolute | urdf_rs::JointType::Continuous => {
				movable_indices.push(current_node_idx);
				let ax = joint.axis.xyz[0] as f32;
				let ay = joint.axis.xyz[1] as f32;
				let az = joint.axis.xyz[2] as f32;
				quote! {
					JointType::Revolute {
						axis: nalgebra::Unit::new_unchecked(nalgebra::Vector3::new(#ax, #ay, #az))
					}
				}
			}
			_ => unimplemented!("Unsupported joint type"),
		};

		let limit_token = {
			let min = joint.limit.lower as f32;
			let max = joint.limit.upper as f32;
			let effort = joint.limit.effort as f32;
			let vel = joint.limit.velocity as f32;
			quote! { JointLimit { min: #min, max: #max, effort: #effort, velocity: #vel } }
		};

		let x = joint.origin.xyz[0] as f32;
		let y = joint.origin.xyz[1] as f32;
		let z = joint.origin.xyz[2] as f32;
		let roll = joint.origin.rpy[0] as f32;
		let pitch = joint.origin.rpy[1] as f32;
		let yaw = joint.origin.rpy[2] as f32;

		let origin_iso = make_isometry_tokens(x, y, z, roll, pitch, yaw);

		node_tokens.push(quote! {
			Node {
				parent: #parent_token,
				joint: Joint {
					joint_type: #joint_type_token,
					pos: 0.0,
					limits: #limit_token,
					origin: #origin_iso,
				},
				world_transform: #identity_iso,
			}
		});
	}

	let num_nodes = node_tokens.len();
	let num_movable = movable_indices.len();

	let path_str = urdf_path.to_str().expect("Valid UTF-8 path");
	let expanded = quote! {
		{
			const _: &[u8] = include_bytes!(#path_str);
			Chain::<#num_movable, #num_nodes, f32>::new(
				[ #(#node_tokens),* ],
				[ #(#movable_indices),* ],
			)
		}
	};

	TokenStream::from(expanded)
}
