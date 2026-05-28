#![cfg(feature = "visuals")]

use simba::scalar::SubsetOf;

use nalgebra::{Isometry3, RealField, Vector3};

#[cfg_attr(feature = "defmt", derive(defmt::Format))]
#[derive(Debug)]
pub enum Geometry<T: RealField> {
	Mesh {
		filename: String,
		scale:    Vector3<T>,
	},
	Cylinder {
		radius: T,
		length: T,
	},
	Sphere {
		radius: T,
	},
	Box {
		depth:  T,
		width:  T,
		height: T,
	},
}

#[cfg_attr(feature = "defmt", derive(defmt::Format))]
#[derive(Debug, Default)]
pub struct Color {
	pub r: f32,
	pub g: f32,
	pub b: f32,
	pub a: f32,
}

#[cfg_attr(feature = "defmt", derive(defmt::Format))]
#[derive(Debug)]
pub struct Visual<T: RealField + SubsetOf<f64>> {
	pub name:     Option<String>,
	pub origin:   Isometry3<T>,
	pub geometry: Geometry<T>,
	pub color:    Color,
}

impl<T: RealField + SubsetOf<f64>> Visual<T> {
	pub(crate) fn from_urdf(v: &urdf_rs::Visual) -> Self {
		use nalgebra::{Translation3, UnitQuaternion};
		use urdf_rs::Geometry as UG;

		let conv = |x: f64| -> T {
			use nalgebra::convert;
			convert(x)
		};
		let k: T = conv(1000.0);

		let t = v.origin.xyz;
		let rpy = v.origin.rpy;
		let origin = Isometry3::from_parts(
			Translation3::new(
				conv(t[0]) * k.clone(),
				conv(t[1]) * k.clone(),
				conv(t[2]) * k.clone(),
			),
			UnitQuaternion::from_euler_angles(conv(rpy[0]), conv(rpy[1]), conv(rpy[2])),
		);

		let geometry = match &v.geometry {
			UG::Mesh { filename, scale } => Geometry::Mesh {
				filename: filename.clone(),
				scale:    if let Some(s) = scale {
					Vector3::new(conv(s[0]), conv(s[1]), conv(s[2]))
				} else {
					Vector3::new(T::one(), T::one(), T::one())
				},
			},
			UG::Cylinder { radius, length } => Geometry::Cylinder {
				radius: conv(*radius) * k.clone(),
				length: conv(*length) * k.clone(),
			},
			UG::Sphere { radius } => Geometry::Sphere {
				radius: conv(*radius) * k.clone(),
			},
			UG::Box { size } => Geometry::Box {
				depth:  conv(size[0]) * k.clone(),
				width:  conv(size[1]) * k.clone(),
				height: conv(size[2]) * k,
			},
			other => panic!(
				"Unsupported geometry type: {}",
				std::any::type_name_of_val(other)
			),
		};

		let color = if let Some(mat) = &v.material
			&& let Some(c) = &mat.color
		{
			let vec = c.rgba;
			Color {
				r: vec[0] as f32,
				g: vec[1] as f32,
				b: vec[2] as f32,
				a: vec[3] as f32,
			}
		} else {
			Color::default()
		};

		Self {
			origin,
			geometry,
			color,
			name: v.name.clone(),
		}
	}
}
