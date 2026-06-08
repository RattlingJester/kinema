use nalgebra::{RealField, SVector};
use simba::scalar::SubsetOf;

/// Closed-form analytical IK solver for 6-DOF manipulators with spherical wrist
#[derive(Debug)]
#[cfg_attr(feature = "defmt", derive(defmt::Format))]
pub struct AnalyticalIK<const DOF: usize, T: RealField + SubsetOf<f64> + Copy> {
	pub solutions: [SVector<T, DOF>; 8],
	pub count:     usize,
}
