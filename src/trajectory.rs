use core::time::Duration;

use nalgebra::{RealField, SVector};
use simba::scalar::SubsetOf;

use crate::trajectory::joint_trapezoidal::TrapProfile;

pub mod cart_trapezoidal;
pub mod joint_trapezoidal;

/// Synchronized joint-space trapezoidal trajectory.
#[cfg_attr(feature = "defmt", derive(defmt::Format))]
#[derive(Debug)]
pub struct Trajectory<const DOF: usize, const JOINTS: usize, T: RealField + SubsetOf<f64>> {
	pub profile:  TrapProfile<DOF, JOINTS, T>,
	pub duration: Duration,
}

impl<const DOF: usize, const JOINTS: usize, T: RealField + SubsetOf<f64> + Copy>
	Trajectory<DOF, JOINTS, T>
{
	pub fn sample(&self, t: T) -> SVector<T, DOF> {
		let half = T::from_subset(&0.5_f64);
		let two = T::from_subset(&2.0_f64);

		let t_ramp = self.profile.t_ramp;
		let t_cruise = self.profile.t_cruise;
		let t_total = t_ramp * two + t_cruise;

		let t = nalgebra::clamp(t, T::zero(), t_total);

		let p0 = &self.profile.start;
		let v_coast = &self.profile.v_coast;

		let s: T = if t <= t_ramp {
			half * t * t / t_ramp
		} else if t <= t_ramp + t_cruise {
			let tau = t - t_ramp;
			half * t_ramp + tau
		} else {
			let tau = t - t_ramp - t_cruise;
			half * t_ramp + t_cruise + tau - half * tau * tau / t_ramp
		};

		p0 + v_coast * s
	}
}
