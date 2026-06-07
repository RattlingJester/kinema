use core::time::Duration;

use nalgebra::{RealField, SVector};
use simba::scalar::SubsetOf;

use crate::{kinematics::Chain, trajectory::Trajectory};

/// Trapezoidal velocity profile for a multi joint move.
#[cfg_attr(feature = "defmt", derive(defmt::Format))]
#[derive(Debug, Clone)]
pub struct TrapProfile<const DOF: usize, const JOINTS: usize, T: RealField + SubsetOf<f64>> {
	/// rad   — starting position
	pub start:    SVector<T, DOF>,
	/// rad   — target position
	pub end:      SVector<T, DOF>,
	/// rad/s — cruise speed (or triangle apex if t_cruise == 0)
	pub v_coast:  SVector<T, DOF>,
	/// s — acceleration phase
	pub t_ramp:   T,
	/// s — constant-velocity phase; 0 for triangular profiles
	pub t_cruise: T,
}

impl<const DOF: usize, const JOINTS: usize, T: RealField + SubsetOf<f64> + Copy> Default
	for TrapProfile<DOF, JOINTS, T>
{
	fn default() -> Self {
		Self {
			start:    SVector::<T, DOF>::zeros(),
			end:      SVector::<T, DOF>::zeros(),
			v_coast:  SVector::<T, DOF>::zeros(),
			t_ramp:   T::zero(),
			t_cruise: T::zero(),
		}
	}
}

impl<const DOF: usize, const JOINTS: usize, T: RealField + SubsetOf<f64> + Copy>
	TrapProfile<DOF, JOINTS, T>
{
	/// Compute a profile given:
	/// `speed_frac` - fraction of joints velocity [0.0..1.0]
	/// `a` - acceleration [rad/s^2]
	fn compute(
		chain: &Chain<DOF, JOINTS, T>,
		start: SVector<T, DOF>,
		end: SVector<T, DOF>,
		speed_frac: T,
		a: T,
	) -> Self {
		let mut p = TrapProfile {
			start,
			end,
			v_coast: SVector::<T, DOF>::zeros(),
			t_ramp: T::zero(),
			t_cruise: T::zero(),
		};

		let distances = end - start;
		let deltas = distances.abs();
		let two = nalgebra::convert::<f64, T>(2.0);

		let mut local_t_ramp = [T::zero(); DOF];
		let mut local_duration = [T::zero(); DOF];

		for (idx, _id, node) in chain.iter_movable() {
			if deltas[idx].is_zero() {
				continue;
			}
			let v_max = node.joint.limits.velocity * speed_frac;

			let d_ramp = (v_max * v_max) / (two * a);
			let d_acc_dec = d_ramp * two;

			if d_acc_dec < deltas[idx] {
				// Trapezoidal
				let t_ramp = v_max / a;
				let d_coast = deltas[idx] - d_acc_dec;
				let t_coast = d_coast / v_max;

				local_t_ramp[idx] = t_ramp;
				local_duration[idx] = t_ramp * two + t_coast;
			} else {
				// Triangular
				let v_peak_tri = (deltas[idx] * a).sqrt();
				let t_ramp = v_peak_tri / a;

				local_t_ramp[idx] = t_ramp;
				local_duration[idx] = t_ramp * two;
			}
		}

		let t_total_max = local_duration
			.iter()
			.fold(T::zero(), |max_t, &t| max_t.max(t));
		let t_ramp_max = local_t_ramp
			.iter()
			.fold(T::zero(), |max_t, &t| max_t.max(t));
		let t_cruise_max = t_total_max - (t_ramp_max * two);

		if t_total_max.is_zero() {
			return p;
		}

		for (idx, _id, _node) in chain.iter_movable() {
			p.t_ramp = t_ramp_max;
			p.t_cruise = t_cruise_max;

			if !deltas[idx].is_zero() {
				let v_scaled = deltas[idx] / (t_total_max - t_ramp_max);

				p.v_coast[idx] = if distances[idx] >= T::zero() {
					v_scaled
				} else {
					-v_scaled
				};
			} else {
				p.v_coast[idx] = T::zero();
			}
		}

		p
	}
}

impl<const DOF: usize, const JOINTS: usize, T: RealField + SubsetOf<f64> + Copy>
	Chain<DOF, JOINTS, T>
{
	///Synchronized joint-space trapezoidal trajectory.
	///speed is a fraction of max angular velocity for each joint in Chain `(0.0..1.0)`
	///acc is `rad/s^2`
	pub fn jplan_trap(
		&self,
		goal: SVector<T, DOF>,
		speed: T,
		acc: T,
	) -> Trajectory<DOF, JOINTS, T> {
		let start = self.joint_positions();

		let profile = TrapProfile::compute(self, start, goal, speed, acc);

		let dur_secs = profile.t_ramp * nalgebra::convert(2.0) + profile.t_cruise;
		let duration = Duration::from_secs_f64(nalgebra::convert(dur_secs));

		Trajectory { profile, duration }
	}
}
