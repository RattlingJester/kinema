use core::time::Duration;

use nalgebra::{RealField, SVector};
use simba::scalar::SubsetOf;

use crate::kinematics::Chain;

/// Synchronized joint-space trapezoidal trajectory.
#[cfg_attr(feature = "defmt", derive(defmt::Format))]
#[derive(Debug)]
pub struct Trajectory<const DOF: usize, T: RealField + SubsetOf<f64>> {
	pub profile:  TrapProfile<DOF, T>,
	pub duration: Duration,
}

/// Trapezoidal velocity profile for a multi joint move.
#[cfg_attr(feature = "defmt", derive(defmt::Format))]
#[derive(Debug, Clone)]
pub struct TrapProfile<const DOF: usize, T: RealField + SubsetOf<f64>> {
	pub start:    SVector<T, DOF>, // rad   — starting position
	pub end:      SVector<T, DOF>, // rad   — target position
	pub v_peak:   SVector<T, DOF>, // rad/s — cruise speed (or triangle apex if t_cruise == 0)
	pub t_ramp:   [Duration; DOF], // — acceleration phase
	pub t_cruise: [Duration; DOF], // — constant-velocity phase; 0 for triangular profiles
	pub duration: [Duration; DOF], // — total: 2 * t_ramp + t_cruise
}

impl<const DOF: usize, T: RealField + SubsetOf<f64> + Copy> Default for TrapProfile<DOF, T> {
	fn default() -> Self {
		Self {
			start:    SVector::<T, DOF>::zeros(),
			end:      SVector::<T, DOF>::zeros(),
			v_peak:   SVector::<T, DOF>::zeros(),
			t_ramp:   [Duration::ZERO; DOF],
			t_cruise: [Duration::ZERO; DOF],
			duration: [Duration::ZERO; DOF],
		}
	}
}

impl<const DOF: usize, T: RealField + SubsetOf<f64> + Copy> TrapProfile<DOF, T> {
	/// Compute a profile given:
	/// `v_max` - max cruise angular velocity [rad/s]
	/// `a` - acceleration [rad/s^2]
	fn compute(start: SVector<T, DOF>, end: SVector<T, DOF>, v_max: T, a: T) -> Self {
		let mut p = Self::default();

		let deltas = end - start;
		p.start = start;
		p.end = end;

		// First pass:
		// compute minimal-time profile for each joint independently
		let mut sync_time = T::zero();

		for i in 0..DOF {
			let d = deltas[i].abs();

			let d_ramp = (v_max * v_max) / a;

			let (v_peak, t_ramp, t_cruise, t_total) = if d <= d_ramp {
				// Triangular profile
				let vp = (d * a).sqrt();
				let tr = vp / a;
				let tc = T::zero();
				let tt = tr + tr;

				(vp, tr, tc, tt)
			} else {
				// Trapezoidal profile
				let tr = v_max / a;
				let dc = d - d_ramp;
				let tc = dc / v_max;
				let tt = tr + tc + tr;

				(v_max, tr, tc, tt)
			};

			let sign = deltas[i].signum();

			p.v_peak[i] = v_peak * sign;

			p.t_ramp[i] = Duration::from_secs_f64(nalgebra::convert(t_ramp));
			p.t_cruise[i] = Duration::from_secs_f64(nalgebra::convert(t_cruise));
			p.duration[i] = Duration::from_secs_f64(nalgebra::convert(t_total));

			if t_total > sync_time {
				sync_time = t_total;
			}
		}

		// Second pass:
		// stretch faster joints to match sync_time
		for i in 0..DOF {
			let d_signed = deltas[i];

			if d_signed == T::zero() {
				continue;
			}

			let d = d_signed.abs();
			let sign = d_signed.signum();

			// Solve:
			//
			// d = v*t_cruise + v^2/a
			// T = 2v/a + t_cruise
			//
			// eliminating t_cruise:
			//
			// d = v*T - v^2/a
			//
			// quadratic:
			//
			// v^2 - a*T*v + a*d = 0

			let b = -(a * sync_time);
			let c = a * d;

			let discriminant = b * b - T::from_f64(4.0).unwrap() * c;

			let sqrt_disc = discriminant.sqrt();

			let v = (-b - sqrt_disc) / T::from_f64(2.0).unwrap();

			let t_ramp = v / a;
			let t_cruise = sync_time - t_ramp - t_ramp;

			p.v_peak[i] = v * sign;

			p.t_ramp[i] = Duration::from_secs_f64(nalgebra::convert(t_ramp));

			p.t_cruise[i] = Duration::from_secs_f64(nalgebra::convert(t_cruise.max(T::zero())));

			p.duration[i] = Duration::from_secs_f64(nalgebra::convert(sync_time));
		}

		p
	}
}

impl<const DOF: usize, const JOINTS: usize, T: RealField + SubsetOf<f64> + Copy>
	Chain<DOF, JOINTS, T>
{
	/// Synchronized joint-space trapezoidal trajectory.
	/// speed is defined as a fraction of max defined for Chain `(0.0..1.0)`
	/// acc is `rad/s^2`
	pub fn jplan_trap(&self, goal: SVector<T, DOF>, speed: T, acc: T) -> Trajectory<DOF, T> {
		let start = self.joints_positions();

		let v_limit = self
			.iter_movable()
			.map(|(_, _, node)| node.joint.limits.velocity)
			.fold(T::zero(), |a, b| a.max(b))
			* speed;

		let profile = TrapProfile::compute(start, goal, v_limit, acc);

		let duration = profile
			.duration
			.iter()
			.fold(Duration::ZERO, |m, &p| m.max(p));

		Trajectory { profile, duration }
	}
}
