use core::time::Duration;

use nalgebra::RealField;
use simba::scalar::SubsetOf;

use crate::kinematics::Chain;

/// Synchronized joint-space trapezoidal trajectory.
#[cfg_attr(feature = "defmt", derive(defmt::Format))]
#[derive(Debug)]
pub struct Trajectory<const DOF: usize, T: RealField + SubsetOf<f64>> {
	pub profiles: [TrapProfile<T>; DOF],
	pub duration: Duration,
}

/// Trapezoidal velocity profile for a single joint.
#[cfg_attr(feature = "defmt", derive(defmt::Format))]
#[derive(Debug, Clone, Copy)]
pub struct TrapProfile<T: RealField + SubsetOf<f64>> {
	pub start:    T,        // rad   — starting position
	pub end:      T,        // rad   — target position
	pub v_peak:   T,        // rad/s — cruise speed (or triangle apex if t_cruise == 0)
	pub t_ramp:   T,        // s     — acceleration phase
	pub t_cruise: T,        // s     — constant-velocity phase; 0 for triangular profiles
	pub duration: Duration, // s     — total: 2 * t_ramp + t_cruise
}

impl<T: RealField + SubsetOf<f64>> Default for TrapProfile<T> {
	fn default() -> Self {
		Self {
			start:    T::zero(),
			end:      T::zero(),
			v_peak:   T::zero(),
			t_ramp:   T::zero(),
			t_cruise: T::zero(),
			duration: Duration::ZERO,
		}
	}
}

impl<T: RealField + SubsetOf<f64> + Copy> TrapProfile<T> {
	/// Compute a profile given physical acceleration `a` [rad/s^2].
	fn compute(start: T, end: T, v_max: T, a: T) -> Self {
		let delta = (end - start).abs();

		// Yeah, funny
		let two: T = nalgebra::convert(2.0);

		if delta <= nalgebra::convert(f64::EPSILON) || v_max <= nalgebra::convert(f64::EPSILON) {
			return Self::default();
		}

		// Accel distance: d = v^2 / (2 * a)
		let d_ramp = v_max * v_max / (two * a);
		if two * d_ramp >= delta {
			// Triangular profile
			// sigma = v_peak^2 / a  →  v_peak = sqrt(a * sigma)
			let v_peak = (a * delta).sqrt();
			let t_ramp = v_peak / a;
			Self {
				start,
				end,
				v_peak,
				t_ramp,
				t_cruise: T::zero(),
				duration: Duration::from_secs_f64(nalgebra::convert(two * t_ramp)),
			}
		} else {
			let t_ramp = v_max / a;
			let t_cruise = (delta - two * d_ramp) / v_max;
			Self {
				start,
				end,
				v_peak: v_max,
				t_ramp,
				t_cruise,
				duration: Duration::from_secs_f64(nalgebra::convert(two * t_ramp + t_cruise)),
			}
		}
	}

	/// Re-plan to fill exactly `target` time by lowering cruise speed.
	/// Closed-form: T = v/a + sigma/v  →  v^2 − T * a * v + a * sigma = 0
	/// Smaller root: v = a * (T − sqrt(T^2 − 4 * sigma / a)) / 2
	fn constrain_to(self, target: Duration, a: T) -> Self {
		let delta = (self.end - self.start).abs();

		// Already constrained or no move
		if self.duration >= target || delta <= nalgebra::convert(f64::EPSILON) {
			return self;
		}

		// Yeah, funny
		let two: T = nalgebra::convert(2.0);
		let four: T = nalgebra::convert(4.0);

		let target: T = nalgebra::convert(target.as_secs_f64());

		let discriminant = (target * target - four * delta / a).max(T::zero());
		let v = a * (target - (discriminant).sqrt()) / two;
		Self::compute(self.start, self.end, v.max(T::zero()), a)
	}
}

impl<const DOF: usize, const JOINTS: usize, T: RealField + SubsetOf<f64> + Copy>
	Chain<DOF, JOINTS, T>
{
	/// Synchronized joint-space trapezoidal trajectory.
	/// speed is defined as a fraction of max defined for Chain `(0.0..1.0)`
	/// acc is `rad/s^2`
	pub fn jplan_trap(&self, goal: &[T], speed: T, acc: T) -> Trajectory<DOF, T> {
		let start = self.joints_positions();

		let mut profiles: [TrapProfile<T>; DOF] = core::array::from_fn(|_| TrapProfile::default());

		for (i, _, node) in self.iter_movable() {
			profiles[i] =
				TrapProfile::compute(start[i], goal[i], node.joint.limits.velocity * speed, acc);
		}

		let duration = profiles
			.iter()
			.map(|p| p.duration)
			.max()
			.unwrap_or(Duration::ZERO);

		for p in &mut profiles {
			*p = p.constrain_to(duration, acc);
		}

		Trajectory { profiles, duration }
	}
}
