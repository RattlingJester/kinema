use nalgebra::{Isometry3, RealField, SVector, SimdPartialOrd};
use simba::scalar::SubsetOf;

use crate::{Error, ik::IkSolver, kinematics::Chain};

/// Trapezoidal profile for interpolated cartesian move.
#[cfg_attr(feature = "defmt", derive(defmt::Format))]
#[derive(Debug, Clone)]
pub struct CartTrap<
	const PATH_LEN: usize,
	const DOF: usize,
	const JOINTS: usize,
	T: RealField + SubsetOf<f64>,
> {
	/// [rad; DOF] - interpolated joint positions for the trajectory
	pub path:          [SVector<T, DOF>; PATH_LEN],
	/// Duration of each segment in seconds, PATH_LEN-1 entries
	pub segment_times: [T; PATH_LEN],
	/// sec - Duration of the trajectory
	pub duration:      T,
}

impl<
	const PATH_LEN: usize,
	const DOF: usize,
	const JOINTS: usize,
	T: RealField + SubsetOf<f64> + Copy,
> Default for CartTrap<PATH_LEN, DOF, JOINTS, T>
{
	fn default() -> Self {
		Self {
			path:          [SVector::zeros(); PATH_LEN],
			segment_times: [T::zero(); PATH_LEN],
			duration:      T::zero(),
		}
	}
}

impl<
	const PATH_LEN: usize,
	const DOF: usize,
	const JOINTS: usize,
	T: RealField + SubsetOf<f64> + Copy,
> CartTrap<PATH_LEN, DOF, JOINTS, T>
{
	pub fn compute(
		chain: &mut Chain<DOF, JOINTS, T>,
		start: Isometry3<T>,
		end: Isometry3<T>,
		speed_frac: T,
		acc: T,
		ik: &impl IkSolver<DOF, JOINTS, T>,
	) -> Result<Self, Error> {
		let orig_pos = chain.joint_positions();

		let mut path = [SVector::zeros(); PATH_LEN];

		for (i, item) in path.iter_mut().enumerate() {
			let s = if PATH_LEN > 1 {
				nalgebra::convert::<f64, T>(i as f64)
					/ nalgebra::convert::<f64, T>((PATH_LEN - 1) as f64)
			} else {
				T::one()
			};

			let target = Self::interpolate(&start, &end, s);

			if let Err(e) = ik.solve(chain, target) {
				#[cfg(feature = "debug")]
				eprintln!("IK failed at waypoint {i}/{PATH_LEN}: {e:?}");
				chain.set_joint_positions(orig_pos)?;
				chain.update_transforms();
				return Err(e);
			}

			*item = chain.joint_positions();
		}

		let mut duration = T::zero();
		let mut segment_times = [T::zero(); PATH_LEN];

		for (i, segment) in path.windows(2).enumerate() {
			let q0 = &segment[0];
			let q1 = &segment[1];

			let mut seg_time = T::zero();

			for (joint_idx, (_, _, node)) in chain.iter_movable().enumerate() {
				let dq = (q1[joint_idx] - q0[joint_idx]).abs();

				let t = Self::move_time(dq, node.joint.limits.velocity * speed_frac, acc);

				seg_time = seg_time.max(t);
			}

			segment_times[i] = seg_time;
			duration += seg_time;
		}

		Ok(Self {
			path,
			segment_times,
			duration,
		})
	}

	pub fn sample(&self, t: T) -> SVector<T, DOF> {
		if PATH_LEN == 0 {
			return SVector::zeros();
		}
		if PATH_LEN == 1 {
			return self.path[0];
		}

		let s = if self.duration > T::zero() {
			(t / self.duration).simd_clamp(T::zero(), T::one())
		} else {
			T::one()
		};

		let idx_f = s * nalgebra::convert::<f64, T>((PATH_LEN - 1) as f64);
		let idx_lo = nalgebra::try_convert::<T, f64>(idx_f).unwrap_or(0.0) as usize;
		let idx_hi = (idx_lo + 1).min(PATH_LEN - 1);
		let frac: T = idx_f - nalgebra::convert::<f64, T>(idx_lo as f64);

		self.path[idx_lo].lerp(&self.path[idx_hi], frac)
	}

	fn move_time(distance: T, vmax: T, amax: T) -> T {
		let two: T = nalgebra::convert(2.0);

		let d_ramp = vmax * vmax / (two * amax);

		if distance <= two * d_ramp {
			two * (distance / amax).sqrt()
		} else {
			let t_ramp = vmax / amax;
			let d_cruise = distance - two * d_ramp;

			two * t_ramp + d_cruise / vmax
		}
	}

	fn interpolate(start: &Isometry3<T>, end: &Isometry3<T>, s: T) -> Isometry3<T> {
		let translation = start.translation.vector.lerp(&end.translation.vector, s);
		let rotation = start.rotation.slerp(&end.rotation, s);
		Isometry3::from_parts(translation.into(), rotation)
	}
}
