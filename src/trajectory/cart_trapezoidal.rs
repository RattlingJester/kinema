use nalgebra::{Isometry3, RealField, SVector, SimdPartialOrd};
use simba::scalar::SubsetOf;

use crate::{
	Error,
	ik::{Constraints, JacobianIK},
	kinematics::Chain,
};

/// Trapezoidal profile for interpolated cartesian move.
#[cfg_attr(feature = "defmt", derive(defmt::Format))]
#[derive(Debug, Clone)]
pub struct CartTrap<
	const PATH_LEN: usize,
	const DOF: usize,
	const JOINTS: usize,
	T: RealField + SubsetOf<f64>,
> {
	/// rad   — starting position
	pub start: Isometry3<T>,
	/// rad   — target position
	pub end:   Isometry3<T>,
	/// [rad; DOF] - interpolated joint positions for the trajectory
	pub path:  [SVector<T, DOF>; PATH_LEN],

	pub duration: T,
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
			start:    Isometry3::identity(),
			end:      Isometry3::identity(),
			path:     [SVector::zeros(); PATH_LEN],
			duration: T::zero(),
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
	#[allow(clippy::too_many_arguments)]
	pub fn compute(
		chain: &mut Chain<DOF, JOINTS, T>,
		start: Isometry3<T>,
		end: Isometry3<T>,
		speed_frac: T,
		acc: T,
		ik: &JacobianIK<JOINTS, T>,
		constraints: Constraints<JOINTS>,
	) -> Result<Self, Error> {
		let orig_pos = chain.joint_positions();

		let two: T = nalgebra::convert(2.0);

		let v_max = chain
			.iter_movable()
			.map(|(_, _, n)| n.joint.limits.velocity)
			.fold(T::zero(), |a, b| a.max(b))
			* speed_frac;

		let path_length = start.rotation.rotation_to(&end.rotation).angle();

		let d_ramp = (v_max * v_max) / (two * acc);

		let (t_ramp, t_cruise, duration) = if path_length <= T::zero() {
			// Pure translation or no motion — fall back to a fixed duration
			// derived from translation distance treated as radians.
			let trans_len = (end.translation.vector - start.translation.vector).norm();
			let d_ramp_t = (v_max * v_max) / (two * acc);
			if two * d_ramp_t >= trans_len {
				let t_ramp = (trans_len / acc).sqrt();
				(t_ramp, T::zero(), two * t_ramp)
			} else {
				let t_ramp = v_max / acc;
				let t_cruise = (trans_len - two * d_ramp_t) / v_max;
				(t_ramp, t_cruise, two * t_ramp + t_cruise)
			}
		} else if two * d_ramp >= path_length {
			let t_ramp = (path_length / acc).sqrt();
			(t_ramp, T::zero(), two * t_ramp)
		} else {
			let t_ramp = v_max / acc;
			let t_cruise = (path_length - two * d_ramp) / v_max;
			(t_ramp, t_cruise, two * t_ramp + t_cruise)
		};

		let mut path = [SVector::<T, DOF>::zeros(); PATH_LEN];

		for (i, item) in path.iter_mut().enumerate() {
			let t = if PATH_LEN > 1 {
				duration * nalgebra::convert::<f64, T>(i as f64 / (PATH_LEN - 1) as f64)
			} else {
				duration
			};

			let t = t.simd_clamp(T::zero(), duration);
			let dist = if t <= t_ramp {
				(v_max / (two * t_ramp)) * t * t
			} else if t <= t_ramp + t_cruise {
				d_ramp + v_max * (t - t_ramp)
			} else {
				let t_dec = t - t_ramp - t_cruise;
				d_ramp + v_max * t_cruise + v_max * t_dec - (v_max / (two * t_ramp)) * t_dec * t_dec
			};
			let s = if path_length > T::zero() {
				(dist / path_length).simd_clamp(T::zero(), T::one())
			} else {
				(dist / (end.translation.vector - start.translation.vector).norm())
					.simd_clamp(T::zero(), T::one())
			};

			let target = Self::interpolate(&start, &end, s);
			#[cfg(feature = "debug")]
			eprintln!(
				"wp {i}: s={s:.6}, joints before IK: {:?}",
				chain.joint_positions()
			);

			if let Err(e) = ik.solve(chain, target, &constraints) {
				chain.set_joint_positions(orig_pos)?;

				#[cfg(feature = "debug")]
				eprintln!("IK failed at waypoint {i}/{PATH_LEN}: {e:?}");

				return Err(e);
			}

			*item = chain.joint_positions();
		}

		Ok(Self {
			start,
			end,
			path,
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

	fn interpolate(start: &Isometry3<T>, end: &Isometry3<T>, s: T) -> Isometry3<T> {
		let translation = start.translation.vector.lerp(&end.translation.vector, s);
		let rotation = start.rotation.slerp(&end.rotation, s);
		Isometry3::from_parts(translation.into(), rotation)
	}
}
