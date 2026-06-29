use core::f64::consts::TAU;

use nalgebra::{Isometry3, RealField, SVector, SimdPartialOrd, Vector3};
use simba::scalar::SubsetOf;

use crate::{Error, ik::IkSolver, kinematics::Chain};

/// Trapezoidal profile for interpolated cartesian move.
#[cfg_attr(feature = "defmt", derive(defmt::Format))]
#[derive(Debug, Clone)]
pub struct CircTrap<
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
> Default for CircTrap<PATH_LEN, DOF, JOINTS, T>
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
> CircTrap<PATH_LEN, DOF, JOINTS, T>
{
	pub fn compute(
		chain: &mut Chain<DOF, JOINTS, T>,
		start: Isometry3<T>,
		via: Isometry3<T>,
		end: Isometry3<T>,
		speed_frac: T,
		acc: T,
		ik: &impl IkSolver<DOF, JOINTS, T>,
	) -> Result<Self, Error> {
		let orig_pos = chain.joint_positions();

		let arc = CircularArc::from_three_points(
			start.translation.vector,
			via.translation.vector,
			end.translation.vector,
		);

		let mut path = [SVector::zeros(); PATH_LEN];

		for (i, item) in path.iter_mut().enumerate() {
			let s = if PATH_LEN > 1 {
				nalgebra::convert::<f64, T>(i as f64)
					/ nalgebra::convert::<f64, T>((PATH_LEN - 1) as f64)
			} else {
				T::one()
			};

			let target = Self::interpolate(&arc, &start, &via, &end, s);
			if let Err(e) = ik.solve(chain, target) {
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

	fn interpolate(
		arc: &CircularArc<T>,
		start: &Isometry3<T>,
		via: &Isometry3<T>,
		end: &Isometry3<T>,
		s: T,
	) -> Isometry3<T> {
		let half: T = nalgebra::convert(0.5);
		let two: T = nalgebra::convert(2.0);

		let translation = arc.point_at(s);
		let rotation = if s <= half {
			start.rotation.slerp(&via.rotation, s * two)
		} else {
			via.rotation.slerp(&end.rotation, (s - half) * two)
		};

		Isometry3::from_parts(translation.into(), rotation)
	}
}

/// Circle defined by three points in 3D space.
struct CircularArc<T: RealField + Copy> {
	/// Circle center
	center:    Vector3<T>,
	/// Vector from center to start point
	r_start:   Vector3<T>,
	/// Normal to the plane of the circle
	normal:    Vector3<T>,
	/// Total arc angle in radians
	arc_angle: T,
}

impl<T: RealField + SubsetOf<f64> + Copy> CircularArc<T> {
	/// Create arc from three points
	fn from_three_points(p1: Vector3<T>, p2: Vector3<T>, p3: Vector3<T>) -> Self {
		let half: T = nalgebra::convert(0.5);

		let v1 = p2 - p1;
		let v2 = p3 - p1;

		let normal = v1.cross(&v2).normalize();

		let mid12: Vector3<T> = (p1 + p2) * half;
		let mid13 = (p1 + p3) * half;

		let perp12 = normal.cross(&v1);
		let perp13 = normal.cross(&v2);

		let diff: Vector3<T> = mid13 - mid12;
		let denom = perp12.cross(&perp13).dot(&normal);

		let t = if denom.abs() > T::default_epsilon() {
			diff.cross(&perp13).dot(&normal) / denom
		} else {
			T::zero()
		};

		let center = mid12 + perp12 * t;

		let r_start: Vector3<T> = p1 - center;
		let r_end: Vector3<T> = p3 - center;
		let r_via: Vector3<T> = p2 - center;

		let x_axis = r_start.normalize();
		let y_axis = normal.cross(&x_axis);

		let angle_via = r_via.dot(&y_axis).atan2(r_via.dot(&x_axis));
		let angle_end = r_end.dot(&y_axis).atan2(r_end.dot(&x_axis));

		let tau: T = nalgebra::convert(TAU);
		let arc_angle = if angle_via >= T::zero() {
			if angle_end >= T::zero() && angle_end >= angle_via {
				angle_end
			} else {
				angle_end + tau
			}
		} else {
			if angle_end <= T::zero() && angle_end <= angle_via {
				angle_end
			} else {
				angle_end - tau
			}
		};

		Self {
			center,
			r_start,
			normal,
			arc_angle,
		}
	}

	/// Position on the arc at parameter s at [0, 1].
	fn point_at(&self, s: T) -> Vector3<T> {
		let angle = self.arc_angle * s;
		let x_axis = self.r_start.normalize();
		let y_axis = self.normal.cross(&x_axis);
		let radius = self.r_start.norm();

		self.center + x_axis * (angle.cos() * radius) + y_axis * (angle.sin() * radius)
	}
}
