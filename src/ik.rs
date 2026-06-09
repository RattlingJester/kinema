use nalgebra::{Isometry3, RealField};
use simba::scalar::SubsetOf;

use crate::{Error, kinematics::Chain};

pub mod jacobian;

/// A trait which should be implemented for any custom solver
pub trait IkSolver<const D: usize, const J: usize, T: RealField + SubsetOf<f64> + Copy> {
	/// Convention is to set joint angles at the end of the solve, not return them
	fn solve(&self, chain: &mut Chain<D, J, T>, target: Isometry3<T>) -> Result<(), Error>;
}
