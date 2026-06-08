use nalgebra::{Isometry3, RealField};
use simba::scalar::SubsetOf;

use crate::{Error, kinematics::Chain};

pub mod jacobian;

pub trait IkSolver<const D: usize, const J: usize, T: RealField + SubsetOf<f64> + Copy> {
	fn solve(&self, chain: &mut Chain<D, J, T>, target: Isometry3<T>) -> Result<(), Error>;
}
