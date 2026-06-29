use core::f32::consts::FRAC_PI_4;

use kinema::SVector;

use kinema_macro::load_urdf;

fn main() {
	let mut robot = load_urdf!("../robot.urdf");

	robot.set_joint_positions_clamped(SVector::from([0.0, FRAC_PI_4, 0.0, 0.0, 0.0, 0.0]));
	robot.update_transforms();

	println!("TCP position: {}", robot.end_transform());
}
