use crate::vectors::{cross_product, dot_product, magnitude};

mod lambert_eqns;
mod newton;
mod stumpff;
mod vectors;

fn main() {
    let _dt: f64 = 5.0 * 3600.; // [s]
    let _r1: [f64; 3] = [1., 2., 3.];
    let _r2: [f64; 3] = [4., 5., 6.];
    let _z_guess = 1.;
}

#[allow(unused)]
enum Direction {
    Prograde,
    Retrograde,
}

#[allow(unused)]
fn lambert(tof: f64, r1_vector: [f64; 3], r2_vector: [f64; 3], direction: Direction, z_guess: f64) {
    let r1 = magnitude(&r1_vector);
    let r2 = magnitude(&r2_vector);

    // Compute dot product of r1 and r2
    let r1_dot_r2 = dot_product(&r1_vector, &r2_vector);

    // Compute the cross product of r1 and r2
    let r1_cross_r2 = cross_product(&r1_vector, &r2_vector);

    // Compute dot product of the cross product of r1 and r2
    let mag_of_r1_r2_cross = magnitude(&r1_cross_r2);

    let dtheta: f64;

    // Howard Curtis Eqn. (5.26) pg. 203
    match direction {
        Direction::Prograde => {
            if r1_cross_r2[2] >= 0. {
                dtheta = (r1_dot_r2 / (r1 * r2)).acos();
            } else {
                dtheta = 360. - (r1_dot_r2 / (r1 * r2)).acos();
            }
        }
        Direction::Retrograde => {
            if r1_cross_r2[2] < 0. {
                dtheta = (r1_dot_r2 / (r1 * r2)).acos();
            } else {
                dtheta = 360. - (r1_dot_r2 / (r1 * r2)).acos();
            }
        }
    }

    // Compute Lambert Parameter A
    // Howard Curtis Eqn. (5.35) pg. 205
    let lambert_param = dtheta.sin() * ((r1 * r2) / (1. - dtheta.cos())).sqrt();

    let z = 0;
}
