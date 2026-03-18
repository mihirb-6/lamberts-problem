use std::f64::consts::SQRT_2;

use crate::stumpff::{stumpff_c, stumpff_s};
use crate::vectors::{cross_product, dot_product, magnitude};

mod stumpff;
mod vectors;

pub const MU: f64 = 3.986004418e5; // [km^3 s^-2]

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

#[allow(unused)]
fn y(v1: f64, v2: f64, a: f64, x: f64) -> f64 {
    v1 + v2 + (a * (x * stumpff_s(x) - 1.) / stumpff_c(x))
}

#[allow(unused)]
fn f(v1: f64, v2: f64, a: f64, x: f64, dt: f64) -> f64 {
    ((y(v1, v2, a, x) / stumpff_c(x)).powf(1.5) * stumpff_s(x)) + (a * y(v1, v2, a, x).sqrt())
        - MU.sqrt() * dt
}

#[allow(unused)]
fn f_prime(v1: f64, v2: f64, a: f64, x: f64) -> f64 {
    if x == 0. {
        SQRT_2 / 40. * y(v1, v2, a, 0.).powf(1.5)
            + (a / 8.) * (y(v1, v2, a, 0.).sqrt() + (a * ((1. / 2.) * y(v1, v2, a, 0.))))
    } else {
        let s_div_c = stumpff_s(x) / stumpff_c(x);
        let y_div_c = y(v1, v2, a, x) / stumpff_c(x);
        let c_div_y = stumpff_c(x) / y(v1, v2, a, x);

        y_div_c.powf(1.5)
            * ((0.5 / x) * (stumpff_c(x) - (1.5 * s_div_c)) + 0.75 * s_div_c * stumpff_s(x))
            + a / 0.125 * (3. * s_div_c * y(v1, v2, a, x).sqrt() + (a * c_div_y.sqrt()))
    }
}

fn newtown(f: fn(f64) -> f64, x_0: f64, max_itrs: u32, tol: f64) {
    let x_old = x_0;
    let x_n = x_0;

    for iteration in 1..max_itrs + 1 {
        let f_x = f(x_n);
    }
}
