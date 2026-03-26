mod lagrange_coeffs;
mod lambert_eqns;
mod newton;
mod stumpff;
mod vectors;

//use crate::lagrange_coeffs::{lagrange_f, lagrange_fdot, lagrange_g, lagrange_gdot};
use crate::lambert_eqns::MU;
use crate::newton::newton;
use crate::vectors::{cross_product, dot_product, magnitude};

fn main() {
    let r1: [f64; 3] = [5000., 10000., 2100.]; // [km]
    let r2: [f64; 3] = [-14600., 2500., 7000.]; // [km]
    let dt: f64 = 5.0 * 3600.; // [s]\
    lambert(r1, r2, Direction::Prograde, dt);
}

#[allow(unused)]
enum Direction {
    Prograde,
    Retrograde,
}

#[allow(unused)]
fn lambert(r1_vector: [f64; 3], r2_vector: [f64; 3], direction: Direction, dt: f64) {
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
    //
    let lambert_a = dtheta.sin() * ((r1 * r2) / (1. - dtheta.cos())).sqrt();

    let z_initial = MU.sqrt() * dt / lambert_a;
    println!("z0 = {}", z_initial);

    let max_iterations = 50;
    let tolerance = 1e-3;

    let z_root = newton(r1, r2, lambert_a, dt, z_initial, max_iterations, tolerance)
        .unwrap_or_else(|e| panic!("{}", e));

    println!("root of z: {}", z_root);

    // Calculate lagrange coefficients
    //let f = lagrange_f(r1, r2, lambert_a, z_root);
    //let fdot = lagrange_fdot(r1, r2, lambert_a, z_root);
    //let g = lagrange_g(r1, r2, lambert_a, z_root);
    //let gdot = lagrange_gdot(r1, r2, lambert_a, z_root);

    // Compute v1 and v2
    // v1 = 1/g * (r2_vector - f*r1_vector)
    // v2 = 1/g * (gdot*r2_vevtor - r1_vector)
}
