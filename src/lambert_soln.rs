use crate::lagrange_coeffs::{lagrange_f, lagrange_fdot, lagrange_g, lagrange_gdot};
use crate::newton::newton;
use nalgebra::Vector3;
use std::f64::consts::PI;

#[allow(unused)]
pub enum Direction {
    Prograde,
    Retrograde,
}

/* Bringing the project altogether, combining multiple functions
 * to generate v1 and v2 from r1, r2, dt, and mu
 */
pub fn lambert(
    r1_vector: Vector3<f64>,
    r2_vector: Vector3<f64>,
    direction: Direction,
    dt: f64,
    mu: f64,
    z_init: f64,
    max_itrs: u32,
) -> (Vector3<f64>, Vector3<f64>) {
    let r1 = r1_vector.magnitude();
    let r2 = r2_vector.magnitude();
    //println!("r1 = {}, r2 = {}", r1, r2);

    // Compute dot product of r1 and r2
    let r1_dot_r2 = r1_vector.dot(&r2_vector);
    //println!("r1 dot r2 = {}", r1_dot_r2);

    // Compute the cross product of r1 and r2
    let r1_cross_r2 = r1_vector.cross(&r2_vector);

    // Compute dot product of the cross product of r1 and r2
    //let mag_of_r1_r2_cross = magnitude(&r1_cross_r2);

    let mut dtheta: f64 = (r1_dot_r2 / (r1 * r2)).acos();

    // Howard Curtis Eqn. (5.26) pg. 203
    match direction {
        Direction::Prograde => {
            if r1_cross_r2.z >= 0. {
                dtheta = (r1_dot_r2 / (r1 * r2)).acos();
            } else {
                dtheta = 2. * PI - dtheta;
            }
        }
        Direction::Retrograde => {
            if r1_cross_r2.z < 0. {
                dtheta = (r1_dot_r2 / (r1 * r2)).acos();
            } else {
                dtheta = 2. * PI - dtheta;
            }
        }
    }

    // Handling dtheta = 0 case:
    // Introducting a slight perturbation
    let eps = 1e-7;

    if dtheta.abs() < eps {
        dtheta = eps;
    }

    if (dtheta - PI).abs() < eps {
        dtheta = PI - eps;
    }

    //println!("dtheta = {:.2 } rad = {:.2} deg",dtheta,dtheta.to_degrees());

    // Compute Lambert Parameter A
    // Howard Curtis Eqn. (5.35) pg. 205
    let lambert_a = dtheta.sin() * ((r1 * r2) / (1. - dtheta.cos())).sqrt();
    //println!("A = {:.2}", lambert_a);

    // Provide an inital estimate for z
    // can also use a plot of f and f prime (lambert_eqns.rs) to estimate
    let z_initial: f64 = z_init; //MU.sqrt() * dt / lambert_a; <- provided by curtis, might implement
    //println!("z0 = {}", z_initial);

    // Set max times to iterate newton
    let max_iterations: u32 = max_itrs;

    // Set a tolerance
    // i.e how close should z_n and z_n+1 be before halting iteration
    // and determining z_n+1 as a root
    let tolerance: f64 = 1e-9;

    // Newton's method to find a root iteratively (see newton.rs)
    let z_root = newton(
        r1,
        r2,
        lambert_a,
        dt,
        z_initial,
        max_iterations,
        tolerance,
        mu,
    )
    .unwrap_or_else(|e| panic!("{}", e));

    // Calculate lagrange coefficients (see langrange_coeffs.rs)
    let f = lagrange_f(r1, r2, lambert_a, z_root);
    #[allow(unused)]
    let fdot = lagrange_fdot(r1, r2, lambert_a, z_root, mu);
    let g = lagrange_g(r1, r2, lambert_a, z_root, mu);
    let gdot = lagrange_gdot(r1, r2, lambert_a, z_root);

    // Lagrange Coeff. Print Statement
    /*
    println!(
        "f = {:.4}\nfdot = {:.4}\ng = {:.4} s\ngdot = {:.4}",
        f, fdot, g, gdot
    );
    */

    let r1v: Vector3<f64> = Vector3::from(r1_vector);
    let r2v: Vector3<f64> = Vector3::from(r2_vector);

    // Compute v1 and v2
    let v1 = 1. / g * (r2v - f * r1v);
    let v2 = 1. / g * (gdot * r2v - r1v);

    // v1, v2 Print Statement
    //println!("v1 = {:.4}\nv2 = {:.4}", v1, v2);

    (v1, v2)
}
