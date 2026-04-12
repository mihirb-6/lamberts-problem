use nalgebra::Vector3;

mod constants;
mod elements;
mod lagrange_coeffs;
mod lambert_eqns;
mod newton;
mod plot;
mod stumpff;
mod vectors;

use crate::constants::MU;
use crate::elements::get_elements;
use crate::lagrange_coeffs::{lagrange_f, lagrange_fdot, lagrange_g, lagrange_gdot};
use crate::newton::newton;
use crate::plot::plot_orbit;
use std::f64::consts::PI;

#[allow(unused)]
pub fn main() {
    let dt: f64 = 500. * 3600.; // [s]
    let r1 = Vector3::new(4673.0, 3623., -6654.0);
    let r2 = Vector3::new(-1230.0, -6825., -3493.);

    let (v1, v2) = lambert(r1, r2, Direction::Prograde, dt);
    let (a, elements, t_1, r_p, r_a) = get_elements(r1, v1);

    let h = elements.x;
    let i = elements.y;
    let raan = elements.z;
    let e = elements.w;
    let w = elements.a;
    let theta = elements.b;

    // Orbital Elements Print Statement
    println!("-----Orbital Elements------");
    println!("Semi-major axis: a = {:.4} [km]", a);
    println!("Periapsis: r_p = {:.4} [km]", r_p);
    println!("Apoapsis: r_a = {:.4} [km]", r_a);
    println!("Eccentricity: e = {:.4}", e);
    println!("Specific Angular Momentum: h = {:.2} [km^2 s^-1]", h);
    println!("Inclination: i = {:.2}°", i.to_degrees());
    println!("RA of Ascending Node (Ω): = {:.2}°", raan.to_degrees());
    println!("Argument of Periapsis (ω) = {:.2}°", w.to_degrees());
    println!("True Anomaly (θ) = {:.2}°", theta.to_degrees());
    println!("----------------------------");
    println!(
        "Perigee encounter in {:.1} s = {:.2} min = {:.2} hr(s)",
        t_1,
        t_1 / 60.,
        t_1 / 3600.
    );

    plot_orbit(e, h, i, raan, w, MU);
}

#[allow(unused)]
pub enum Direction {
    Prograde,
    Retrograde,
}

fn lambert(
    r1_vector: Vector3<f64>,
    r2_vector: Vector3<f64>,
    direction: Direction,
    dt: f64,
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

    let z_initial = 0.; //MU.sqrt() * dt / lambert_a;
    //println!("z0 = {}", z_initial);

    let max_iterations = 500;
    let tolerance = 1e-9;

    let z_root = newton(r1, r2, lambert_a, dt, z_initial, max_iterations, tolerance)
        .unwrap_or_else(|e| panic!("{}", e));

    println!("-> Root of z: {:.5}", z_root);

    // Calculate lagrange coefficients
    let f = lagrange_f(r1, r2, lambert_a, z_root);
    #[allow(unused)]
    let fdot = lagrange_fdot(r1, r2, lambert_a, z_root);
    let g = lagrange_g(r1, r2, lambert_a, z_root);
    let gdot = lagrange_gdot(r1, r2, lambert_a, z_root);

    // Lagrange Coeff. Print Statement
    /*
    println!(
        "f = {:.4}\nfdot = {:.4}\ng = {:.4} s\ngdot = {:.4}",
        f, fdot, g, gdot
    );
    */

    let r1v = Vector3::from(r1_vector);
    let r2v = Vector3::from(r2_vector);

    // Compute v1 and v2
    let v1 = 1. / g * (r2v - f * r1v);
    let v2 = 1. / g * (gdot * r2v - r1v);

    // v1, v2 Print Statement
    //println!("v1 = {:.4}\nv2 = {:.4}", v1, v2);

    (v1, v2)
}
