use nalgebra::Vector3;

mod constants;
mod elements;
mod lagrange_coeffs;
mod lambert_eqns;
mod newton;
mod stumpff;
mod vectors;

use crate::constants::MU;
use crate::elements::get_elements;
use crate::lagrange_coeffs::{lagrange_f, lagrange_fdot, lagrange_g, lagrange_gdot};
use crate::newton::newton;
use std::f64::consts::PI;

#[allow(unused)]
pub fn main() {
    let t = 2.0 * std::f64::consts::PI * (7000.0_f64.powi(3) / MU).sqrt();
    let dt = t / 4.0;
    //let dt: f64 = 1.0 * 3600.; // [s]
    let r1 = Vector3::new(8000.0, 0.0, 0.0);
    let r2 = Vector3::new(12000.0, 0.0, 0.0);

    let (v1, v2) = lambert(r1, r2, Direction::Prograde, dt);
    let elements = get_elements(r1, v1);

    let h = elements.x;
    let i = elements.y.to_degrees();
    let raan = elements.z.to_degrees();
    let e = elements.w;
    let w = elements.a.to_degrees();
    let theta = elements.b.to_degrees();

    // Perigee and Apogee Radii
    let r_p = (h.powi(2) / MU) * 1. / (1. + e * 0_f64.cos());
    let r_a = (h.powi(2) / MU) * 1. / (1. + e * 180_f64.to_radians().cos());

    // Semimajor Axis
    let a = 0.5 * (r_p + r_a);

    // Period
    let period = 2. * PI / MU.sqrt() * a.powf(1.5);

    // Eccentric Anomaly
    let e1 = 2. * (((1. - e) / (1. + e)).sqrt() * (theta.to_radians() / 2.).tan()).atan();

    // Mean Anomaly
    let me1 = e1 - (e * e1.sin());

    // Time since periapsis
    let t_1 = (h.powi(3) / MU.powi(2)) * 1. / (1. - e.powi(2)).powf(1.5) * me1;

    // Orbital Elements Print Statement
    println!(
        "
        a = {a:.4} [km]\n
        e = {e:.4}\n
        h = {h:.2} [km^2 s^-1]\n
        i = {i:.2}°\n
        RAAN = {raan:.2}°\n
        w = {w:.2}°\n
        theta = {theta:.2}°\n"
    );

    println!("Perigee encounter in {t_1:.1} s");
}

#[allow(unused)]
pub enum Direction {
    Prograde,
    Retrograde,
}

pub fn lambert(
    r1_vector: Vector3<f64>,
    r2_vector: Vector3<f64>,
    direction: Direction,
    dt: f64,
) -> (Vector3<f64>, Vector3<f64>) {
    let r1 = r1_vector.magnitude();
    let r2 = r2_vector.magnitude();
    println!("r1 = {}, r2 = {}", r1, r2);

    // Compute dot product of r1 and r2
    let r1_dot_r2 = r1_vector.dot(&r2_vector);
    println!("r1 dot r2 = {}", r1_dot_r2);

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

    println!(
        "dtheta = {:.2 } rad = {:.2} deg",
        dtheta,
        dtheta.to_degrees()
    );

    // Compute Lambert Parameter A
    // Howard Curtis Eqn. (5.35) pg. 205
    let lambert_a = dtheta.sin() * ((r1 * r2) / (1. - dtheta.cos())).sqrt();
    println!("A = {:.2}", lambert_a);

    let z_initial = 1.5; //MU.sqrt() * dt / lambert_a;
    println!("z0 = {}", z_initial);

    let max_iterations = 100;
    let tolerance = 1e-6;

    let z_root = newton(r1, r2, lambert_a, dt, z_initial, max_iterations, tolerance)
        .unwrap_or_else(|e| panic!("{}", e));

    println!("Root of z: {:.5}", z_root);

    // Calculate lagrange coefficients
    let f = lagrange_f(r1, r2, lambert_a, z_root);
    #[allow(unused)]
    let fdot = lagrange_fdot(r1, r2, lambert_a, z_root);
    let g = lagrange_g(r1, r2, lambert_a, z_root);
    let gdot = lagrange_gdot(r1, r2, lambert_a, z_root);

    // Lagrange Coeff. Print Statement
    println!(
        "f = {:.4}\nfdot = {:.4}\ng = {:.4} s\ngdot = {:.4}",
        f, fdot, g, gdot
    );

    let r1v = Vector3::from(r1_vector);
    let r2v = Vector3::from(r2_vector);

    // Compute v1 and v2
    let v1 = 1. / g * (r2v - f * r1v);
    let v2 = 1. / g * (gdot * r2v - r1v);

    // v1, v2 Print Statement
    println!("v1 = {:.4}\nv2 = {:.4}", v1, v2);

    (v1, v2)
}

#[cfg(test)]
mod main_tests {
    use nalgebra::Vector3;
    use std::f64::consts::PI;

    use crate::{Direction, lambert};

    const MU: f64 = 398600.4418;

    fn propagate(r: Vector3<f64>, v: Vector3<f64>, dt: f64) -> Vector3<f64> {
        let steps = 10_000;
        let h = dt / steps as f64;

        let mut r_curr = r;
        let mut v_curr = v;

        for _ in 0..steps {
            let r_norm = r_curr.norm();
            let acc = -MU / r_norm.powi(3) * r_curr;

            v_curr += acc * h;
            r_curr += v_curr * h;
        }

        r_curr
    }

    #[test]
    fn prograde_quarter_orbit() {
        let r1 = Vector3::new(7000.0, 0.0, 0.0);
        let r2 = Vector3::new(0.0, 7000.0, 0.0);

        let t = 2.0 * PI * (7000.0_f64.powi(3) / MU).sqrt();
        let dt = t / 4.0;

        let (v1, _) = lambert(r1, r2, Direction::Prograde, dt);

        let r2_calc = propagate(r1, v1, dt);

        assert!((r2_calc - r2).norm() < 1e-1);

        // Direction check: prograde => +z angular momentum
        let h = r1.cross(&v1);
        assert!(h.z > 0.0);
    }

    #[test]
    fn retrograde_direction_check() {
        let r1 = Vector3::new(7000.0, 0.0, 0.0);
        let r2 = Vector3::new(0.0, 7000.0, 0.0);

        let t = 2.0 * PI * (7000.0_f64.powi(3) / MU).sqrt();
        let dt = t / 4.0;

        let (v1, _) = lambert(r1, r2, Direction::Retrograde, dt);

        let h = r1.cross(&v1);

        // Retrograde should flip angular momentum
        assert!(h.z < 0.0);
    }

    #[test]
    fn time_scaling_behavior() {
        let r1 = Vector3::new(5000.0, 10000.0, 2100.0);
        let r2 = Vector3::new(-14600.0, 2500.0, 7000.0);

        let (v_fast, _) = lambert(r1, r2, Direction::Prograde, 1000.0);
        let (v_slow, _) = lambert(r1, r2, Direction::Prograde, 10000.0);

        assert!(v_fast.norm() > v_slow.norm());
    }

    #[test]
    fn endpoint_accuracy_general_case() {
        let r1 = Vector3::new(5000.0, 10000.0, 2100.0);
        let r2 = Vector3::new(-14600.0, 2500.0, 7000.0);
        let dt = 3600.0;

        let (v1, _) = lambert(r1, r2, Direction::Prograde, dt);

        let r2_calc = propagate(r1, v1, dt);

        assert!((r2_calc - r2).norm() < 1.0, "Propagation error too large");
    }

    #[test]
    fn collinear_edge_case() {
        let r1 = Vector3::new(8000.0, 0.0, 0.0);
        let r2 = Vector3::new(12000.0, 0.0, 0.0);

        let dt = 2000.0;

        let (v1, _) = lambert(r1, r2, Direction::Prograde, dt);

        assert!(v1.iter().all(|x| x.is_finite()));
        assert!(v1.y.abs() < 1e-6 && v1.z.abs() < 1e-6);
    }

    #[test]
    fn energy_consistency() {
        let r1 = Vector3::new(7000.0, 0.0, 0.0);
        let r2 = Vector3::new(8000.0, 1000.0, 0.0);
        let dt = 3000.0;

        let (v1, _) = lambert(r1, r2, Direction::Prograde, dt);

        let energy = v1.norm_squared() / 2.0 - MU / r1.norm();

        assert!(energy.is_finite());
    }
}
