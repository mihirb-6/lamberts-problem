use std::f64::consts::PI;

use nalgebra::{Vector3, Vector6};

use crate::constants::MU;

#[allow(unused)]
pub fn get_elements(r: Vector3<f64>, v: Vector3<f64>) -> Vector6<f64> {
    // Distance (r)
    let mag_r = r.magnitude();

    // Speed (v)
    let mag_v = v.magnitude();

    // Radial Velocity (v_r)
    let vr = r.dot(&v) / mag_r;

    match vr {
        vr if vr > 0. => println!("Satellite is flying away from periapsis"),
        vr if vr < 0. => println!("Satellite is flying towards periapsis"),
        _ => println!("vr = 0"),
    }

    // Specific Angular Momentum (h):
    let h = r.cross(&v);

    // Magnitude of h: 1st Element
    let mag_h = h.magnitude();

    // Inclination (i): 2nd Element
    let i = (h.z / mag_h).acos();

    // Node line Vector (N)
    let k: Vector3<f64> = Vector3::new(0., 0., 1.); // z-axis (K-hat) unit vector
    let n = k.cross(&h);

    // Magnitude of N
    let mag_n = n.magnitude();

    // Right Ascension of the Ascending Node (Omega) (RAAN): 3rd Element
    let mut raan = (n.x / mag_n).acos();

    if n.y < 0. {
        raan = 360. * PI / 180. - (n.x / mag_n).acos();
    }

    // Eccentricity Vecotor (e)
    let e = (1. / MU) * ((mag_v.powi(2) - (MU / mag_r)) * &r - (mag_r * vr * &v));

    // Eccentricity: 4th Element
    let mag_e = e.magnitude();

    // Argument of periapsis: 5th Element
    let mut w = (n.dot(&e) / (mag_n * mag_e)).acos();

    if e.z < 0. {
        w = 360. * PI / 180. - (n.dot(&e) / (mag_n * mag_e)).acos();
    }

    // True Anomaly: 6th Element
    let mut theta = (e.dot(&r) / (mag_e * mag_r)).acos();

    if vr < 0. {
        theta = 360. * PI / 180. - (e.dot(&r) / (mag_e * mag_r)).acos();
    }

    Vector6::new(mag_h, i, raan, mag_e, w, theta)
}
