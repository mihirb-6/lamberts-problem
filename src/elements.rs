use std::f64::consts::PI;

use nalgebra::Vector3;

pub struct Elements {
    pub angular_momentum: f64,      // kg * m^2 / s
    pub inclination: f64,           //rad
    pub raan: f64,                  // rad
    pub eccentricity: f64,          // dimensionless
    pub argument_of_periapsis: f64, // rad
    pub true_anomaly: f64,          // rad
}

#[allow(unused)]
impl Elements {
    pub fn period(&self, mu: f64) -> f64 {
        let a = self.seminajor_axis(mu);
        2. * PI / mu.sqrt() * a.powf(1.5)
    }
    pub fn seminajor_axis(&self, mu: f64) -> f64 {
        let r_p = self.periapsis(mu);
        let r_a = self.apoapsis(mu);
        0.5 * (r_p + r_a)
    }
    pub fn eccentric_anomaly(&self) -> f64 {
        let mag_e = self.eccentricity;
        let theta = self.true_anomaly;
        2. * (((1. - mag_e) / (1. + mag_e)).sqrt() * (theta / 2.).tan()).atan()
    }
    pub fn mean_anomaly(&self) -> f64 {
        let e1 = self.eccentric_anomaly();
        let mag_e = self.eccentricity;

        e1 - (mag_e * e1.sin())
    }
    pub fn time_since_perapsis(&self, mu: f64) -> f64 {
        let mag_h = self.angular_momentum;
        let mag_e = self.eccentricity;
        let me1 = self.mean_anomaly();

        (mag_h.powi(3) / mu.powi(2)) * 1. / (1. - mag_e.powi(2)).powf(1.5) * me1
    }
    pub fn periapsis(&self, mu: f64) -> f64 {
        let mag_h = self.angular_momentum;
        let mag_e = self.eccentricity;
        (mag_h.powi(2) / mu) * 1. / (1. + mag_e * 0_f64.cos())
    }
    pub fn apoapsis(&self, mu: f64) -> f64 {
        let mag_h = self.angular_momentum;
        let mag_e = self.eccentricity;
        (mag_h.powi(2) / mu) * 1. / (1. + mag_e * 180_f64.to_radians().cos())
    }
}

// ------- get_elements --------
// Inputs:
//         r: position vector at time t [m]
//         v: velocity vector at time t [m/s]
// Outputs:
//         (p): period [s]
//         (elements): 6 orbital elements [h, i, raan, e, w, theta]
//         (t_1): time since pariapsis
//         (r_p): Periapsis [m]
//         (r_a): Apoapsis [m]
pub fn get_elements(r: Vector3<f64>, v: Vector3<f64>, mu: f64) -> Elements {
    // Distance (r)
    let mag_r = r.magnitude();

    // Speed (v)
    let mag_v = v.magnitude();

    // Radial Velocity (v_r)
    let vr = r.dot(&v) / mag_r;

    /*
    match vr {
        vr if vr > 0. => println!("-> Object is flying away from periapsis"),
        vr if vr < 0. => println!("-> Object is flying towards periapsis"),
        _ => println!("vr = 0"),
    }
    */

    // Specific Angular Momentum (h):
    let h = r.cross(&v);

    // Magnitude of h                                       =>> 1st Element
    let mag_h = h.magnitude();

    // Inclination (i)                                      =>> 2nd Element
    let i = (h.z / mag_h).acos();

    // Node line Vector (N)
    let k: Vector3<f64> = Vector3::new(0., 0., 1.); // z-axis (K-hat) unit vector
    let n = k.cross(&h);

    // Magnitude of N
    let mag_n = n.magnitude();

    // Right Ascension of the Ascending Node (Omega) (RAAN) =>> 3rd Element
    let mut raan = (n.x / mag_n).acos();

    if n.y < 0. {
        raan = 360. * PI / 180. - (n.x / mag_n).acos();
    }

    // Eccentricity Vecotor (e)
    let e = (1. / mu) * ((mag_v.powi(2) - (mu / mag_r)) * &r - (mag_r * vr * &v));

    // Eccentricity                                         =>> 4th Element
    let mag_e = e.magnitude();

    // Argument of periapsis                                =>> 5th Element
    let mut w = (n.dot(&e) / (mag_n * mag_e)).acos();

    if e.z < 0. {
        w = 360. * PI / 180. - (n.dot(&e) / (mag_n * mag_e)).acos();
    }

    // True Anomaly:                                        =>>6th Element
    let mut theta = (e.dot(&r) / (mag_e * mag_r)).acos();

    if vr < 0. {
        theta = 360. * PI / 180. - (e.dot(&r) / (mag_e * mag_r)).acos();
    }

    let orbital_elements = Elements {
        angular_momentum: mag_h,
        inclination: i,
        raan: raan,
        eccentricity: mag_e,
        argument_of_periapsis: w,
        true_anomaly: theta,
    };
    // Return a tuple of a vector containing elements
    orbital_elements
}
