use crate::stumpff::{stumpff_c, stumpff_s};
use std::f64::consts::SQRT_2;
pub static MU: f64 = 3.986004418e5; // [km^3 s^-2]

#[allow(unused)]
pub fn y(v1: f64, v2: f64, a: f64, x: f64) -> f64 {
    v1 + v2 + (a * (x * stumpff_s(x) - 1.) / stumpff_c(x))
}

#[allow(unused)]
pub fn f(v1: f64, v2: f64, a: f64, x: f64, dt: f64) -> f64 {
    ((y(v1, v2, a, x) / stumpff_c(x)).powf(1.5) * stumpff_s(x)) + (a * y(v1, v2, a, x).sqrt())
        - MU.sqrt() * dt
}

#[allow(unused)]
pub fn f_prime(v1: f64, v2: f64, a: f64, x: f64) -> f64 {
    let s_div_c = stumpff_s(x) / stumpff_c(x);
    let y_div_c = y(v1, v2, a, x) / stumpff_c(x);
    let c_div_y = stumpff_c(x) / y(v1, v2, a, x);

    if x == 0. {
        (SQRT_2 / 40.) * y(v1, v2, a, 0.).powf(1.5)
            + (a / 8.) * (y(v1, v2, a, 0.).sqrt() + (a * (0.5 * (1. / y(v1, v2, a, 0.))).sqrt()))
    } else {
        y_div_c.powf(1.5)
            * ((0.5 / x) * (stumpff_c(x) - (1.5 * s_div_c)) + 0.75 * s_div_c * stumpff_s(x))
            + a * 0.125 * (3. * s_div_c * y(v1, v2, a, x).sqrt() + (a * c_div_y.sqrt()))
    }
}
