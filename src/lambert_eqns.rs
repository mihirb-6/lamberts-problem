use crate::stumpff::{stumpff_c, stumpff_s};
use std::f64::{NAN, consts::SQRT_2};

/* Characteristic equations taken from the Curtis orbital mechanics
 * textbook. Functions f and fprime are used in newton.rs (an iterative root finder)
 * Function y is called when f or fprime is called newton.rs
 */
pub fn y(r1: f64, r2: f64, a: f64, x: f64) -> Result<f64, String> {
    if x != NAN {
        return Ok(r1 + r2 + (a * (x * stumpff_s(x) - 1.) / stumpff_c(x).sqrt()));
    }
    Err("Bad value passed into y(z)".to_string())
}

pub fn f(r1: f64, r2: f64, a: f64, x: f64, dt: f64, mu: f64) -> f64 {
    let y_uw = y(r1, r2, a, x).unwrap(); // unwrap b/c fn y returns a Result

    ((y_uw / stumpff_c(x)).powf(1.5) * stumpff_s(x)) + (a * y_uw.sqrt()) - (mu.sqrt() * dt)
}

pub fn f_prime(r1: f64, r2: f64, a: f64, x: f64) -> f64 {
    let y_uw = y(r1, r2, a, x).unwrap(); // unwrap b/c fn y returns a Result

    if x == 0. {
        (SQRT_2 / 40.) * y_uw.powf(1.5)
            + (a / 8.) * (y_uw.sqrt() + (a * (0.5 * (1. / y_uw)).sqrt()))
    } else {
        let s_div_c = stumpff_s(x) / stumpff_c(x); // easier variables to work with
        let y_div_c = y_uw / stumpff_c(x);
        let c_div_y = stumpff_c(x) / y_uw;
        y_div_c.powf(1.5)
            * ((0.5 / x) * (stumpff_c(x) - (1.5 * s_div_c)) + 0.75 * s_div_c * stumpff_s(x))
            + a * 0.125 * (3. * s_div_c * y_uw.sqrt() + (a * c_div_y.sqrt()))
    }
}
