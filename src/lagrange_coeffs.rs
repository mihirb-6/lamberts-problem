use crate::lambert_eqns::y;
use crate::stumpff::{stumpff_c, stumpff_s};

/* f, g, fdot*, gdot are used to generate v1 and v2
 * after a root (see newton.rs, lambert_eqns.rs) is found
 * (fdot is not necessary for the solution)
 */

pub fn lagrange_f(r1: f64, r2: f64, a: f64, z_root: f64) -> f64 {
    1. - (y(r1, r2, a, z_root).unwrap() / r1)
}

pub fn lagrange_g(r1: f64, r2: f64, a: f64, z_root: f64, mu: f64) -> f64 {
    a * (y(r1, r2, a, z_root).unwrap() / mu).sqrt()
}

#[allow(unused)]
pub fn lagrange_fdot(r1: f64, r2: f64, a: f64, z_root: f64, mu: f64) -> f64 {
    (mu.sqrt() / (r1 * r2))
        * (y(r1, r2, a, z_root).unwrap() / stumpff_c(z_root)).sqrt()
        * (z_root * stumpff_s(z_root) - 1.)
}

pub fn lagrange_gdot(r1: f64, r2: f64, a: f64, z_root: f64) -> f64 {
    1. - (y(r1, r2, a, z_root).unwrap() / r2)
}
