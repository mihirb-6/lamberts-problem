use crate::lambert_eqns::{MU, y};
use crate::stumpff::{stumpff_c, stumpff_s};

#[allow(unused)]
pub fn lagrange_f(r1: f64, r2: f64, a: f64, z_root: f64) -> f64 {
    1. - (y(r1, r2, a, z_root) / r1)
}

#[allow(unused)]
pub fn lagrange_g(r1: f64, r2: f64, a: f64, z_root: f64) -> f64 {
    a * (y(r1, r2, a, z_root) / MU).sqrt()
}

#[allow(unused)]
pub fn lagrange_fdot(r1: f64, r2: f64, a: f64, z_root: f64) -> f64 {
    (MU.sqrt() / (r1 * r2))
        * (y(r1, r2, a, z_root) / stumpff_c(z_root)).sqrt()
        * (z_root * stumpff_s(z_root) - 1.)
}

#[allow(unused)]
pub fn lagrange_gdot(r1: f64, r2: f64, a: f64, z_root: f64) -> f64 {
    1. - (y(r1, r2, a, z_root) / r2)
}
