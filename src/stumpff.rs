#[allow(unused)]
pub fn stumpff_c(x: f64) -> f64 {
    match x {
        x if x > 0. => (1. - x.sqrt().cos()) / x,
        x if x == 0. => 1. / 2.,
        x if x < 0. => (1. - -x.sqrt().cosh()) / x,
        _ => panic!(),
    }
}

#[allow(unused)]
pub fn stumpff_s(x: f64) -> f64 {
    match x {
        x if x > 0. => (x.sqrt() - x.sqrt().sin()) / x.powi(3).sqrt(),
        x if x == 0. => 1. / 6.,
        x if x < 0. => (-x.sqrt().sinh() - -x.sqrt()) / -x.powi(3).sqrt(),
        _ => panic!(),
    }
}
