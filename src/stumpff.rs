#[allow(unused)]
fn stumpff_c(z: f64) -> f64 {
    match z {
        z if z > 0. => (1. - z.sqrt().cos()) / z,
        z if z == 0. => 1. / 2.,
        z if z < 0. => (1. - -z.sqrt().cosh()) / z,
        _ => panic!(),
    }
}

#[allow(unused)]
fn stumpff_s(z: f64) -> f64 {
    match z {
        z if z > 0. => (z.sqrt() - z.sqrt().sin()) / z.powi(3).sqrt(),
        z if z == 0. => 1. / 6.,
        z if z < 0. => (-z.sqrt().sinh() - -z.sqrt()) / -z.powi(3).sqrt(),
        _ => panic!(),
    }
}
