#[allow(unused)]
pub fn cross_product(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    let product: [f64; 3] = [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ];

    product
}

#[allow(unused)]
pub fn dot_product(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let product = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

    product
}

#[allow(unused)]
pub fn magnitude(a: &[f64; 3]) -> f64 {
    (a[0].powi(2) + a[1].powi(2) + a[2].powi(2)).sqrt()
}
