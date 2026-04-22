use crate::lambert_eqns::{f, f_prime};

pub fn newton(
    r1: f64,
    r2: f64,
    a: f64,
    dt: f64,
    x_0: f64,
    max_itrs: u32,
    tol: f64,
    mu: f64,
) -> Result<f64, String> {
    println!("------------NEWTON'S METHOD OUTPUT------------");
    let mut x_i = x_0;

    for iteration in 1..max_itrs + 1 {
        let f_x = f(r1, r2, a, x_i, dt, mu);
        let fprime = f_prime(r1, r2, a, x_i);

        let x_next = x_i - (f_x / fprime);

        if (x_next - x_i).abs() < tol {
            println!("-> Root of z: {:.5}", x_next);
            println!("---------CONVERGED AFTER {} ITERATIONS---------", iteration);
            return Ok(x_next);
        }

        println!(
            "Itr {}: z = {:.7}  |  f(z) = {:.2}  |   f'(z) = {:.2}",
            iteration, x_i, f_x, fprime
        );

        x_i = x_next;
    }

    Err("-----EXCEEDED ITERATION LIMIT WITHOUT CONVERGENCE-----".to_string())
}
