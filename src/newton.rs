use crate::lambert_eqns::{f, f_prime};

pub fn newton(
    r1: f64,
    r2: f64,
    a: f64,
    dt: f64,
    x_0: f64,
    max_itrs: u32,
    tol: f64,
) -> Result<f64, String> {
    println!("-----Newton's Method Output-----");
    let mut x_i = x_0;

    for iteration in 1..max_itrs + 1 {
        let f_x = f(r1, r2, a, x_i, dt);
        let fprime = f_prime(r1, r2, a, x_i);

        let x_next = x_i - (f_x / fprime);

        if (x_next - x_i).abs() < tol {
            println!("-----Converged after {} iterations-----", iteration);
            return Ok(x_next);
        }

        println!(
            "Itr {}: z = {:.5}  |  f(z) = {:.2}  |   f'(z) = {:.2}",
            iteration, x_i, f_x, fprime
        );

        x_i = x_next;
    }

    Err("-----Exceeded iteration limit without convergence-----".to_string())
}

// -------------------------------
//             TESTS
// -------------------------------
