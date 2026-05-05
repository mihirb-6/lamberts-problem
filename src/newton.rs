use crate::lambert_eqns::{f, f_prime};

/* Iterative method to find a root
 * specifically for f and fprime from lambert_eqns.rs,
 * but can generalize if necessary
 */
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

    // for loop, start iterating until our defined max iterations
    for iteration in 1..max_itrs + 1 {
        let f_x = f(r1, r2, a, x_i, dt, mu); // easier variable to work with
        let fprime = f_prime(r1, r2, a, x_i); // same here

        // estimate the next root
        let x_next = x_i - (f_x / fprime); // characteristic equation - from Curtis

        // if z_n+1 - z_n is within our defined tolerance
        if (x_next - x_i).abs() < tol {
            println!("FOUND ROOT OF z: {:.5}", x_next);
            println!("---------CONVERGED AFTER {} ITERATIONS---------", iteration);
            // return z_n+1
            return Ok(x_next);
        }

        if (x_next) < 0. {
            panic!("******HYPERBOLIC ORBIT, CHOOSE DIFFERENT VECTORS/TIMEFRAME/BODY******");
        }

        println!("Itr {}: z = {:.7}", iteration, x_i);

        // if z_n+1 - z_n is not within specified tolerance
        // z_n+1 becomes z_n and the process repeats until specified # of max iterations
        x_i = x_next;
    }

    // TODO! -> HYPERBOLIC ORBITS (e>1) WILL EXCEED ITERATION!!!!!!!!!!!!! <- RAISE ERROR
    // if went above max iteration threshold, raise an error
    Err("-----EXCEEDED ITERATION LIMIT WITHOUT CONVERGENCE-----".to_string())
}
