#[allow(unused)]
pub fn newton(
    v1: f64,
    v2: f64,
    a: f64,
    dt: f64,
    f: fn(v1: f64, v2: f64, a: f64, x: f64, dt: f64) -> f64,
    f_prime: fn(v1: f64, v2: f64, a: f64, x: f64) -> f64,
    x_0: f64,
    max_itrs: u32,
    tol: f64,
) -> Result<f64, String> {
    //const H: f64 = 1e-8; // small step for numerical differentiation

    let mut x_old = x_0;
    let mut x_n = x_0;

    for iteration in 1..max_itrs + 1 {
        //let f_x = f(x_n);
        //let fprime = (f(x_n + H) - f(x_n)) / H;
        let f_x = f(v1, v2, a, x_n, dt);
        let fprime = f_prime(v1, v2, a, x_n);

        x_old = x_n;
        x_n -= f_x / fprime;

        println!("Iteration {}: x = {}, f(x) = {}", iteration, x_n, f_x);

        if (x_n - x_old).abs() < tol {
            println!("Converged after {} iterations", iteration);
            return Ok(x_n);
        }
    }

    return Err("Exceeded iteration limit without convergence.".to_string());
}

// -------------------------------
//             TESTS
// -------------------------------

/*
#[cfg(test)]
mod newton_tests {
    use super::*;

    // -----------------------------------------------------------------------
    // Helpers
    // -----------------------------------------------------------------------

    /// Absolute-error float comparison.
    fn approx_eq(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    // Plain `fn` pointers required by the signature
    // (non-capturing closures coerce automatically).
    fn f_quadratic(x: f64) -> f64 {
        x * x - 2.0
    } // root: √2
    fn f_linear(x: f64) -> f64 {
        x - 3.0
    } // root: 3.0
    fn f_cubic(x: f64) -> f64 {
        x * x * x - x - 2.0
    } // root: ≈ 1.5214
    fn f_trig(x: f64) -> f64 {
        x.cos() - x
    } // root: ≈ 0.7391 (Dottie number)
    fn f_no_real_root(x: f64) -> f64 {
        x * x + 1.0
    } // no real root → diverges

    // -----------------------------------------------------------------------
    // BUG-EXPOSURE TESTS  (pass today, document the premature-return defect)
    // -----------------------------------------------------------------------

    /// The convergence guard fires before any update on iteration 1,
    /// so the function always returns x_start unchanged.
    #[test]
    #[ignore = "bug fixed"]
    fn bug_returns_x_start_immediately() {
        let result = newton(f_quadratic, 1.5, 100, 1e-6);
        // Should be √2 ≈ 1.4142, but the bug returns 1.5.
        assert_eq!(
            result, 1.5,
            "Bug: newton returned x_start without iterating"
        );
    }

    #[test]
    fn bug_zero_iterations_returns_zero() {
        // max_itrs = 0 → loop body never runs → falls through to the 0.0 sentinel.
        let result = newton(f_linear, 5.0, 0, 1e-6);
        assert_eq!(result, 0.0);
    }

    // -----------------------------------------------------------------------
    // INTENDED-BEHAVIOUR TESTS
    // -----------------------------------------------------------------------

    /// f(x) = x² – 2, root at √2. Classic Newton showcase.
    #[test]
    fn test_sqrt_2_quadratic() {
        let root = newton(f_quadratic, 1.5, 100, 1e-9);
        assert!(
            approx_eq(root, 2f64.sqrt(), 1e-6),
            "Expected √2 ≈ {:.9}, got {:.9}",
            2f64.sqrt(),
            root
        );
    }

    /// f(x) = x – 3 should converge in a single Newton step.
    #[test]
    fn test_linear_single_step() {
        let root = newton(f_linear, 0.0, 10, 1e-9);
        assert!(approx_eq(root, 3.0, 1e-9), "Expected 3.0, got {root}");
    }

    /// f(x) = x³ – x – 2, real root ≈ 1.5214.
    #[test]
    fn test_cubic_root() {
        let expected = 1.5213797068045676_f64;
        let root = newton(f_cubic, 2.0, 100, 1e-9);
        assert!(
            approx_eq(root, expected, 1e-6),
            "Expected ≈ {expected:.9}, got {root:.9}"
        );
    }

    /// Dottie number: fixed point of cos, f(x) = cos(x) – x ≈ 0.7391.
    #[test]
    fn test_dottie_trig() {
        let expected = 0.7390851332151607_f64;
        let root = newton(f_trig, 0.5, 100, 1e-9);
        assert!(
            approx_eq(root, expected, 1e-6),
            "Expected ≈ {expected:.9}, got {root:.9}"
        );
    }

    /// Starting exactly on the root should return it immediately (f(x) ≈ 0).
    #[test]
    fn test_start_at_root() {
        let root = newton(f_quadratic, 2f64.sqrt(), 100, 1e-9);
        assert!(
            approx_eq(root, 2f64.sqrt(), 1e-12),
            "Should return the root when started there"
        );
    }

    /// A function with no real root should exhaust iterations and return 0.0.
    #[test]
    fn test_no_convergence_returns_sentinel() {
        let result = newton(f_no_real_root, 1.0, 50, 1e-9);
        assert_eq!(result, 0.0, "Expected 0.0 sentinel on non-convergence");
    }

    /// Tighter tolerance should not degrade accuracy vs. loose tolerance.
    #[test]
    fn test_tight_tolerance_more_accurate() {
        let loose = newton(f_quadratic, 1.5, 100, 1e-4);
        let tight = newton(f_quadratic, 1.5, 100, 1e-12);
        let truth = 2f64.sqrt();
        assert!(
            (tight - truth).abs() <= (loose - truth).abs(),
            "Tighter tolerance should yield equal or better accuracy"
        );
    }
}
*/
