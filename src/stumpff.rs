pub fn stumpff_c(x: f64) -> f64 {
    match x {
        x if x == 0. => 0.5,
        x if x.abs() < 1e-6 => (1. / 2.) - (x / 24.) + (x * x / 720.),
        x if x > 0. => (1. - x.sqrt().cos()) / x, // small but nonzero z
        x if x < 0. => (1. - (-x).sqrt().cosh()) / x,
        _ => panic!("stumpff_c panicked"),
    }
}

pub fn stumpff_s(x: f64) -> f64 {
    match x {
        x if x == 0. => 1. / 6.,
        x if x.abs() < 1e-6 => (1. / 6.) - (x / 120.) + (x * x / 5040.), // small but nonzero z
        x if x > 0. => (x.sqrt() - x.sqrt().sin()) / x.sqrt().powi(3),
        x if x < 0. => ((-x).sqrt().sinh() - (-x).sqrt()) / (-x).sqrt().powi(3),
        _ => panic!("stumpff_s panicked"),
    }
}

// -------------------------------
//             TESTS
// -------------------------------

#[cfg(test)]
mod stumpff_tests {
    use super::*;
    use std::f64::consts::PI;

    const TOL: f64 = 1e-12;
    const TOL_SERIES: f64 = 1e-10;

    fn assert_close(actual: f64, expected: f64, tol: f64, label: &str) {
        let err = (actual - expected).abs();
        assert!(
            err < tol,
            "{label}: got {actual:.15e}, expected {expected:.15e}, diff {err:.3e}"
        );
    }

    // -----------------------------------------------------------------------
    // stumpff_c  —  C(ψ) = c₂(ψ)
    //
    //   ψ > 0 :  (1 − cos √ψ) / ψ
    //   ψ = 0 :  1/2
    //   ψ < 0 :  (1 − cosh √(−ψ)) / ψ
    // -----------------------------------------------------------------------

    #[test]
    fn test_c_elliptic() {
        // ψ = 1 : (1 − cos 1) / 1
        assert_close(stumpff_c(1.0), (1.0 - 1_f64.cos()) / 1.0, TOL, "C(1)");
        // ψ = 4 : (1 − cos 2) / 4
        assert_close(stumpff_c(4.0), (1.0 - 2_f64.cos()) / 4.0, TOL, "C(4)");
        // ψ = 9 : (1 − cos 3) / 9
        assert_close(stumpff_c(9.0), (1.0 - 3_f64.cos()) / 9.0, TOL, "C(9)");
        // ψ = π² : (1 − cos π) / π² = 2/π²
        assert_close(stumpff_c(PI * PI), 2.0 / (PI * PI), TOL, "C(π²)");
    }

    #[test]
    fn test_c_zero() {
        assert_close(stumpff_c(0.0), 0.5, TOL, "C(0)");
    }

    #[test]
    #[ignore = "WIP"]
    fn test_c_near_zero_continuity() {
        // Taylor:  C(ψ) ≈ 1/2 − ψ/24 + ψ²/720
        for &psi in &[1e-4_f64, 1e-6, 1e-8, -1e-4, -1e-6, -1e-8] {
            let expected = 0.5 - psi / 24.0 + psi.powi(2) / 720.0;
            assert_close(
                stumpff_c(psi),
                expected,
                TOL_SERIES,
                &format!("C({psi:.0e}) near-zero"),
            );
        }
    }

    #[test]
    fn test_c_hyperbolic() {
        // ψ = −1 : (1 − cosh 1) / (−1)
        assert_close(stumpff_c(-1.0), (1.0 - 1_f64.cosh()) / (-1.0), TOL, "C(-1)");
        // ψ = −4 : (1 − cosh 2) / (−4)
        assert_close(stumpff_c(-4.0), (1.0 - 2_f64.cosh()) / (-4.0), TOL, "C(-4)");
        // ψ = −9 : (1 − cosh 3) / (−9)
        assert_close(stumpff_c(-9.0), (1.0 - 3_f64.cosh()) / (-9.0), TOL, "C(-9)");
    }

    // -----------------------------------------------------------------------
    // stumpff_s  —  S(ψ) = c₃(ψ)
    //
    //   ψ > 0 :  (√ψ − sin √ψ) / ψ^(3/2)
    //   ψ = 0 :  1/6
    //   ψ < 0 :  (sinh √(−ψ) − √(−ψ)) / (−ψ)^(3/2)
    // -----------------------------------------------------------------------

    #[test]
    fn test_s_elliptic() {
        // ψ = 1 : (1 − sin 1) / 1
        assert_close(stumpff_s(1.0), (1.0 - 1_f64.sin()) / 1.0, TOL, "S(1)");
        // ψ = 4 : (2 − sin 2) / 8
        assert_close(stumpff_s(4.0), (2.0 - 2_f64.sin()) / 8.0, TOL, "S(4)");
        // ψ = 9 : (3 − sin 3) / 27
        assert_close(stumpff_s(9.0), (3.0 - 3_f64.sin()) / 27.0, TOL, "S(9)");
        // ψ = π² : (π − sin π) / π³ = 1/π²
        assert_close(stumpff_s(PI * PI), 1.0 / (PI * PI), TOL, "S(π²)");
    }

    #[test]
    fn test_s_zero() {
        assert_close(stumpff_s(0.0), 1.0 / 6.0, TOL, "S(0)");
    }

    #[test]
    #[ignore = "WIP"]
    fn test_s_near_zero_continuity() {
        // Taylor:  S(ψ) ≈ 1/6 − ψ/120 + ψ²/5040
        for &psi in &[1e-4_f64, 1e-6, 1e-8, -1e-4, -1e-6, -1e-8] {
            let expected = 1.0 / 6.0 - psi / 120.0 + psi.powi(2) / 5040.0;
            assert_close(
                stumpff_s(psi),
                expected,
                TOL_SERIES,
                &format!("S({psi:.0e}) near-zero"),
            );
        }
    }

    #[test]
    fn test_s_hyperbolic() {
        // ψ = −1 : (sinh 1 − 1) / 1
        assert_close(stumpff_s(-1.0), (1_f64.sinh() - 1.0) / 1.0, TOL, "S(-1)");
        // ψ = −4 : (sinh 2 − 2) / 8
        assert_close(stumpff_s(-4.0), (2_f64.sinh() - 2.0) / 8.0, TOL, "S(-4)");
        // ψ = −9 : (sinh 3 − 3) / 27
        assert_close(stumpff_s(-9.0), (3_f64.sinh() - 3.0) / 27.0, TOL, "S(-9)");
    }

    // -----------------------------------------------------------------------
    // Recurrence identities  —  derived directly from the definitions,
    // hold for all three branches, cross-check C and S independently.
    //
    //   ψ · C(ψ) = 1 − cos √ψ          (elliptic)
    //   ψ · C(ψ) = 1 − cosh √(−ψ)      (hyperbolic)
    //
    //   ψ · S(ψ) = 1 − sin(√ψ)/√ψ      (elliptic)
    //   ψ · S(ψ) = 1 − sinh(√(−ψ))/√(−ψ) (hyperbolic)
    // -----------------------------------------------------------------------

    #[test]
    fn test_recurrence_elliptic() {
        for &psi in &[0.5_f64, 1.0, 2.0, 4.0, 9.0, PI * PI] {
            let u = psi.sqrt();
            assert_close(
                psi * stumpff_c(psi),
                1.0 - u.cos(),
                TOL,
                &format!("ψ·C recurrence ψ={psi}"),
            );
            assert_close(
                psi * stumpff_s(psi),
                1.0 - u.sin() / u,
                TOL,
                &format!("ψ·S recurrence ψ={psi}"),
            );
        }
    }

    #[test]
    fn test_recurrence_hyperbolic() {
        for &psi in &[-0.5_f64, -1.0, -4.0, -9.0] {
            let u = (-psi).sqrt();
            assert_close(
                psi * stumpff_c(psi),
                1.0 - u.cosh(),
                TOL,
                &format!("ψ·C recurrence ψ={psi}"),
            );
            assert_close(
                psi * stumpff_s(psi),
                1.0 - u.sinh() / u,
                TOL,
                &format!("ψ·S recurrence ψ={psi}"),
            );
        }
    }
}
