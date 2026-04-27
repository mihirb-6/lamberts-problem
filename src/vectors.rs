/*
 * Manual vector magntide and dot and cross products
 * but nalgebra has 'easy to work with' vectors, so this file is useless
 */
#[allow(unused)]
pub fn cross_product(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        (a[1] * b[2]) - (a[2] * b[1]),
        (a[2] * b[0]) - (a[0] * b[2]),
        (a[0] * b[1]) - (a[1] * b[0]),
    ]
}

#[allow(unused)]
pub fn dot_product(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    (a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2])
}

#[allow(unused)]
pub fn magnitude(a: &[f64; 3]) -> f64 {
    (a[0].powi(2) + a[1].powi(2) + a[2].powi(2)).sqrt()
}

// -------------------------------
//             TESTS by claude so double check
// -------------------------------

#[cfg(test)]
mod vector_tests {
    use super::*;

    const EPSILON: f64 = 1e-10;

    fn assert_vec_eq(result: [f64; 3], expected: [f64; 3]) {
        for i in 0..3 {
            assert!(
                (result[i] - expected[i]).abs() < EPSILON,
                "Component {i}: expected {}, got {}",
                expected[i],
                result[i]
            );
        }
    }

    // --- cross_product ---

    #[test]
    fn test_cross_product_basis_vectors_xy() {
        // x̂ × ŷ = ẑ
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        assert_vec_eq(cross_product(&a, &b), [0.0, 0.0, 1.0]);
    }

    #[test]
    fn test_cross_product_basis_vectors_yz() {
        // ŷ × ẑ = x̂
        let a = [0.0, 1.0, 0.0];
        let b = [0.0, 0.0, 1.0];
        assert_vec_eq(cross_product(&a, &b), [1.0, 0.0, 0.0]);
    }

    #[test]
    fn test_cross_product_basis_vectors_zx() {
        // ẑ × x̂ = ŷ
        let a = [0.0, 0.0, 1.0];
        let b = [1.0, 0.0, 0.0];
        assert_vec_eq(cross_product(&a, &b), [0.0, 1.0, 0.0]);
    }

    #[test]
    fn test_cross_product_anticommutative() {
        // a × b = -(b × a)
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        let ab = cross_product(&a, &b);
        let ba = cross_product(&b, &a);
        assert_vec_eq(ab, [-ba[0], -ba[1], -ba[2]]);
    }

    #[test]
    fn test_cross_product_parallel_vectors_zero() {
        // Parallel vectors → zero vector
        let a = [1.0, 2.0, 3.0];
        let b = [2.0, 4.0, 6.0];
        assert_vec_eq(cross_product(&a, &b), [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_cross_product_self_is_zero() {
        let a = [3.0, -1.0, 4.0];
        assert_vec_eq(cross_product(&a, &a), [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_cross_product_general() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        // (2*6 - 3*5, 3*4 - 1*6, 1*5 - 2*4) = (-3, 6, -3)
        assert_vec_eq(cross_product(&a, &b), [-3.0, 6.0, -3.0]);
    }

    #[test]
    fn test_cross_product_result_orthogonal_to_inputs() {
        // a × b must be orthogonal to both a and b
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        let c = cross_product(&a, &b);
        assert!(dot_product(&a, &c).abs() < EPSILON, "c not orthogonal to a");
        assert!(dot_product(&b, &c).abs() < EPSILON, "c not orthogonal to b");
    }

    #[test]
    fn test_cross_product_with_zero_vector() {
        let a = [1.0, 2.0, 3.0];
        let zero = [0.0, 0.0, 0.0];
        assert_vec_eq(cross_product(&a, &zero), [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_cross_product_negative_components() {
        let a = [-1.0, 2.0, -3.0];
        let b = [4.0, -5.0, 6.0];
        // (2*6 - (-3)*(-5), (-3)*4 - (-1)*6, (-1)*(-5) - 2*4)
        // = (12 - 15, -12 + 6, 5 - 8) = (-3, -6, -3)
        assert_vec_eq(cross_product(&a, &b), [-3.0, -6.0, -3.0]);
    }

    // --- dot_product ---

    #[test]
    fn test_dot_product_orthogonal_vectors() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        assert!((dot_product(&a, &b) - 0.0).abs() < EPSILON);
    }

    #[test]
    fn test_dot_product_parallel_unit_vectors() {
        let a = [1.0, 0.0, 0.0];
        assert!((dot_product(&a, &a) - 1.0).abs() < EPSILON);
    }

    #[test]
    fn test_dot_product_general() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        // 1*4 + 2*5 + 3*6 = 4 + 10 + 18 = 32
        assert!((dot_product(&a, &b) - 32.0).abs() < EPSILON);
    }

    #[test]
    fn test_dot_product_commutative() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        assert!((dot_product(&a, &b) - dot_product(&b, &a)).abs() < EPSILON);
    }

    #[test]
    fn test_dot_product_with_zero_vector() {
        let a = [1.0, 2.0, 3.0];
        let zero = [0.0, 0.0, 0.0];
        assert!((dot_product(&a, &zero) - 0.0).abs() < EPSILON);
    }

    #[test]
    fn test_dot_product_negative_result() {
        let a = [1.0, 0.0, 0.0];
        let b = [-1.0, 0.0, 0.0];
        assert!((dot_product(&a, &b) - (-1.0)).abs() < EPSILON);
    }

    #[test]
    fn test_dot_product_equals_magnitude_squared() {
        // a · a = |a|²
        let a = [3.0, 4.0, 0.0];
        let mag = magnitude(&a);
        assert!((dot_product(&a, &a) - mag * mag).abs() < EPSILON);
    }

    // --- magnitude ---

    #[test]
    fn test_magnitude_unit_x() {
        let a = [1.0, 0.0, 0.0];
        assert!((magnitude(&a) - 1.0).abs() < EPSILON);
    }

    #[test]
    fn test_magnitude_3_4_0() {
        // Classic 3-4-5 right triangle in xy-plane
        let a = [3.0, 4.0, 0.0];
        assert!((magnitude(&a) - 5.0).abs() < EPSILON);
    }

    #[test]
    fn test_magnitude_3d_pythagorean() {
        // sqrt(1² + 2² + 2²) = sqrt(9) = 3
        let a = [1.0, 2.0, 2.0];
        assert!((magnitude(&a) - 3.0).abs() < EPSILON);
    }

    #[test]
    fn test_magnitude_zero_vector() {
        let a = [0.0, 0.0, 0.0];
        assert!((magnitude(&a) - 0.0).abs() < EPSILON);
    }

    #[test]
    fn test_magnitude_non_negative() {
        let a = [-3.0, -4.0, 0.0];
        assert!(magnitude(&a) >= 0.0);
        assert!((magnitude(&a) - 5.0).abs() < EPSILON);
    }

    #[test]
    fn test_magnitude_consistent_with_dot_product() {
        // |a| = sqrt(a · a)
        let a = [1.0, 2.0, 3.0];
        let expected = dot_product(&a, &a).sqrt();
        assert!((magnitude(&a) - expected).abs() < EPSILON);
    }
}
