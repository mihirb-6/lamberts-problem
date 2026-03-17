mod stumpff;
mod vectors;

fn main() {
    let _tof: f64 = 5.0 * 3600.; // [s]
    let _r1: [f64; 3] = [1., 2., 3.];
    let _r2: [f64; 3] = [4., 5., 6.];
}

#[allow(unused)]
enum Direction {
    Prograde,
    Retrograde,
}

#[allow(unused)]
fn lambert(tof: f64, r1_vector: [f64; 3], r2_vector: [f64; 3], direction: Direction) {
    const MU: f64 = 3.986004418e5; // [km^3 s^-2]

    let r1 = vectors::magnitude(&r1_vector);
    let r2 = vectors::magnitude(&r2_vector);

    // Compute dot product of r1 and r2
    let r1_dot_r2 = vectors::dot_product(&r1_vector, &r2_vector);

    // Compute the cross product of r1 and r2
    let r1_cross_r2 = vectors::cross_product(&r1_vector, &r2_vector);

    // Compute dot product of the cross product of r1 and r2
    let mag_of_r1_r2_cross = vectors::magnitude(&r1_cross_r2);

    // Compute tranfer angle using atan2: atan2( mag(cross(r1_v, r2_v)), dot(r1_v,r2_v) )
    // Raw output lies between 0 and pi: 0 <= dtheta <= pi
    let mut dtheta = mag_of_r1_r2_cross.atan2(r1_dot_r2);

    /*
     * If motion is PROGRADE:
     *                  - keep dtheta as it is
     * If motion is RETROGRADE or long-way:
     *                  - define dtheta = 2pi - dtheta
     */

    match direction {
        Direction::Prograde => {}
        Direction::Retrograde => {
            dtheta = 2. * std::f64::consts::PI - dtheta;
        }
    }
}
