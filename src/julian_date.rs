use nalgebra::Vector2;
// julian: converts date (MM/DD/YYYY HH:MM:SS UTC) to julian day
// where:
// 1901 ≤ y ≤ 2099
// 1 ≤ m ≤ 12
// 1 ≤ d ≤ 31
#[allow(unused)]
pub fn julian(yint: i32, mint: i32, dint: i32, hhint: i32, mmint: i32, ssint: i32) -> Vector2<f64> {
    let y = if yint >= 1901 && yint <= 2099 {
        yint as f64
    } else {
        panic!("YEAR value is out of bounds")
    };

    let m = if mint >= 1 && mint <= 12 {
        mint as f64
    } else {
        panic!("MONTH value is out of bounds")
    };

    let d = if dint >= 1 && dint <= 31 {
        dint as f64
    } else {
        panic!("DAY value is out of bounds")
    };

    let hh = if hhint >= 0 {
        hhint as f64
    } else {
        panic!("HOUR(S) value cannot be negative")
    };

    let mm = if mmint >= 0 && mmint <= 60 {
        mmint as f64
    } else {
        panic!("MINUTE(S) value is out of bounds")
    };

    let ss = if ssint >= 0 && ssint <= 60 {
        ssint as f64
    } else {
        panic!("SECOND(S) value is out of bounds")
    };

    let ut = hh + mm / 60. + ss / 3600.;

    println!("ut = {} hr(s)", ut);

    let j0 = (367. * y) - (((((m + 9.) / 12.) as i32 + yint) as f64 * 1.75) as i32) as f64
        + ((275. * m / 9.) as i32) as f64
        + d
        + 1721013.5;

    let jd = j0 + ut / 24.;

    println!("j0 = {} hr(s)", j0);

    Vector2::new(j0, jd)
}

//julian_centuries: the time in Julian centuries between and J2000 and the date in fn JULIAN
#[allow(unused)]
pub fn julian_centuries(j0: f64) -> f64 {
    // [Julian centuries / Cy]
    (j0 - 2451545.) / 36525.
}

// : the Greenwich sidereal time (θ_g0) at 0 hr UT
/*
pub fn local_sidereal_time(t0: f64, ut: f64, eastlong: f64) -> f64 {
    // the Greenwich sidereal time (θ_g0) at 0 hr UT [deg]
    let mut theta0 = 100.4606184 + (36000.77004 * t0) + (0.000387933 * t0.powi(2))
        - (2.583 * 10e-8 * t0.powi(3));
    while theta0 > 360. {
        theta0 = theta0 - 360.
    }

    // at any other universal time
    let thetag = theta0 + 360.98564724 * ut / 24.;

    // incorporate east longitude in deg to get local sidereal time
    let mut lst = thetag + eastlong;
    while lst > 360. {
        lst = lst - 360.
    }

    lst
}
*/
