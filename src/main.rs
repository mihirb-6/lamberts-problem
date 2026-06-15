// nalgebra - helpful for vector math
use nalgebra::Vector3;
// for planetary ephemeris
#[allow(unused)]
use satkit::prelude::*;
// serde - helpful for JSON reading and writing
use serde::Deserialize;
use serde::Serialize;
// standard library - file management, helpful io functions, etc.
use std::collections::HashMap;
use std::error::Error;
use std::fs;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};
// clap - command-line parsing
use clap::Parser;

/* Mod [file] allows rust-analyzer to recognize them */
mod elements;
mod julian_date;
mod lagrange_coeffs;
mod lambert_eqns;
mod lambert_soln;
mod newton;
mod stumpff;
mod vectors;

// functions in other files that are used in main
use crate::elements::{Elements, get_elements};
use crate::lambert_soln::{Direction, lambert};

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    // prograde/retrograde option, default set to Prograde
    #[arg(short, long, default_value_t = String::from("Prograde"))]
    direction: String,

    // initial guess for z
    #[arg(short, long, default_value_t = 0.)]
    z_init: f64,

    // maximum iterations to run newton's method
    #[arg(short, long, default_value_t = 100)]
    max_itrs: u32,

    // input file containing vectors, time, and central body
    #[arg(short, long)]
    infile: PathBuf,

    // option to save orbital elements to JSON 'N'/'Y'
    #[arg(short, long, default_value_t = 'N')]
    json: char,

    // option to select central body
    // see grav_param for possible central bodies
    #[arg(short, long, default_value_t = String::from("Earth"))]
    body: String,
}

// the current implementation for an input data structure
// will likely change it once I figure out how to make it intuitive to use
#[derive(Deserialize, Debug)]
struct PositionVector {
    x1: f64,
    y1: f64,
    z1: f64,
    x2: f64,
    y2: f64,
    z2: f64,
    time: f64, // [s]
}

// a struct to neatly package the resulting data and export JSON later on
#[derive(Serialize, Debug)]
struct OrbitalElements {
    semimajor_axis: f64,
    periapsis: f64,
    apoapsis: f64,
    eccentricity: f64,
    angular_momentum: f64,
    inclination: f64,
    raan: f64,
    argument_of_periapsis: f64,
    true_anomaly: f64,
}

pub fn main() {
    // read command line parse
    let args = Args::parse();

    // title + banner

    // dictionary/hashmap storing key-value pairs of
    // central bodies and their gravitational parameter values
    let grav_param: HashMap<String, f64> = HashMap::from([
        (String::from("Sun"), 132712000000.),
        (String::from("Mercury"), 22030.),
        (String::from("Venus"), 324900.),
        (String::from("Earth"), 3.986004418e5),
        (String::from("Moon"), 4903.),
        (String::from("Mars"), 42828.),
        (String::from("Jupiter"), 126686000.),
        (String::from("Saturn"), 37931000.),
        (String::from("Uranus"), 5794000.),
        (String::from("Neptune"), 6835100.),
        (String::from("Pluto"), 830.),
    ]);

    // read the input JSON file
    let u = read_vectors_from_file(args.infile).unwrap();
    //println!("{:#?}", u);

    /* Input values using a JSON file */
    let central_body: String = args.body;
    let dt: f64 = u.time; // [s]
    let r1 = Vector3::new(u.x1, u.y1, u.z1);
    let r2 = Vector3::new(u.x2, u.y2, u.z2);

    let direction = if args.direction == "Retrograde" {
        Direction::Retrograde
    } else {
        Direction::Prograde
    };

    let mu = grav_param[&central_body];
    // Call lambert function run solver and obtain v1 and v2
    #[allow(unused)]
    let (v1, v2) = lambert(r1, r2, direction, dt, mu, args.z_init, args.max_itrs);
    // Using r1 & v1 OR r2 & v2, obtain a set of orbital elements
    let elements = get_elements(r1, v1, mu);

    // Orbital Elements Print Statement
    println!("--------------ORBITAL ELEMENTS-----------------");
    println!(
        "ECCENTRICITY (e):                {:.4}",
        elements.eccentricity
    );
    println!(
        "SPECIFIC ANGULAR MOMENTUM (h):   {:.2} [km^2 s^-1]",
        elements.angular_momentum
    );
    println!(
        "SEMIMAJOR AXIS (a):              {:.4} [km]",
        elements.seminajor_axis(mu)
    );
    println!(
        "PERIAPSIS (r_p):                 {:.4} [km]",
        elements.periapsis(mu)
    );
    println!(
        "APOAPSIS (r_a):                  {:.4} [km]",
        elements.apoapsis(mu)
    );
    println!(
        "INCLINATION (i):                 {:.2} [°]",
        elements.inclination.to_degrees()
    );
    println!(
        "RA OF ASCENDING NODE (Ω):        {:.2} [°]",
        elements.raan.to_degrees()
    );
    println!(
        "ARGUMENT OF PERIAPSIS (ω)        {:.2} [°]",
        elements.argument_of_periapsis.to_degrees()
    );
    println!(
        "TRUE ANOMALY (θ)                 {:.2} [°]",
        elements.true_anomaly.to_degrees()
    );
    println!("-----------------------------------------------");
    println!(
        "Perigee encounter in {:.1} s = {:.2} min = {:.2} hr(s)",
        elements.time_since_perapsis(mu),
        elements.time_since_perapsis(mu) / 60.,
        elements.time_since_perapsis(mu) / 3600.
    );

    // JSON output handling
    if args.json == 'Y' {
        get_json(
            elements.angular_momentum,
            elements.periapsis(mu),
            elements.apoapsis(mu),
            elements,
            mu,
        )
        .expect("Unable to produce JSON files");
    }
}

// Function responsible for exporting output to a JSON file
#[allow(unused)]
fn get_json(a: f64, r_p: f64, r_a: f64, elements: Elements, mu: f64) -> std::io::Result<()> {
    let values = OrbitalElements {
        semimajor_axis: elements.seminajor_axis(mu),
        periapsis: elements.periapsis(mu),
        apoapsis: elements.apoapsis(mu),
        eccentricity: elements.eccentricity,
        angular_momentum: elements.angular_momentum,
        inclination: elements.inclination.to_degrees(),
        raan: elements.raan.to_degrees(),
        argument_of_periapsis: elements.argument_of_periapsis.to_degrees(),
        true_anomaly: elements.true_anomaly.to_degrees(),
    };

    let j = match serde_json::to_string(&values) {
        Ok(v) => v,
        Err(_) => String::from("Unable to load data"),
    };

    println!("Successfully wrote to JSON.");

    fs::write("orbital_elements.json", j)?;
    Ok(())
}

// Title + banner
#[allow(unused)]
fn print_intro() {
    println!("*************************************************");
    println!("                  Lambert Solver                 ");
    println!("*************************************************");
}

// Function to read input JSON using serde and std lib
#[allow(unused)]
fn read_vectors_from_file<P: AsRef<Path>>(path: P) -> Result<PositionVector, Box<dyn Error>> {
    // Open the file in read-only mode with buffer.
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    // Read the JSON contents of the file as an instance of `User`.
    let u = serde_json::from_reader(reader)?;

    // Return the `User`.
    Ok(u)
}
