// nalgebra - helpful for vector math
use nalgebra::{Vector3, Vector6};
// serde - helpful for JSON reading and writing
use serde::Deserialize;
use serde::Serialize;
// standard library - file management, helpful io functions, etc.
use std::collections::HashMap;
use std::error::Error;
use std::fs;
use std::fs::File;
use std::io::{self, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::exit;
// structOpt - command-line parsing
use structopt::StructOpt;

/* Mod [file] allows rust-analyzer to recognize them */
mod elements;
mod julian_date;
mod lagrange_coeffs;
mod lambert_eqns;
mod lambert_soln;
mod newton;
mod plot;
mod stumpff;
mod vectors;

// functions in other files that are used in main
use crate::elements::get_elements;
use crate::lambert_soln::{Direction, lambert};

#[allow(unused)]
use crate::plot::plot_orbit; // for plotting, might use Python instead

// structOpt - commandline parsing struct
// infile - input file that will be read by serde-json
#[derive(Debug, StructOpt)]
struct Opts {
    #[structopt(parse(from_os_str))]
    infile: PathBuf,

    // DEAD CODE FOR NOW ==> TODO!
    #[structopt(short, long, parse(from_os_str))]
    outfile: Option<PathBuf>,
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
    hours: f64,
    planet: String,
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
    let opts = Opts::from_args();
    println!("{:?}", opts);

    // title + banner
    print_intro();

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
    let u = read_vectors_from_file(opts.infile).unwrap();
    println!("{:#?}", u);

    /* Input values using a JSON file */
    let planet: String = u.planet;
    let dt: f64 = u.hours * 3600.; // [s]
    let r1 = Vector3::new(u.x1, u.y1, u.z1);
    let r2 = Vector3::new(u.x2, u.y2, u.z2);

    /* Input values one by one using the terminal
    let (planet, x1, y1, z1, x2, y2, z2, hours) = query_values();
    let planet: String = planet;
    let dt: f64 = hours * 3600.; // [s]
    let r1 = Vector3::new(x1, y1, z1);
    let r2 = Vector3::new(x2, y2, z2);
    */

    /* Input values within the code/IDE
    let planet: String = String::from("Mars");
    let dt: f64 = 500. * 3600.; // [s]
    let r1 = Vector3::new(4343.0, 3653., -6344.0);
    let r2 = Vector3::new(-1340.0, 6325., -7333.);
    */

    // Call lambert function run solver and obtain v1 and v2
    #[allow(unused)]
    let (v1, v2) = lambert(r1, r2, Direction::Retrograde, dt, grav_param[&planet]);
    // Using r1 & v1 OR r2 & v2, obtain a set of orbital elements
    let (a, elements, t_1, r_p, r_a) = get_elements(r1, v1, grav_param[&planet]);

    // Not necessary but its cleaner
    let h = elements.x;
    let i = elements.y;
    let raan = elements.z;
    let e = elements.w;
    let w = elements.a;
    let theta = elements.b;

    //plot_orbit(e, h, i, raan, w, grav_param[planet], a);

    // Orbital Elements Print Statement
    println!("--------------ORBITAL ELEMENTS----------------");
    println!("SEMIMAJOR AXIS (a): = {:.4} [km]", a);
    println!("PERIAPSIS (r_p): = {:.4} [km]", r_p);
    println!("APOAPSIS (r_a): = {:.4} [km]", r_a);
    println!("ECCENTRICITY (e): = {:.4}", e);
    println!("SPECIFIC ANGULAR MOMENTUM (h): = {:.2} [km^2 s^-1]", h);
    println!("INCLINATION (i): = {:.2}°", i.to_degrees());
    println!("RA OF ASCENDING NODE (Ω): = {:.2}°", raan.to_degrees());
    println!("ARGUMENT OF PERIAPSIS (ω) = {:.2}°", w.to_degrees());
    println!("TRUE ANOMALY (θ) = {:.2}°", theta.to_degrees());
    println!("----------------------------------------------");
    println!(
        "Perigee encounter in {:.1} s = {:.2} min = {:.2} hr(s)",
        t_1,
        t_1 / 60.,
        t_1 / 3600.
    );

    /* Asking if use wants to save the elements/output to a JSON file */
    print!("EXTRACT TO ELEMENTS TO JSON? (Y/N): ");
    io::stdout().flush().expect("Failed to flush stdout.");
    let mut choice = String::new();
    io::stdin()
        .read_line(&mut choice)
        .expect("Read line error.");
    let choice: char = choice.trim().parse().expect("Unable to parse");

    if !matches!(choice, 'Y' | 'N') {
        println!("'Y' or 'N' only");
        exit(000001);
    }

    if choice == 'Y' {
        get_json(a, r_p, r_a, elements).expect("Unable to produce JSON files");
    } else {
    }
}

// Function responsible for exporting output to a JSON file
#[allow(unused)]
fn get_json(a: f64, r_p: f64, r_a: f64, elements: Vector6<f64>) -> std::io::Result<()> {
    let values = OrbitalElements {
        semimajor_axis: a,
        periapsis: r_p,
        apoapsis: r_a,
        eccentricity: elements.w,
        angular_momentum: elements.x,
        inclination: elements.y,
        raan: elements.z,
        argument_of_periapsis: elements.a,
        true_anomaly: elements.b,
    };

    let j = match serde_json::to_string(&values) {
        Ok(v) => v,
        Err(_) => String::from("Unable to load data"),
    };

    println!("Wrote to json: {}", j);

    fs::write("orbital_elements.json", j)?;
    Ok(())
}

// Title + banner
fn print_intro() {
    println!("*************************************************");
    println!("                  Lambert Solver                 ");
    println!("*************************************************");
}

// Function asking user to input values in the terminal
#[allow(unused)]
fn query_values() -> (String, f64, f64, f64, f64, f64, f64, f64) {
    println!(
        "Choose main body -> Sun/Mercury/Venus/Earth/Moon/Mars/Jupiter/Saturn/Uranus/Neptune/Pluto"
    );
    io::stdout().flush().expect("Failed to flush stdout.");
    let mut planet = String::new();
    io::stdin()
        .read_line(&mut planet)
        .expect("Read line error.");
    let planet: String = planet
        .trim()
        .parse()
        .expect("Enter one of the objects listed in the prompt.");

    println!("Position Vector R1 in km...");
    print!("Input x1: ");
    io::stdout().flush().expect("Failed to flush stdout.");
    let mut x1 = String::new();
    io::stdin().read_line(&mut x1).expect("Read line error.");
    let x1: f64 = x1.trim().parse().expect("Enter a floating point number");

    print!("Input y1: ");
    io::stdout().flush().expect("Failed to flush stdout.");
    let mut y1 = String::new();
    io::stdin().read_line(&mut y1).expect("Read line error.");
    let y1: f64 = y1.trim().parse().expect("Enter a floating point number");

    print!("Input z1: ");
    io::stdout().flush().expect("Failed to flush stdout.");
    let mut z1 = String::new();
    io::stdin().read_line(&mut z1).expect("Read line error.");
    let z1: f64 = z1.trim().parse().expect("Enter a floating point number");

    println!("Position Vector R2 in km...");
    print!("Input x2: ");
    io::stdout().flush().expect("Failed to flush stdout.");
    let mut x2 = String::new();
    io::stdin().read_line(&mut x2).expect("Read line error.");
    let x2: f64 = x2.trim().parse().expect("Enter a floating point number");

    print!("Input y2: ");
    io::stdout().flush().expect("Failed to flush stdout.");
    let mut y2 = String::new();
    io::stdin().read_line(&mut y2).expect("Read line error.");
    let y2: f64 = y2.trim().parse().expect("Enter a floating point number");

    print!("Input z2: ");
    io::stdout().flush().expect("Failed to flush stdout.");
    let mut z2 = String::new();
    io::stdin().read_line(&mut z2).expect("Read line error.");
    let z2: f64 = z2.trim().parse().expect("Enter a floating point number");

    println!("Time between R1 and R2 (dt) in hours...");

    print!("Input dt: ");
    io::stdout().flush().expect("Failed to flush stdout.");
    let mut hours = String::new();
    io::stdin().read_line(&mut hours).expect("Read line error.");
    let hours: f64 = hours.trim().parse().expect("Enter a floating point number");

    (planet, x1, y1, z1, x2, y2, z2, hours)
}

// Function to read input JSON using serde and std lib
fn read_vectors_from_file<P: AsRef<Path>>(path: P) -> Result<PositionVector, Box<dyn Error>> {
    // Open the file in read-only mode with buffer.
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    // Read the JSON contents of the file as an instance of `User`.
    let u = serde_json::from_reader(reader)?;

    // Return the `User`.
    Ok(u)
}
