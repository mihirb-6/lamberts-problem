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
use crate::elements::get_elements;
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
    //let opts = Opts::from_args();
    //println!("{:?}", opts);

    // read command line parse
    let args = Args::parse();

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
    let u = read_vectors_from_file(args.infile).unwrap();
    //println!("{:#?}", u);

    /* Input values using a JSON file */
    let body: String = args.body;
    let dt: f64 = u.time; // [s]
    let r1 = Vector3::new(u.x1, u.y1, u.z1);
    let r2 = Vector3::new(u.x2, u.y2, u.z2);
    let direction = if args.direction == "Retrograde" {
        Direction::Retrograde
    } else {
        Direction::Prograde
    };

    // Call lambert function run solver and obtain v1 and v2
    #[allow(unused)]
    let (v1, v2) = lambert(
        r1,
        r2,
        direction,
        dt,
        grav_param[&body],
        args.z_init,
        args.max_itrs,
    );
    // Using r1 & v1 OR r2 & v2, obtain a set of orbital elements
    let (a, elements, t_1, r_p, r_a) = get_elements(r1, v1, grav_param[&body]);

    // Not necessary but its cleaner
    let h = elements.x;
    let i = elements.y;
    let raan = elements.z;
    let e = elements.w;
    let w = elements.a;
    let theta = elements.b;

    // Orbital Elements Print Statement
    println!("--------------ORBITAL ELEMENTS-----------------");
    println!("ECCENTRICITY (e):                {:.4}", e);
    println!("SPECIFIC ANGULAR MOMENTUM (h):   {:.2} [km^2 s^-1]", h);
    println!("SEMIMAJOR AXIS (a):              {:.4} [km]", a);
    println!("PERIAPSIS (r_p):                 {:.4} [km]", r_p);
    println!("APOAPSIS (r_a):                  {:.4} [km]", r_a);
    println!("INCLINATION (i):                 {:.2} [°]", i.to_degrees());
    println!(
        "RA OF ASCENDING NODE (Ω):        {:.2} [°]",
        raan.to_degrees()
    );
    println!("ARGUMENT OF PERIAPSIS (ω)        {:.2} [°]", w.to_degrees());
    println!(
        "TRUE ANOMALY (θ)                 {:.2} [°]",
        theta.to_degrees()
    );
    println!("-----------------------------------------------");
    println!(
        "Perigee encounter in {:.1} s = {:.2} min = {:.2} hr(s)",
        t_1,
        t_1 / 60.,
        t_1 / 3600.
    );

    // JSON output handling
    if args.json == 'Y' {
        get_json(a, r_p, r_a, elements).expect("Unable to produce JSON files");
    }
}

// Function responsible for exporting output to a JSON file
fn get_json(a: f64, r_p: f64, r_a: f64, elements: Vector6<f64>) -> std::io::Result<()> {
    let values = OrbitalElements {
        semimajor_axis: a,
        periapsis: r_p,
        apoapsis: r_a,
        eccentricity: elements.w,
        angular_momentum: elements.x,
        inclination: elements.y.to_degrees(),
        raan: elements.z.to_degrees(),
        argument_of_periapsis: elements.a.to_degrees(),
        true_anomaly: elements.b.to_degrees(),
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
fn print_intro() {
    println!("*************************************************");
    println!("                  Lambert Solver                 ");
    println!("*************************************************");
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
