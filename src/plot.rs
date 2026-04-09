use nalgebra::{Matrix3, Vector3};
use plotters::prelude::*;

pub fn plot_orbit(e: f64, h: f64, i: f64, raan: f64, w: f64, mu: f64) {
    let theta_sweep: Vec<f64> = (0..500)
        .map(|k| k as f64 * 2.0 * std::f64::consts::PI / 500.0)
        .collect();

    // Compute radius
    let r_vals: Vec<f64> = theta_sweep
        .iter()
        .map(|&th| (h * h / mu) / (1.0 + e * th.cos()))
        .collect();

    // Perifocal coordinates
    let points_pf: Vec<Vector3<f64>> = theta_sweep
        .iter()
        .zip(r_vals.iter())
        .map(|(&th, &r)| Vector3::new(r * th.cos(), r * th.sin(), 0.0))
        .collect();

    // Rotation matrices
    let r3_raan = Matrix3::new(
        raan.cos(),
        -raan.sin(),
        0.0,
        raan.sin(),
        raan.cos(),
        0.0,
        0.0,
        0.0,
        1.0,
    );

    let r1_i = Matrix3::new(1.0, 0.0, 0.0, 0.0, i.cos(), -i.sin(), 0.0, i.sin(), i.cos());

    let r3_w = Matrix3::new(w.cos(), -w.sin(), 0.0, w.sin(), w.cos(), 0.0, 0.0, 0.0, 1.0);

    let transform = r3_raan * r1_i * r3_w;

    // Transform to ECI
    let points_eci: Vec<Vector3<f64>> = points_pf.iter().map(|r| transform * r).collect();

    // Set up plot
    let root = BitMapBackend::new("orbit.png", (1500, 1500)).into_drawing_area();
    root.fill(&BLACK).unwrap();

    let mut chart = ChartBuilder::on(&root)
        .margin(10)
        .caption("Orbit Plot", ("sans-serif", 20).into_font().color(&WHITE))
        .build_cartesian_2d(-1e6f64..1e6f64, -1e6f64..1e6f64)
        .unwrap();

    // Configure mesh with white grid lines and labels
    chart
        .configure_mesh()
        .x_labels(4)
        .y_labels(4)
        .axis_style(&WHITE) // axes lines
        .light_line_style(&RGBAColor(255, 255, 255, 0.1)) // grid lines
        .label_style(("sans-serif", 30).into_font().color(&WHITE)) // labels
        .draw()
        .unwrap();

    // Plot orbit (projected onto XY plane)
    chart
        .draw_series(LineSeries::new(
            points_eci.iter().map(|v| (v.x, v.y)),
            &GREEN,
        ))
        .unwrap();

    // main body
    chart
        .draw_series(std::iter::once(Circle::new((0.0, 0.0), 5, BLUE.filled())))
        .unwrap();

    root.present().unwrap();
}
