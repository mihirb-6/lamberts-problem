# Lambert’s Problem Solver (Rust)

A lightweight and focused implementation of **Lambert’s Problem** in Rust, designed to compute orbital transfer velocities between two position vectors given a time of flight.

---

## Overview

Lambert’s Problem is a cornerstone of astrodynamics: determining the orbit connecting two points in space within a specified time interval. This repository provides a foundational implementation, emphasizing ease of use via terminal.

This project reflects a foundational step toward trajectory design, orbit determination, and mission analysis.

---

## Features

* Deterministic solution for two-point boundary value problems
* Support for **prograde and retrograde trajectories**
* Vector-based implementation using modern Rust practices
* Clear separation of geometric and numerical components

---

## Installation

```bash
git clone https://github.com/mihirb-6/lamberts-problem.git
cd lamberts-problem # navigate to where you downloaded the code
cargo build # executable can be found in ~/target/debug/
# can also do 'cargo build --release' for a more optimized version (~/target/release/lamberts-problem)
```

---

## Example Command-line Interface Usage
(zsh since I'm working from a mac, for Windows it should be similar)
```zsh
# To access the full list of flags and inputs you can use:
./target/release/lamberts-problem --help

# The most basic execution with only an input JSON file with two position vectors and time of flight
./target/release/lamberts-problem -i [input_filename].json

# A more customized way to run the program (-j Y saves orbital elements to a JSON file)
./target/release/lamberts-problem -i [input_filename].json -b Mars -d Retrograde -z 1.5 -j Y
```

---

## Technical Approach

The implementation follows a classical formulation of Lambert’s problem, solving for transfer geometry and iteratively converging on a valid trajectory using Newton's method.

Key components include:

* Transfer angle determination via cross and dot products
* Handling edge cases (e.g., collinear vectors, 0°/360° transfer angles)
* Root-finding for convergence on orbital parameters


---

## Motivation

This project serves as both:

* A **technical exercise in orbital mechanics**, and
* A **foundation for more advanced astrodynamics simulations**, including N-body systems and trajectory optimization

It reflects a broader commitment to developing tools for scientific computing and spaceflight analysis.

---

## Future Work

* Improved error handling
* Improved numerical stability and convergence methods
* Integration with ephemeris data
* Visualization of standard + transfer orbits

---

## License

MIT License

---

## Acknowledgments

Inspired (and largely an implementation of) by Howard D. Curtis' Orbital Mechanics textbook and his chapter on Lambert's problem
---

[1]: Curtis, H.D. (2013) Orbital Mechanics for Engineering Students. Butterworth-Heinemann, Oxford.
https://doi.org/10.1016/B978-0-08-097747-8.00006-2
