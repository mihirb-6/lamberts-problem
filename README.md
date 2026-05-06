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

## Usage

```bash
./target/release/lamberts-problem -i [input_filename].json
```

---

## Technical Approach

The implementation follows a classical formulation of Lambert’s problem, solving for transfer geometry and iteratively converging on a valid trajectory.

Key components include:

* Transfer angle determination via cross and dot products
* Handling edge cases (e.g., collinear vectors, 0°/360° transfer angles)
* Root-finding for convergence on orbital parameters

The design prioritizes transparency over abstraction, making the algorithm accessible for further experimentation and research.

---

## Motivation

This project serves as both:

* A **technical exercise in orbital mechanics**, and
* A **foundation for more advanced astrodynamics simulations**, including N-body systems and trajectory optimization

It reflects a broader commitment to developing tools for scientific computing and spaceflight analysis.

---

## Future Work

* Multi-revolution Lambert solutions
* Improved numerical stability and convergence methods
* Integration with ephemeris data
* Visualization of transfer orbits

---

## Contributing

Contributions are welcome. Please open an issue or submit a pull request with clear documentation and test cases.

---

## License

MIT License

---

## Acknowledgments

Inspired by classical astrodynamics literature and modern implementations of Lambert solvers used in trajectory design and mission planning.

---

[1]: https://mcaneff.github.io/pdfs/Orbital%20Mechanics%20Lamberts%20Problem.pdf?utm_source=chatgpt.com "ORBITAL DYNAMICS"
