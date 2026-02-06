---
title: Overview
nav_order: 2
---

# Overview

AFiD_ELIFPPS is a numerical solver for
Eulerian–Lagrangian simulations of incompressible flows
with dispersed point particles.

---

## Problem Scope

AFiD_ELIFPPS targets problems involving:
- Incompressible Newtonian fluids
- Dispersed, dilute particle phases with large density ratios
- Two-way coupled momentum exchange
- Periodic or wall-bounded box domains with open boundaries

Typical applications include:
- Particle-laden turbulent flows
- Sedimentation and settling problems

---

## Modeling Approach

### Fluid Phase
- Eulerian formulation
- Incompressible Navier–Stokes equations

### Particle Phase
- Lagrangian point-particle model
- Maxey-Riley equations (Only steady state drag currently implemented)

### Coupling
- Two-way momentum coupling between phases. Particle forces computed from slip, and modelled as explicit source terms in the flow momentum equations

Details of governing equations and force models
are provided in the [Physics](physics/index.md) section.

---

## Numerical Approach (High Level)

- Finite difference discretization using second order central difference scheme
- Structured Cartesian grids
- Semi-implicit time integration 
- MPI pencil decomposition based on 2DECOMP&FFT library

Details are documented in the
[Numerical Methods](numerics/index.md) section.

---

## Assumptions and Limitations

The solver assumes:
- Box (cuboidal) domains
- Cartesian grid with uniform grid in two directions
- Moderate Reynolds numbers
- Low particle volume fraction
- Point-particle approximation

Known limitations and planned extensions are documented in
[Limitations and Roadmap](roadmap.md).

---

## Relationship to Other Codes

The core of this solver borrows heavily from several influences ([AFiD_Ventilation](https://github.com/ygsvanshu/AFiD_Ventilation), [AFiD-MuRPhFi](https://github.com/chowland/AFiD-MuRPhFi), [AFiD-2.0](https://github.com/ygsvanshu/AFiD_2.0) and [AFiD-1.0](https://github.com/PhysicsofFluids/AFiD) to mention a few). This implementation is, in essence, a modified version of the original high-performance solver developed at Physics of Fluids group, University of Twente in collaboration with University of Rome (Tor Vergata) and SURF for simulating Rayleigh–Bénard convection. 

---

## How to Get Started

For build instructions and first runs, see
[Getting Started](getting-started.md).


 