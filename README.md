# N-Body Simulation

Simulation of planetary motion for the top 100 massive celestial bodies using Verlet's method.

## Project Overview
This project simulates the motion of celestial bodies, such as planets and moons, in a gravitational system. The goal is to validate and verify numerical methods for solving coupled second-order differential equations using Verlet integration. The project includes analysis of:
- Elliptical orbits
- Conservation of energy (kinetic, potential, and total)
- Conservation of angular momentum

### Features
- **Simulation Method**: Verlet integration for numerical integration.
- **Key Outputs**:
  - Position and velocity data for celestial bodies over time.
  - Kinetic energy, potential energy, and total energy evolution.
  - Angular momentum components in Cartesian coordinates.

---

## Compilation Instructions
To compile and run the program:

1. Compile the supporting module file (`model.f90`) first:
   ```bash
   gfortran -c model.f90 -fcheck=all -g
