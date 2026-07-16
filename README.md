# Conjugate Heat Transfer Simulation Tool (C++)

A Finite Volume Method (FVM) solver written in C++ to model steady-state, axisymmetric conjugate heat transfer in a cylindrical pipe carrying a hot fluid that loses heat to the surrounding air through the pipe wall.

## Problem Overview

The tool solves a fully coupled, non-linear heat transfer problem involving three simultaneous physical mechanisms:

- **Internal forced convection**: heat exchange between the fluid and the inner pipe wall, with a convection coefficient that depends on the local flow regime (laminar or turbulent), computed via Hausen, Dittus-Boelter, or Sieder-Tate correlations.
- **Axial conduction through the pipe wall**, governed by the wall's thermal conductivity.
- **External natural convection**: heat dissipation from the outer pipe surface to still ambient air, with a coefficient derived from the Rayleigh number for a horizontal cylinder.

Fluid and air properties (density, viscosity, specific heat, conductivity) are **not assumed constant** — they are computed as temperature-dependent polynomial/exponential correlations, which introduces a non-linear coupling between wall temperature, fluid temperature, and convection coefficients that must be solved iteratively.

## Methodology

- **Spatial discretization**: the pipe wall is divided into `n` control volumes along its length, using the Finite Volume Method.
- **Fluid march**: at each axial node, fluid temperature, pressure, and velocity are obtained by simultaneously integrating the energy, momentum, and continuity equations, solved locally via fixed-point iteration (since the convection coefficient depends on Reynolds, Prandtl, and Graetz numbers, which in turn depend on the not-yet-converged temperature).
- **Wall energy balance**: assembly of the `aW`, `aE`, `aP`, `bP` coefficients at each node, coupling both convection coefficients.
- **Linear system solution**: the resulting tridiagonal system is solved using the **TDMA algorithm** (Thomas algorithm).
- **Outer iterative loop**: since properties and convection coefficients depend on temperature, the fluid, external convection, and wall blocks are solved in a segregated loop until the maximum wall temperature change between iterations falls below a set tolerance.

## Verification

Two independent verification studies were carried out to validate the implementation:

**1. Global energy balance**

The heat lost by the fluid and the heat gained by the ambient air are computed independently and compared:

| Magnitude | Value (W) |
|---|---|
| Q̇_air | 975.180 |
| Q̇_water | -975.180 |

The exact match (in magnitude) confirms the FVM implementation is conservative by construction, with no spurious numerical generation or loss of energy.

**2. Comparison against an analytical reference solution**

For the base case (water, v_in = 1 m/s), compared against a reference solution assuming constant thermophysical properties and convection coefficients:

| Variable | Numerical | Analytical | Relative Error |
|---|---|---|---|
| T_out (°C) | 94.24 | 94.23 | 0.01% |
| P_out (Pa) | 190352.6 | 190425 | 0.04% |
| v_out (m/s) | 0.999 | 1.00 | 0.10% |
| Q̇_w (W) | 975.2 | 985 | 1.00% |

## Physical Studies

The tool was also used to explore the sensitivity of the wall temperature profile to key physical parameters:

- **Wall conductivity**: increasing it drastically has negligible effect on the profile, confirming axial conduction is not the limiting heat transfer mechanism in the base case.
- **Fluid inlet velocity**: reducing it 10x (to 0.1 m/s) causes the wall temperature to drop from ~94.3°C to ~87.3°C along the pipe, versus less than 1°C in the base case — driven by the drastically lower mass flow rate and longer residence time.
- **Ambient temperature**: lowering it to -40°C increases the driving force for external convection, though the effect is moderated by internal convection remaining the dominant thermal resistance.
- **Combined effect**: applying all three changes simultaneously produces a temperature drop far greater than the sum of the individual effects, revealing a synergistic interaction between the three mechanisms.

## Numerical verification: mesh independence

The solution was tested with 5, 10, and 50 control volumes. The resulting temperature profiles were nearly identical, confirming the solution is mesh-independent — consistent with the smooth, gradual nature of the temperature profile in this problem.

## Tech Stack

- **Language**: C++
- **Numerical methods**: Finite Volume Method, TDMA (tridiagonal matrix algorithm), fixed-point iteration
- **No external libraries** — implemented from scratch using the standard library only

## Context

Developed as part of the *Gas Dynamics and Heat & Mass Transfer* course (Numerical Methods module), BSc in Aerospace Vehicle Engineering, UPC-ESEIAAT.
