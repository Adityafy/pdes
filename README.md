# PDEs & Dynamical Systems  
## PDEs like Generalized Swift–Hohenberg equation, Lyapunov Analysis, and CLVs

This repository contains numerical solvers and analysis tools for studying **spatiotemporal chaos**, **pattern formation**, and **tangent space dynamics** in nonlinear PDEs, with a focus on the Generalized Swift–Hohenberg (GSH) equation.

The framework supports:

- Fourier pseudospectral spatial discretization (periodic domains)
- Exponential time integration (ETD schemes)
- Tangent space evolution
- Lyapunov analysis
- Covariant Lyapunov Vector (CLV) computation

---

## Numerical Methods

### Spatial Discretization
- Fourier pseudospectral method  
- Periodic boundary conditions  
- FFT-based derivative evaluation  
- Optional dealiasing  

### Time Integration
- ETD1 and predictor–corrector exponential schemes  
- Unified integrator for flow and tangent dynamics  

### Tangent Space & Lyapunov Analysis
- Simultaneous evolution of multiple perturbations  
- QR orthogonalization and Gram-Schmidt vectors
- Lyapunov exponents computation  
- Dynamic algorithm for CLVs

---


## Dependencies

- MATLAB  
- FFT support (built-in)  
- No external toolboxes required  

---

## Goals of the Repository

This project is designed not just to simulate PDEs, but to study the **chaotic dynamics** in spatially extended systems.

It provides a complete computational pipeline from PDE integration to Lyapunov and CLV analysis, suitable for research in nonlinear dynamics and pattern formation.

---

## Author

Aditya Raj  
PhD Candidate
Virginia Tech
