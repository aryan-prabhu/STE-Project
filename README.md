# Joint Source Localisation and Probabilistic Boundary Tracking in Atmospheric Dispersion

**M.Tech Research Project — Aerospace Dynamics and Control, IIT Bombay**  
**Author:** Aryan Prabhu | **Supervisor:** Dr. Rohit Nanavati

---

## Overview

This repository implements a Bayesian framework for Source Term Estimation (STE) in 2D atmospheric dispersion. The core scientific contribution is the **Dual Information Gain (DIG)** motion policy for UAV-based source localisation, which exploits the duality between the plume boundary (a level set of the concentration field) and the source parameters θ = (x_s, y_s, Q) as dual representations of the same advection-diffusion process.

## Research Questions

- **RQ1:** Is the coupled Bayesian formulation well-posed?
- **RQ2:** Does DIG outperform Entrotaxis under matched sensor budgets?
- **RQ3:** Does an optimal trade-off parameter λ* exist, and how does it scale with noise, wind variability, and prior width?

## Repository Structure

```
STE_project/
├── src/                  # Core classes (OOP, MATLAB R2025b)
│   ├── AdvectionDiffusion2D.m   # Eq. 1: 2D advection-diffusion PDE solver
│   ├── PointSensor.m            # Eq. 2: Bilinear interpolation + Gaussian noise
│   ├── LevelSetExtractor.m      # Eq. 3: Marching-squares boundary extraction
│   └── ParticleFilterSIR.m      # Eq. 4: Sequential Importance Resampling filter
├── lib/                  # Helper functions
│   ├── gaussianPlume.m          # Analytical Gaussian plume (validation reference)
│   └── plumeConc.m              # Per-particle concentration prediction
├── tests/                # Validation scripts
└── startup.m             # Path setup
```

## Current Status

Equations 1–4 implemented and validated at the single-step level:
- Forward model: first-order upwind advection, central-difference diffusion, CFL-guarded
- Sensor model: bilinear interpolation, Gaussian noise, validated statistically
- Level-set extraction: Marching Squares via `contourc`, physical-space conversion
- Particle filter: SIR with systematic resampling, Shannon entropy tracking

## Dependencies

- MATLAB R2025b
- No external toolboxes required

## How to Run

1. Clone the repository
2. Open MATLAB and run `startup.m` to configure the path
3. Run any script in `tests/` to reproduce validation results