# Lorenz-System-
The Lorenz system is a system of ordinary differential equations first studied by mathematician and meteorologist Edward Lorenz. It is notable for having chaotic solutions for certain parameter values and initial conditions.
# Lorenz Toolbox

**Author:** Derek Martinez  
**Last updated:** 2025-06-17

---

## Overview

`LorenzAllInOne.m` is an **all-in-one MATLAB script** for exploring the Lorenz system. It combines an interactive baseline run with a full research suite:

1. **Baseline integration**  
   - Numeric solution of the Lorenz ODEs  
   - 3-D phase plot, animated trajectory, 2-D projections, and time-series  
   - Lyapunov exponent estimate (with two-trajectory fallback if needed)  
   - Optional correlation dimension (skipped if `corrDimension.m` is missing)  

2. **LHS parameter sweep**  
   - Latin-Hypercube sampling over σ, β, ρ  
   - Parallelized Lyapunov exponent map in a 3-D scatter  

3. **ρ-bifurcation scan**  
   - Sweep ρ over a user-defined range  
   - Plot of max/min z-values (bifurcation diagram)  

4. **Uncertainty quantification (UQ)**  
   - Ensemble of noisy initial conditions  
   - Bootstrap 95% CI on the mean Lyapunov exponent  

5. **Proper Orthogonal Decomposition (POD/SVD)**  
   - Compute dominant modes of the attractor  
   - Compare 2-mode reduced trajectory to the full system  

6. **Finite-difference sensitivities**  
   - Numerical gradient ∂λ/∂σ, ∂λ/∂β, ∂λ/∂ρ  

7. **Automatic output**  
   - Saves all figures as PNG in `Lorenz_Elite_Out/`  
   - Exports `workspace.mat` for post-processing  

---

## Requirements

- **MATLAB R2018b** or later  
- **Parallel Computing Toolbox** (for `parpool`)  
- **Statistics and Machine Learning Toolbox** (for `lhsdesign`)  
- (Optional) A `corrDimension.m` function on your path, if you want correlation-dimension estimates  

---

## Installation

1. Place `LorenzAllInOne.m` in a folder on your MATLAB path.  
2. (Optional) Add any helper functions (e.g. `corrDimension.m`) to the same folder.  
3. Start MATLAB and navigate to that folder.  

---

## Usage

1. **Run**  
   ```matlab
   >> LorenzAllInOne
