# Inverse Aerodynamic Design of a Converging-Diverging Nozzle

**Author:** Danial Moghadam  
**Discipline:** Computational Fluid Dynamics (CFD) & Aerodynamic Shape Optimization

## Overview
This repository contains an automated inverse-design computational pipeline for optimizing the internal wall geometry of a converging-diverging nozzle. Utilizing a custom-built 1D, incompressible, inviscid flow solver, a Genetic Algorithm (GA) is deployed to manipulate a 7th-degree Bézier curve boundary. 

The optimization objective is to achieve a targeted pressure distribution designed to mimic boundary layer protection mechanisms via a Stratford curve approximation. The framework bypasses non-physical local minima to converge on a fully optimized aerodynamic shape subject to a 2.0-meter packaging constraint.

## Key Features & Methodology

### 1. Custom Flow Solver (Karimian & Schneider Method)
Rather than relying on commercial black-box software, a custom numerical solver was developed from the ground up.
* Solves the 1D Euler momentum and continuity equations using a fully implicit, $2N \times 2N$ block tri-diagonal matrix assembly.
* Implements the **Karimian & Schneider method** to distinguish between the convecting mass-conserving velocity and the convected momentum velocity, successfully suppressing the non-physical "checkerboard problem" typical of colocated grids.
* Validated against exact analytical formulations (Bernoulli's equation and mass conservation) to ensure absolute physical accuracy.

### 2. Physics Anchoring & Target Synthesis
To ensure the inviscid solver generates a geometry capable of surviving in a real, viscous fluid, the target pressure curve was heavily constrained:
* **Bernoulli Anchoring:** The inlet static pressure is analytically pre-calculated to respect global energy conservation with the atmospheric exhaust.
* **Stratford Approximation:** The diverging section targets an exponential decay, allowing the flow to ride the absolute limit of wall shear stress ($\tau_w = 0$) for maximum pressure recovery without boundary layer separation.

### 3. Stochastic Optimization & Geometric Regularization
A Genetic Algorithm operates as the global search optimizer, evaluating the Root Mean Square Error (RMSE) of guessed geometries.
* **Bézier Parameterization:** The physical wall is governed by a parameter-independent 7th-degree Bézier curve, allowing for maximum geometric flexibility.
* **Roughness Penalty:** To prevent the algorithm from exploiting the inviscid solver via artificial high-frequency wall oscillations, a geometric roughness penalty—based on the second derivative of the control points—is applied.

## Results
The algorithm actively explores the design space and successfully avoids localized minima. The final optimized geometry dynamically constricts the throat to meet localized pressure requirements before flawlessly hugging the Stratford pressure recovery target. 

By embedding viscous constraints directly into the objective function, the inviscid framework successfully generates a design that promotes boundary layer separation avoidance with extreme computational efficiency.

## Repository Structure
* **`/Report`**: Contains the final detailed project report (PDF).
* **`/1D_Inviscid_Nozzle_Toolkit`**: The isolated 1D Karimian & Schneider solver and the analytical piecewise-parabolic validation script.
* **`/MATLAB_Optimization`**: The full inverse-design suite, including the Genetic Algorithm wrapper, the Bézier geometry generator, and the parameter independence study.
* **`/Images`**: Convergence plots, Pareto fronts, and final geometric overlays.
