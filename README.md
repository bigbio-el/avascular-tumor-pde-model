# Avascular Tumor PDE Model

MATLAB implementation of the mathematical models presented in:

**Lavrenteva, E., Theodoropoulos, C., & Binns, M. (2023).  
Analytical Models of Intra- and Extratumoral Cell Interactions at Avascular Stage of Growth in the Presence of Targeted Chemotherapy.  
Bioengineering, 10(3), 385.**  
https://doi.org/10.3390/bioengineering10030385

## Overview

This repository contains the MATLAB code associated with a system of
partial differential equation (PDE) models describing avascular tumor
growth and interactions among tumor cells, surrounding tissue, nutrient
availability, innate immune response, and chemotherapy.

Three progressively extended model variants are provided.

## Repository structure

### 1. Basic model

`models/01_basic_model/`

Describes interactions among:

- proliferating tumor cells (P)
- quiescent tumor cells (Q)
- necrotic cells (N)
- surrounding healthy tissue cells (S)
- nutrient supply (C)

Main MATLAB script:

`Dynamics_main_Final.m`

### 2. Model with innate immune response

`models/02_immune_response/`

Extends the basic model by including the innate immune cell population (L).

Main MATLAB script:

`Dynamics_main_Final_Immune.m`

### 3. Model with immune response and chemotherapy

`models/03_immune_response_chemotherapy/`

Extends the immune-response model by including chemotherapy drug dynamics (D)
and drug effects on tumor, surrounding tissue, and immune cells.

Main MATLAB script:

`ChemoTumorDynamicsChemoImmune.m`

## Numerical implementation

The PDE systems are spatially discretized using finite-difference
approximations and integrated in time using the explicit Euler method.

The helper files `G.m`, `H.m`, and `Z.m` define nutrient-dependent model
functions. `myfun.m` contains the discretized system equations for each
model variant. The chemotherapy model additionally uses `V.m`.

## Running the models

Open MATLAB and set the current folder to one of the individual model
directories.

For example, for the basic model:

```matlab
cd models/01_basic_model
Dynamics_main_Final

cd models/02_immune_response
Dynamics_main_Final_Immune

cd models/02_immune_response
Dynamics_main_Final_Immune

## Citation

If you use this code, please cite:

Lavrenteva, E.; Theodoropoulos, C.; Binns, M. Analytical Models of Intra-
and Extratumoral Cell Interactions at Avascular Stage of Growth in the
Presence of Targeted Chemotherapy. Bioengineering 2023, 10, 385.

https://doi.org/10.3390/bioengineering10030385

## Requirements

MATLAB

No toolbox-specific functions are evident in the supplied implementation.
