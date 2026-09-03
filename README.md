# Avascular Tumor PDE Model

MATLAB implementation of the mathematical models presented in:

**Lavrenteva, E., Theodoropoulos, C., & Binns, M. (2023).**  
**Analytical Models of Intra- and Extratumoral Cell Interactions at Avascular Stage of Growth in the Presence of Targeted Chemotherapy.**  
*Bioengineering*, **10**(3), 385.  
https://doi.org/10.3390/bioengineering10030385

## Overview

This repository contains the MATLAB code associated with a system of partial differential equation (PDE) models describing avascular tumor growth and interactions among tumor cells, surrounding healthy tissue, nutrient availability, innate immune response, and chemotherapy.

Three progressively extended model variants are provided:

1. Basic avascular tumor model
2. Model with innate immune response
3. Model with innate immune response and chemotherapy

## Repository Structure

```text
models/
├── 01_basic_model/
├── 02_immune_response/
└── 03_immune_response_chemotherapy/
```

### 1. Basic Model

`models/01_basic_model/`

The basic model describes interactions among:

- proliferating tumor cells (`P`)
- quiescent tumor cells (`Q`)
- necrotic cells (`N`)
- surrounding healthy tissue cells (`S`)
- nutrient supply (`C`)

Main MATLAB script:

`Dynamics_main_Final.m`

### 2. Model with Innate Immune Response

`models/02_immune_response/`

This model extends the basic formulation by including the innate immune cell population (`L`).

Main MATLAB script:

`Dynamics_main_Final_Immune.m`

### 3. Model with Immune Response and Chemotherapy

`models/03_immune_response_chemotherapy/`

This model extends the immune-response formulation by including chemotherapy drug dynamics (`D`) and the effects of chemotherapy on tumor cells, surrounding healthy tissue, and immune cells.

Main MATLAB script:

`ChemoTumorDynamicsChemoImmune.m`

## Numerical Implementation

The models are formulated as systems of partial differential equations describing the spatial and temporal dynamics of the interacting cell populations and other model components.

The PDE systems are spatially discretized using finite-difference approximations and integrated in time using an explicit Euler method.

The individual model directories also contain the MATLAB functions required by the corresponding simulation.

The helper files `G.m`, `H.m`, and `Z.m` define functions used in the model equations, while `myfun.m` contains the discretized equations for each model variant.

The chemotherapy model additionally contains `V.m`, which defines the chemotherapy input used in the simulation.

## Running the Models

Download or clone the repository and open MATLAB.

Run each model from its corresponding directory so that MATLAB can access the associated helper functions.

### Basic Model

```matlab
cd models/01_basic_model
Dynamics_main_Final
```

### Model with Innate Immune Response

```matlab
cd models/02_immune_response
Dynamics_main_Final_Immune
```

### Model with Immune Response and Chemotherapy

```matlab
cd models/03_immune_response_chemotherapy
ChemoTumorDynamicsChemoImmune
```

The scripts perform the numerical simulations and generate spatial profiles of the modeled cell populations and other system components at selected time points.

## Model Variables

The complete model includes the following variables:

- `P` — proliferating tumor cell population density
- `Q` — quiescent tumor cell population density
- `N` — necrotic cell population density
- `S` — surrounding healthy tissue cell population density
- `C` — nutrient supply
- `L` — innate immune cell population density
- `D` — chemotherapy drug concentration

The three model variants progressively introduce these components, allowing the effects of immune response and chemotherapy on avascular tumor dynamics to be investigated.

## Requirements

- MATLAB

## Publication

The mathematical formulation, model assumptions, parameter definitions, numerical methodology, validation, and interpretation of the simulations are described in the associated publication:

**Lavrenteva, E., Theodoropoulos, C., & Binns, M. (2023).**  
Analytical Models of Intra- and Extratumoral Cell Interactions at Avascular Stage of Growth in the Presence of Targeted Chemotherapy.  
*Bioengineering*, **10**(3), 385.

https://doi.org/10.3390/bioengineering10030385

## Citation

If you use this code or the associated mathematical model in your research, please cite:

> Lavrenteva, E.; Theodoropoulos, C.; Binns, M. Analytical Models of Intra- and Extratumoral Cell Interactions at Avascular Stage of Growth in the Presence of Targeted Chemotherapy. *Bioengineering* **2023**, *10*, 385. https://doi.org/10.3390/bioengineering10030385

## License

This repository is distributed under the MIT License. See the `LICENSE` file for details.
https://doi.org/10.3390/bioengineering10030385

## Requirements

MATLAB

No toolbox-specific functions are evident in the supplied implementation.
