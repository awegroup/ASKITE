# ASKITE: Aero-structural coupled kite simulator

[![Python package testing](https://github.com/jellepoland/kitesim/actions/workflows/testing.yml/badge.svg?branch=main)](https://github.com/jellepoland/kitesim/actions/workflows/testing.yml)


## Purpose

This code is designed for aero-structural coupled simulations, which can be utilized in the kite design process to determine the performance characteristics of soft-wing kites without the need for actual flight testing.
It includes the [TU Delft V3 Kite](https://awegroup.github.io/TUDELFT_V3_KITE/) as an example, and relies primarly on:


# Project Structure

| File/Directory         | Content |
|-----------------------|---------|
| data/                 | Input data (YAML, geometry, etc.) and simulation results. |
| docs/                 | Sphinx documentation and images. |
| examples/             | Main entry points for running and analyzing simulations. See below. |
| requirements.txt      | Project dependencies. |
| settings/             | Configuration files (e.g., config.yaml). |
| src/                  | Main Python package source code (see details below). |
| tests/                | Test scripts and utilities. |
| venv/                 | Python virtual environment (not tracked in git). |
| README.md             | This file. |

## Main Simulation Workflow

The typical workflow consists of:

1. **Running a simulation**: Use `examples/run_simulation.py` to set up, run, and save a new simulation. This script:
	- Loads configuration and geometry from YAML files in `data/`.
	- Initializes aerodynamic and structural models (using modules in `src/kitesim/`).
	- Runs the coupled aero-structural solver.
	- Saves results to an HDF5 file in the results directory.

2. **Loading and analyzing results**: Use `examples/load_simulation.py` to load a previously saved simulation output and visualize or post-process the results. This script:
	- Loads the HDF5 results file.
	- Prints summary information about the run.
	- Provides interactive 3D visualization of the kite structure and its evolution.

## Details for `src/` Directory

The `src/kitesim/` directory contains the main simulation modules:

| Module                        | Purpose |
|-------------------------------|---------|
| aerodynamic.py                | Aerodynamic model setup and calculations |
| aerostructural_coupled_solver.py | Main coupled solver for aero-structural simulation |
| structural_pss.py             | Structural solver using the Particle System Simulator (PSS) |
| structural_pyfe3d.py          | Structural solver using pyfe3d (alternative) |
| aero2struc.py                 | Mapping between aerodynamic and structural meshes |
| read_struc_geometry_yaml.py   | Utilities for reading geometry and configuration from YAML |
| plotting.py                   | 3D and interactive plotting utilities |
| tracking.py                   | Tracking and post-processing of simulation data |
| utils.py                      | General utilities (file I/O, config, etc.) |
| logging_config.py             | Logging setup for the package |
| solver.py, struc2aero.py      | Additional solver and mapping utilities |
| __init__.py                   | Package initialization |

## Example Usage

To run a new simulation:

```bash
python examples/run_simulation.py
```

To load and interactively analyze a previous simulation:

```bash
python examples/load_simulation.py
```

See the docstrings in those scripts for more details on their workflow and configuration.

---

## More details per folder
[data](doc/data/data.md)

**Detailed List for `src` Directory**
The `src` directory contains various Python packages involved in the project:
<!-- - [coconut](doc/coconut.md) - Description or purpose of the `coconut` package. -->
- [aerodynamic](doc/src/aerodynamic/aerodynamic.md) - Description or purpose of the `aerodynamic` package.
- [coupling](doc/src/coupling/coupling.md) - Description or purpose of the `coupling` package.
- [initialisation](doc/src/initialisation/initialisation.md) - Description or purpose of the `initialisation` package.
- [particleSystem](doc/src/particleSystem/particleSystem.md) - Description or purpose of the `particleSystem` package.
- [post_processing](doc/src/post_processing/post_processing.md) - Description or purpose of the `post_processing` package.
- [solver](doc/src/solver/solver.md) - Description or purpose of the `solver` package.
- [structural](doc/src/structural/structural.md) - Description or purpose of the `structural` package.


## Installation Instructions
1. Clone the repository:
    ```bash
    git clone https://github.com/ocayon/Vortex-Step-Method
    ```

2. Navigate to the repository folder:
    ```bash
    cd Vortex-Step-Method
    ```
    
3. Create a virtual environment:
   
   Linux or Mac:
    ```bash
    python3 -m venv venv
    ```
    
    Windows:
    ```bash
    python -m venv venv
    ```
    
5. Activate the virtual environment:

   Linux or Mac:
    ```bash
    source venv/bin/activate
    ```

    Windows
    ```bash
    .\venv\Scripts\activate
    ```

6. Install the required dependencies:

   For users:
    ```bash
    pip install .
    ```
        
   For developers:
    ```bash
    pip install -e .[dev]
    ```
    
    For ubuntu add:
    ```
    pip install pyqt5
    sudo apt install cm-super
    sudo apt install dvipng
   ```

7. To deactivate the virtual environment:
    ```bash
    deactivate
    ```

## Dependencies
See [pyproject.toml](pyproject.toml)

**Machine Learning**

The code base is adapted to work with a machine learning model trained on more than a hundred thousands Reynolds-average Navier Stokes (RANS) Computational Fluid Dynamics (CFD) simulations made for leading-edge inflatable airfoils, documented in the MSc. thesis of [K.R.G. Masure](https://resolver.tudelft.nl/uuid:865d59fc-ccff-462e-9bac-e81725f1c0c9), the [code base is also open-source accessible](https://github.com/awegroup/Pointwise-Openfoam-toolchain).

As the three trained models, for Reynolds number = 1e6, 5e6 and 1e7 are too large (~2.3GB) for GitHub, they have to be downloaded separately, and added to the `data/ml_models` folder. They are accessible trough [Zenodo](https://doi.org/10.5281/zenodo.16925758), and so is the [CFD data](https://doi.org/10.5281/zenodo.16925833) on which the models are trained. More description on its usages is found in [Airfoil Aerodynamics](docs/AirfoilAerodynamics.md).

## Contributing Guide
Please report issues and create pull requests using the URL:
```
https://github.com/ocayon/Vortex-Step-Method/
```

This is required because you cannot/should not do it using the URL
```
https://github.com/awegroup/Vortex-Step-Method
```

We welcome contributions to this project! Whether you're reporting a bug, suggesting a feature, or writing code, here’s how you can contribute:

1. **Create an issue** on GitHub
2. **Create a branch** from this issue
   ```bash
   git checkout -b issue_number-new-feature
   ```
3. --- Implement your new feature---
4. Verify nothing broke using **pytest**
```
  pytest
```
5. **Commit your changes** with a descriptive message
```
  git commit -m "#<number> <message>"
```
6. **Push your changes** to the github repo:
   git push origin branch-name
   
7. **Create a pull-request**, with `base:develop`, to merge this feature branch
8. Once the pull request has been accepted, **close the issue**

## Citation
If you use this project in your research, please consider citing it. 
Citation details can be found in the [CITATION.cff](CITATION.cff) file included in this repository.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

### Copyright

Copyright (c) 2023 Jelle Poland, TU Delft
Copyright (c) 2023 Jelle Poland, Patrick Roeleveld, TU Delft
