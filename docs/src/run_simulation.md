# Running a KiteSim Simulation: Overview

This document explains the workflow and logic of the `run_simulation.py` script in the KiteSim project. It details how configuration files are loaded, how the simulation is executed, and how results are saved and can be post-processed.

---

## 1. Purpose of `run_simulation.py`

The script `run_simulation.py` is the main entry point for running a full aero-structural simulation of a kite using the KiteSim framework. It orchestrates the loading of configuration files, runs the coupled solver, and saves the results for further analysis or visualization.

---

## 2. Workflow Breakdown

### Step 1: Project Directory Setup
- The script determines the project root directory using the location of the script itself.

### Step 2: Load and Save Configuration Files
- Calls `load_and_save_config_files(PROJECT_DIR)` from `kitesim.utils`.
- Loads two YAML files:
  - `data/config.yaml` (main simulation settings)
  - `data/<kite_name>/config_kite.yaml` (kite geometry and physical properties)
- Copies these config files into a timestamped results directory for reproducibility.

### Step 3: Run the Aero-Structural Simulation
- Calls `aerostructural_coupled_solver.main(config, config_kite)`.
- This function:
  - Initializes the kite structure and bridle system (via `initialisation.py`)
  - Sets up the aerodynamic and structural solvers
  - Runs the coupled simulation loop, exchanging forces and displacements between the solvers
  - Tracks relevant data (positions, forces, etc.) and meta information (timing, convergence)

### Step 4: Save Results
- Results are saved as an HDF5 file (`sim_output.h5`) in the results directory using `save_results()`.
- The file contains both tracking arrays (time histories) and meta data (attributes).

### Step 5: Load Results for Post-Processing
- The script demonstrates how to reload the results using `load_sim_output()`.
- This enables further analysis, plotting, or animation of the simulation output.

---

## 3. File Redirection and Data Flow

- **Configuration**: All simulation parameters and geometry are defined in YAML files and loaded at runtime.
- **Results**: All outputs are written to a timestamped results folder, ensuring traceability and reproducibility.
- **Post-Processing**: The script is structured to allow easy extension for plotting, animation, or further analysis.

---

## 4. Extending the Script

- You can add custom post-processing steps after the simulation, such as:
  - Plotting force histories (`f_int`, `f_ext`, `f_residual`)
  - Animating the kite's motion over time
  - Exporting results to other formats

---

## 5. Summary of Main Functions Used

- `load_and_save_config_files(PROJECT_DIR)`: Loads and archives config files
- `aerostructural_coupled_solver.main(config, config_kite)`: Runs the main simulation
- `save_results(tracking_data, meta, h5_path)`: Saves results to disk
- `load_sim_output(h5_path)`: Loads results for analysis

---

## 6. Example Usage

To run a simulation, simply execute:

```bash
python examples/run_simulation.py
```

Results will be saved in a new folder under `results/<kite_name>/YYYY_MM_DD_HHMMh/`.

---

For more details on each module, see the documentation in the corresponding `docs/src/` subfolders.
