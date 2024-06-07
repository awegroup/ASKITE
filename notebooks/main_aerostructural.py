"""
### Info

Author: Jelle Poland \
Citing: https://doi.org/10.3390/en16145264 \
License: ... \
Github: ...
"""

### Initialisation

# Making things autoreload - needed for Jupyter Kernel/Interactive env.
# %load_ext autoreload
# %autoreload 2
# %matplotlib widget

import os
import sys

# from IPython.display import display, Latex

# TODO: can we remove this?!?
# Define the right path
# project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
# sys.path.append(f"{project_root}")  # Needed for running in terminal
# sys.path.insert(0, f"{project_root}")  # Needed for running in terminal
# os.chdir(f"{project_root}")  # Needed for running in interactive python environment

from pathlib import Path

from kitesim.initialisation import (
    InputVSM,
    InputBridleAero,
    setup_config,
    get_mutable_variables,
)
from kitesim.solver import solver_main
from kitesim.post_processing import main_post_processing


# Import modules
def main():
    """Main function"""
    # Find the root directory of the repository
    root_dir = os.path.abspath(os.path.dirname(__file__))
    while not os.path.isfile(os.path.join(root_dir, ".gitignore")):
        root_dir = os.path.abspath(os.path.join(root_dir, ".."))
        if root_dir == "/":
            raise FileNotFoundError(
                "Could not find the root directory of the repository."
            )
    # defining paths
    # path_config = "../data/config.yaml"
    path_config = Path(root_dir) / "data" / "config.yaml"
    # underlying mechanism assumes specific folder structure inside processed_data
    ## kite config files in folder: processed_data/kite_name
    ## kite data files in folder: processed_data/kite_name/processed_design_files
    path_processed_data_folder = Path(root_dir) / "processed_data"
    path_results_folder = Path(root_dir) / "results"

    config = setup_config(
        path_config,
        path_processed_data_folder,
    )
    input_VSM = InputVSM.create(config)
    input_bridle_aero = InputBridleAero.create(config)

    # Get mutable variables
    points_ini, vel_app, params_dict, psystem = get_mutable_variables(config)

    # Defining input_dict
    sim_input = {
        "points": points_ini,
        "vel_app": vel_app,
        "psystem": psystem,
        "params": params_dict,
        "config": config,
        "input_VSM": input_VSM,
        "input_bridle_aero": input_bridle_aero,
    }

    # Create results folder
    path_results_folder_run = main_post_processing.create_results_folder(
        config, path_results_folder
    )

    # Save inputs
    main_post_processing.saving_all_dict_entries(
        sim_input, "input", path_results_folder_run
    )

    # AeroStructural Simulation
    sim_output = solver_main.run_aerostructural_solver(
        sim_input,
    )

    # Save outputs
    main_post_processing.saving_all_dict_entries(
        sim_output, "output", path_results_folder_run
    )

    # Create interpretable results
    main_post_processing.processing_output(
        path_results_folder_run,
    )


if __name__ == "__main__":
    main()
