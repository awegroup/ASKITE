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
from kitesim.post_processing import post_processing_main


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

    # TODO: Save inputs

    # AeroStructural Simulation
    points, df_position, post_processing_data = solver_main.run_aerostructural_solver(
        points_ini,
        vel_app,
        psystem,
        params_dict,
        config,
        input_VSM,
        input_bridle_aero,
    )

    ## Save outputs (same folder)
    # TODO: Should this be placed inside the solver_main loop?
    # Saving non-interpretable results
    path_run_results_folder = post_processing_main.save_non_interpretable_results(
        config,
        input_VSM,
        input_bridle_aero,
        points_ini,
        vel_app,
        params_dict,
        psystem,
        points,
        df_position,
        post_processing_data,
        path_results_folder,
    )
    # Saving interpretable results (generated from the saved non-interpretable results)
    post_processing_main.save_interpretable_results(path_run_results_folder)


if __name__ == "__main__":
    main()
