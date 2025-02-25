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
from pathlib import Path

from kitesim.initialisation import initialisation_main
from kitesim.solver import solver_main_pss_vsm
from kitesim.post_processing import post_processing_main
from kitesim.logging_config import *


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

    # Loading the input
    sim_input = initialisation_main.get_sim_input(
        path_config,
        path_processed_data_folder,
    )

    # Create results folder
    path_results_folder_run = post_processing_main.create_results_folder(
        sim_input, path_results_folder
    )

    # Save inputs
    post_processing_main.saving_all_dict_entries(
        sim_input, "input", path_results_folder_run
    )

    # AeroStructural Simulation
    sim_output = solver_main_pss_vsm.run_aerostructural_solver(
        sim_input,
    )

    # Save outputs
    post_processing_main.saving_all_dict_entries(
        sim_output, "output", path_results_folder_run
    )

    # Create interpretable results
    loaded_data_input = post_processing_main.processing_output(
        path_results_folder_run,
    )


if __name__ == "__main__":
    main()
