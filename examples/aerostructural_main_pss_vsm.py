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

from kitesim.solver import main_pss_vsm
from kitesim.post_processing import post_processing_main
from kitesim.logging_config import *
from kitesim.utils import load_and_save_config_files, load_sim_output

# TODO: could add a function somewhere that always finds the right root dir
# # Find the root directory of the repository
# root_dir = os.path.abspath(os.path.dirname(__file__))
# while not os.path.isfile(os.path.join(root_dir, ".gitignore")):
#     root_dir = os.path.abspath(os.path.join(root_dir, ".."))
#     if root_dir == "/":
#         raise FileNotFoundError(
#             "Could not find the root directory of the repository."
#         )


# Import modules
def main():
    """Main function"""
    PROJECT_DIR = Path(__file__).resolve().parents[1]

    # load and save config files
    config, config_kite, results_dir = load_and_save_config_files(PROJECT_DIR)
    print(f"config files saved in {results_dir}\n")

    # run AeroStructural simulation
    sim_output = main_pss_vsm.run_aerostructural_solver(
        config, config_kite, PROJECT_DIR, results_dir
    )

    # Load results
    meta_data_dict, track = load_sim_output(Path(results_dir) / "sim_output.h5")

    print(f"n_iter: {meta_data_dict['n_iter']}")
    # print(f"track: {track}")
    f_internal = track["f_residual"] - track["f_external"]

    # # Save outputs
    # post_processing_main.saving_all_dict_entries(
    #     sim_output, "output", path_results_folder_run
    # )

    # # Create interpretable results
    # loaded_data_input = post_processing_main.processing_output(
    #     path_results_folder_run,
    # )


if __name__ == "__main__":
    main()
