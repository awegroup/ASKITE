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
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import yaml
import importlib
import pytest
import pandas as pd
import dill
from IPython.display import display, Latex

# Define the right path
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(f"{project_root}")  # Needed for running in terminal
sys.path.insert(0, f"{project_root}")  # Needed for running in terminal
os.chdir(f"{project_root}")  # Needed for running in interactive python environment

from src.initialisation.input_classes import input_VSM, input_bridle_aero
from src.initialisation.yaml_loader import config
from src.initialisation.mutable_variables import get_mutable_variables
from src.solver import solver_main
from src.post_processing import post_processing_main


# Import modulesgit
def main():
    """Main function"""

    ##TODO: stream-line this better Data loading from another run
    date = "2024_04_25"
    folder_name = f"data/output/V3_25/torque_paper/{date}_v1/"
    data_names = [
        "points",
        "psystem",
        "print_data",
        "plot_data",
        "animation_data",
        "config",
        "vel_app",
        "sim_name",
        "input_VSM",
    ]
    loaded_data = {}
    for name in data_names:
        with open(f"{folder_name}/{name}.pkl", "rb") as f:
            loaded_data[name] = dill.load(f)

    points = loaded_data["points"]
    psystem = loaded_data["psystem"]
    print_data = loaded_data["print_data"]
    plot_data = loaded_data["plot_data"]
    animation_data = loaded_data["animation_data"]
    # config = loaded_data["config"]
    vel_app = loaded_data["vel_app"]
    sim_name = loaded_data["sim_name"]
    input_VSM = loaded_data["input_VSM"]

    # Get mutable variables
    points, vel_app, params_dict, psystem = get_mutable_variables()

    # AeroStructural Simulation
    points, print_data, plot_data, animation_data = (
        solver_main.run_aerostructural_solver(
            points,
            vel_app,
            psystem,
            params_dict,
            config,
            input_VSM,
            input_bridle_aero,
        )
    )

    ##TODO:  Data-saving procedure for the torque paper, should be stream-lined
    date = "2024_04_25"
    folder_name = f"data/output/{config.kite_name}/torque_paper/{date}_v2/"
    is_with_save = True
    if is_with_save:
        to_be_saved_data = [
            [points, "points"],
            [psystem, "psystem"],
            [print_data, "print_data"],
            [plot_data, "plot_data"],
            [animation_data, "animation_data"],
            [config, "config"],
            [vel_app, "vel_app"],
            [config.sim_name, "sim_name"],
            [input_VSM, "input_VSM"],
        ]

        # Ensure the folder exists
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)

        # Serialize and save all data with dill
        for item in to_be_saved_data:
            data, name = item
            with open(f"{folder_name}{name}.pkl", "wb") as f:
                dill.dump(data, f)

    # Post-processing output
    np.save(
        f"data/output/{config.kite_name}/points/{config.sim_name}_power_{1e3*config.depower_tape_final_extension:.0f}_steer_{1e3*np.abs(config.steering_tape_final_extension):.0f}.npy",
        points,
    )
    if config.is_with_printing:
        post_processing_main.print_results(
            points,
            print_data,
            config,
        )
    if config.is_with_plotting:
        post_processing_main.plot(
            plot_data,
            points,
            vel_app,
            config,
        )
    plt.show()
    if config.is_with_animation:
        print(f"")
        print("--> Generating ANIMATION \{*_*}/")
        post_processing_main.animate(
            animation_data,
            vel_app,
            config,
            input_VSM,
        )


if __name__ == "__main__":
    main()
