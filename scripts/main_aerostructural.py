"""
### Info

Author: Jelle Poland \
Citing: https://doi.org/10.3390/en16145264 \
License: ... \
Github: ...
"""

### Initialisation

# Making things autoreload - needed for Jupyter Kernel/Interactive env.
%load_ext autoreload
%autoreload 2
%matplotlib widget

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
sys.path.append(project_root)  # Needed for running in terminal
os.chdir(project_root)  # Needed for running in interactive python environment

# Import modules
from src.initialisation import (
    load_surfplan,
    pulley_connectivity,
    input_classes,
    input_particleSystem,
    particles_with_rotational_resistance,
    actuation_relations,
    yaml_loader,
)
from src.initialisation.input_classes import (
    input_VSM,
    input_bridle_aero,
    input_structural_solver,
)
from src.initialisation.yaml_loader import config
from src.particleSystem import ParticleSystem
from src.coupling import coupling_aero2struc, coupling_struc2aero
from src.structural import structural_model, structural_mesher
from src.solver import solver_main, solver_utils
from src.post_processing import (
    functions_print,
    functions_plot,
    post_processing_utils,
    post_processing_main,
)
from src.aerodynamic import (
    VSM,
    breukels_2D,
    plate_aero,
    bridle_line_system_aero,
    tether_aero,
)

# Get mutable variables
from src.initialisation.mutable_variables import get_variables

points, vel_app, params_dict, initial_conditions, psystem = get_variables()

######################################
##### AeroStructural Simulation ######
######################################

## straight symmetric case
sim_name = "straight_symmetric"
points, print_data, plot_data, animation_data = solver_main.run_aerostructural_solver(
    points,
    vel_app,
    psystem,
    params_dict,
    config,
    input_VSM,
    input_bridle_aero,
    is_with_vk_optimization=False,
    is_circular_case=False,
    is_run_only_1_time_step=False,
    is_print_intermediate_results=True,
    is_with_gravity=True,
    sim_name="straight_symmetric",
    is_with_velocity_initialization=False,
)

# TODO:  Data-saving procedure for the torque paper, should be stream-lined
# data_names = [
#     "points",
#     "psystem",
#     "print_data",
#     "plot_data",
#     "animation_data",
#     "config",
#     "vel_app",
#     "sim_name",
#     "input_VSM",
# ]
# loaded_data = {}
# for name in data_names:
#     with open(
#         f"results/{config.kite_name}/torque_paper/2023_26_14_straight/{name}.pkl", "rb"
#     ) as f:
#         loaded_data[name] = dill.load(f)

# points = loaded_data["points"]
# psystem = loaded_data["psystem"]
# print_data = loaded_data["print_data"]
# plot_data = loaded_data["plot_data"]
# animation_data = loaded_data["animation_data"]
# config = loaded_data["config"]
# vel_app = loaded_data["vel_app"]
# sim_name = loaded_data["sim_name"]
# input_VSM = loaded_data["input_VSM"]
# folder_name = f"results/{config.kite_name}/torque_paper/2023_26_14_v1/"

# # circular case
# sim_name = "steady_circular"
# points, print_data, plot_data, animation_data = solver_main.run_aerostructural_solver(
#     points,
#     vel_app,
#     psystem,
#     params_dict,
#     config,
#     input_VSM,
#     input_bridle_aero,
#     is_with_vk_optimization=True,
#     is_circular_case=True,
#     is_run_only_1_time_step=False,
#     is_print_intermediate_results=True,
#     is_with_gravity=False,
#     sim_name="steady_circular",
#     is_with_velocity_initialization=False,
# )

# Output settings
is_with_printing = True
is_with_plotting = True
is_with_animation = True
is_with_save = False

if is_with_printing:
    post_processing_main.print_results(
        points,
        print_data,
        config,
    )
if is_with_plotting:
    post_processing_main.plot(
        plot_data,
        points,
        vel_app,
        config,
    )
if is_with_animation:
    post_processing_main.animate(
        animation_data,
        vel_app,
        sim_name,
        config,
        input_VSM,
    )
# if is_with_save:
#     to_be_saved_data = [
#         [points, "points"],
#         [psystem, "psystem"],
#         [print_data, "print_data"],
#         [plot_data, "plot_data"],
#         [animation_data, "animation_data"],
#         [config, "config"],
#         [vel_app, "vel_app"],
#         [sim_name, "sim_name"],
#         [input_VSM, "input_VSM"],
#     ]

#     # Ensure the folder exists
#     if not os.path.exists(folder_name):
#         os.makedirs(folder_name)

#     # Serialize and save all data with dill
#     for item in to_be_saved_data:
#         data, name = item
#         with open(f"{folder_name}{name}.pkl", "wb") as f:
#             dill.dump(data, f)
