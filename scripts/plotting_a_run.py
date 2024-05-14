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
from src.initialisation.input_classes import input_VSM, input_bridle_aero
from src.initialisation.yaml_loader import config
from src.initialisation.mutable_variables import get_mutable_variables
from src.solver import solver_main
from src.post_processing import post_processing_main

# loading data

# variables to tune
date = "2024_03_25"
version = "v3"
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
    with open(
        f"data/output/{config.kite_name}/torque_paper/{date}_{version}/{name}.pkl",
        "rb",
    ) as f:
        loaded_data[name] = dill.load(f)

points = loaded_data["points"]
psystem = loaded_data["psystem"]
print_data = loaded_data["print_data"]
plot_data = loaded_data["plot_data"]
animation_data = loaded_data["animation_data"]
config = loaded_data["config"]
vel_app = loaded_data["vel_app"]
sim_name = loaded_data["sim_name"]
input_VSM = loaded_data["input_VSM"]

if config.is_with_plotting:
    post_processing_main.plot(
        plot_data,
        points,
        vel_app,
        config,
    )

plt.show()
