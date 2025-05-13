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
from kitesim.post_processing import post_processing_main

# from initialisation.yaml_loader import config
# .input_classes #import input_VSM, input_bridle_aero
# from kitesim import initialisation.yaml_loader #import config
# from kitesim import initialisation.mutable_variables #import get_mutable_variables
# from kitesim import solver #import solver_main
# from kitesim import post_processing #import post_processing_main

# loading data

# variables to tune
date = "2024_03_25"
version = "v3"
data_names = [
    "points",
    # "psystem",
    # "print_data",
    # "plot_data",
    # "animation_data",
    # "config",
    # "vel_app",
    # "sim_name",
    # "input_VSM",
]
loaded_data = {}
for name in data_names:
    with open(
        f"data/output/TUDELFT_V3_KITE/{date}_{version}/{name}.pkl",
        "rb",
    ) as f:
        print(f"Loading {name} from file")
        loaded_data[name] = dill.load(f)

points = loaded_data["points"]
# psystem = loaded_data["psystem"]
# print_data = loaded_data["print_data"]
# plot_data = loaded_data["plot_data"]
# animation_data = loaded_data["animation_data"]
# config = loaded_data["config"]
# vel_app = loaded_data["vel_app"]
# sim_name = loaded_data["sim_name"]
# input_VSM = loaded_data["input_VSM"]

print(f"points: {points}")

## now also load points

file_path = f"data/output/TUDELFT_V3_KITE/points_up_0.npy"
points_up_0 = np.load(file_path)

# plot points in 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.scatter(
    points[:, 0],
    points[:, 1],
    points[:, 2],
    c="r",
    marker="o",
)
ax.scatter(
    points_up_0[:, 0],
    points_up_0[:, 1],
    points_up_0[:, 2],
    c="b",
    marker="o",
)
plt.title("3D Scatter plot of points")
plt.xlabel("X-axis")
plt.ylabel("Y-axis")
ax.set_zlabel("Z-axis")


# Function to set equal axis
def set_axes_equal(ax):
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    plot_radius = 0.5 * max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


set_axes_equal(ax)
plt.show()

# if config.is_with_plotting:
#     post_processing_main.plot(
#         plot_data,
#         points,
#         vel_app,
#         config,
#     )

# plt.show()
