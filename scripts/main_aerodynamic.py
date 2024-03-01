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
sys.path.append(project_root)  # Needed for running in terminal
os.chdir(project_root)  # Needed for running in interactive python environment

# Import modules
from src.initialisation.input_classes import input_VSM, input_bridle_aero
from src.initialisation.yaml_loader import config
from src.coupling import coupling_aero2struc, coupling_struc2aero
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
from src.initialisation.mutable_variables import get_mutable_variables


def run_aerodynamic(points, vel_app, input_VSM, input_bridle_aero, config):
    """Run the aerodynamic calculations"""

    # Struc --> aero
    (
        points_wing_segment_corners_aero_orderded,
        index_transformation_struc_to_aero,
    ) = coupling_struc2aero.extract_wingpanel_corners_aero_orderded(
        points, config.kite.connectivity.plate_point_indices
    )
    # creating a dict with key value pairs, to transform from aero to struc
    index_transformation_aero_to_struc_dict = {}
    for i, value in enumerate(index_transformation_struc_to_aero):
        index_transformation_aero_to_struc_dict[value] = i

    # Struc --> aero
    points_wing_segment_corners_aero_orderded = points[
        index_transformation_struc_to_aero
    ]
    # Wing Aerodynamic
    (
        force_aero_wing_VSM,
        moment_aero_wing_VSM,
        F_rel,
        ringvec,
        controlpoints,
        wingpanels,
        rings,
        coord_L,
        coord_refined,
    ) = VSM.calculate_force_aero_wing_VSM(
        points_wing_segment_corners_aero_orderded, vel_app, input_VSM
    )
    # Aero --> struc
    if config.coupling_method == "NN":
        force_aero_wing = coupling_aero2struc.aero2struc_NN(
            config.aero.n_chordwise_aero_nodes,
            wingpanels,
            force_aero_wing_VSM,
            points_wing_segment_corners_aero_orderded,
            index_transformation_aero_to_struc_dict,
            points,
        )
    elif config.coupling_method == "MSc_Oriol":
        force_aero_wing = coupling_aero2struc.aero2struc(
            points,
            config.kite.connectivity.wing_ci,
            config.kite.connectivity.wing_cj,
            config.kite.connectivity.plate_point_indices,
            force_aero_wing_VSM,
            moment_aero_wing_VSM,
            ringvec,
            controlpoints,
        )
    else:
        raise ValueError("Coupling method not recognized; wrong name or typo")

    # Bridle Aerodynamics
    if config.is_with_aero_bridle:
        force_aero_bridle = (
            bridle_line_system_aero.calculate_force_aero_bridle_thedens2022(
                points, vel_app, input_bridle_aero
            )
        )
    else:
        force_aero_bridle = [0]
    force_aero = force_aero_wing + force_aero_bridle

    return (
        force_aero_wing,
        force_aero_bridle,
        force_aero,
        F_rel,
        wingpanels,
        controlpoints,
        rings,
        coord_L,
    )


def main():
    """main function"""

    # initialisation of mutable variables
    points, vel_app, params_dict, psystem = get_mutable_variables()
    # Running the aerodynamic calculations
    (
        force_aero_wing,
        force_aero_bridle,
        force_aero,
        F_rel,
        wingpanels,
        controlpoints,
        rings,
        coord_L,
    ) = run_aerodynamic(points, vel_app, input_VSM, input_bridle_aero, config)

    # Printing the results
    functions_print.print_settings(vel_app, config)
    functions_print.print_initial_kite_dimensions(config)
    functions_print.print_aero(
        points, vel_app, force_aero_wing, force_aero_bridle, config
    )

    # Plotting
    post_processing_main.plot(
        [
            [wingpanels, controlpoints, rings, coord_L, F_rel],
            config.kite.wing_rest_lengths_initial,
            config.kite.bridle_rest_lengths_initial,
        ],
        points,
        vel_app,
        config,
    )


if __name__ == "__main__":
    main()
