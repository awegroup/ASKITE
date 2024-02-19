# making things autorelad - needed for Jupyter Kernel
# %load_ext autoreload
# %autoreload 2
# %matplotlib widget


# Add the parent directory of `src` to sys.path
import os, sys

script_dir = os.path.dirname(__file__)  # Directory of the script
parent_dir = os.path.dirname(script_dir)  # Parent directory (project_root)
sys.path.append(parent_dir)

# Time to import the modules
from src.initialisation import load_surfplan, pulley_connectivity, actuation_relations

## All Immutatbles are stored in the dataclass config, i.e. simulation settings and configuration of the kite
from src.initialisation.yaml_loader import config
from src.initialisation.input_classes import (
    input_VSM,
    input_bridle_aero,
    input_structural_solver,
)
from src.initialisation import (
    input_particleSystem,
    particles_with_rotational_resistance,
)
from src.particleSystem.ParticleSystem import ParticleSystem
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

from test import test_main
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import time
import scipy.optimize
import yaml
import importlib
import pytest
import pandas as pd
import dill
from IPython.display import (
    display,
    Latex,
)  # {} gives units. {{}} is normal {} in Latexs


## Mutable variables
# Initializing Mutable Variables
points = config.kite.points_ini
# defining vel_app (vector of vel_app_norm)
vel_app = config.vel_wind - config.vel_kite

# Should be the same for each kite
(
    connectivity_matrix,
    wing_connectivity,
) = input_particleSystem.define_connectivity_matrix(config)
params_dict = input_particleSystem.define_params(
    config, wing_connectivity, connectivity_matrix
)
initial_conditions = input_particleSystem.define_initial_conditions_kite(config)
points_between_dict = particles_with_rotational_resistance.extract_points_between_dict(
    config
)
if config.is_with_initial_plot:
    functions_plot.plot_initial_geometry(config, points_between_dict)

is_with_rotational_resistance = False
if config.kite_name == "V9_60C":
    is_with_rotational_resistance = True

if is_with_rotational_resistance:
    (
        leadingedge_rotational_resistance_dict,
        strut_rotational_resistance_dict,
    ) = particles_with_rotational_resistance.extract_rotational_resistances_dicts(
        points_between_dict, config
    )
    # first do the struts
    k_bend_strut = 1e10
    params_dict = particles_with_rotational_resistance.initialize_bending_spring(
        k_bend_strut,
        initial_conditions,
        params_dict,
        connectivity_matrix,
        strut_rotational_resistance_dict,
    )
    # secondly do the leading-edge
    k_bend_leadingedge = 1e4
    params_dict = particles_with_rotational_resistance.initialize_bending_spring(
        k_bend_leadingedge,
        initial_conditions,
        params_dict,
        connectivity_matrix,
        leadingedge_rotational_resistance_dict,
    )
# Should be the same for each kite
psystem = ParticleSystem(connectivity_matrix, initial_conditions, params_dict)

# Printing initial dimensions
print(f"scaling-factor: {config.geometric_scaling_factor}")
print(f"ref_chord: {config.kite.ref_chord:.2f}m")
print(f"wing_span: {config.kite.span:.2f}m")
print(f"wing height: {config.kite.height:.2f}m")
print(f"wing area: {config.kite.area_surface:.2f}m2")
print(f"projected_area: {config.kite.area_projected:.2f}m")


def main():
    print("Starting the simulation")


if __name__ == "__main__":
    main()
