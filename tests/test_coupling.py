import numpy as np

from src.coupling import coupling_aero2struc, coupling_struc2aero
from src.aerodynamic import (
    VSM,
    breukels_2D,
    plate_aero,
    bridle_line_system_aero,
    tether_aero,
)
from src.initialisation.yaml_loader import config

# Struc --> aero
(
    points_wing_segment_corners_aero_orderded,
    index_transformation_struc_to_aero,
) = coupling_struc2aero.extract_wingpanel_corners_aero_orderded(
    config.kite.points_ini, config.kite.connectivity.plate_point_indices
)
points_old_method = points_wing_segment_corners_aero_orderded
points_new_method = config.kite.points_ini[index_transformation_struc_to_aero]
print(f"If true its great! {np.allclose(points_old_method, points_new_method)}")

# Aero --> struc Nearest Neighbour Method

## CAN ALSO TEST COMPLETE FUNCTION AGAINST URI HIS METHOD?

# 1. check if correct midpoints are calculated by computing them differently for a given panel
# def calculate_midpoints(wingpanels):

# 2. check if all points are between the midpoints and length matches the number of chordwise nodes
# def interpolate_chordwise_points(midpoints, n_chordwise_aero_nodes):

# 3. distribution ?
# def generate_distribution(n_chordwise_aero_nodes):


# 4. check if the mapping is correct, maybe defining a super simple case?
# def map_aerodynamic_to_structural(

# 5. Check complete function, how?
# def aero2struc_NN(
#     n_chordwise_aero_nodes,
#     wingpanels,
#     force_aero_wing_VSM,
#     points_wing_segment_corners_aero_orderded,
#     index_transformation_aero_to_struc_dict,
#     points_structural_nodes,
# ):
