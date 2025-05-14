import logging
import numpy as np
from PSS.particleSystem.SpringDamper import SpringDamperType
from PSS.particleSystem import ParticleSystem
from pathlib import Path
import yaml
import attr
import os
import numpy as np
from pathlib import Path

import importlib.util
from scipy.spatial import ConvexHull
import numpy as np
import sys


def string_to_springdampertype(link_type: str) -> SpringDamperType:
    """
    Convert a string representation of a link type to a SpringDamperType enum value.

    Args:
        link_type (str): String representation of the link type.

    Returns:
        SpringDamperType: Corresponding enum value.

    Raises:
        ValueError: If the input string doesn't match any SpringDamperType.
    """
    try:
        return SpringDamperType(link_type.lower())
    except ValueError:
        raise ValueError(f"Invalid link type: {link_type}")


def define_kite_connectivity(config):
    wing_connectivity = np.column_stack(
        (config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj)
    )

    bridle_connectivity = np.column_stack(
        (config.kite.connectivity.bridle_ci, config.kite.connectivity.bridle_cj)
    )
    print(f"wing_ci {config.kite.connectivity.wing_ci}")
    print(f"wing_cj {config.kite.connectivity.wing_cj}")
    print(f"bridle_ci {config.kite.connectivity.bridle_ci}")
    print(f"bridle_cj {config.kite.connectivity.bridle_cj}")
    breakpoint()
    kite_connectivity = np.vstack((wing_connectivity, bridle_connectivity))

    return kite_connectivity, wing_connectivity


def define_initial_conditions_kite(config):
    # defining parameters
    points_ini = np.array(config.kite.points_ini)
    if config.is_with_initial_point_velocity:
        print("Error: initial point velocity has never been defined")
    else:
        vel_ini = np.zeros(points_ini.shape)
    m_array = config.kite.mass_points
    fixed_nodes = np.array(config.kite.bridle.bridle_point_index)
    # fill with: position, initial velocity?, mass, fixed boolean
    conditions = []
    n = config.kite.n_points
    for i in range(n):
        if i in fixed_nodes:
            conditions.append([points_ini[i], vel_ini[i], m_array[i], True])
        else:
            conditions.append([points_ini[i], vel_ini[i], m_array[i], False])

    return conditions


def define_param_dict_input_to_pss(config, wing_connectivity, kite_connectivity):
    bridle_rest_lengths = config.kite.bridle_rest_lengths_initial
    wing_rest_lengths = config.kite.wing_rest_lengths_initial
    rest_lengths = np.concatenate((wing_rest_lengths, bridle_rest_lengths))

    n_wing_elements = len(wing_connectivity)

    # Transform pulley.other_line_pair to a dict
    # data_struc is [["3",value],["5", value], ...]
    other_line_pair_dict = {}
    for entry in config.kite.pulley.other_line_pair:
        other_line_pair_dict[entry[0]] = entry[1]

    # Correct the key
    other_line_pair_corrected_dict = {}
    for key_i in other_line_pair_dict.keys():
        corrected_key = str(int(key_i) + n_wing_elements)
        other_line_pair_corrected_dict[corrected_key] = other_line_pair_dict[key_i]

    # initializing connectivity lists
    canopy_indices = []
    tube_indices = [i for i, conn in enumerate(wing_connectivity)]

    # initializing empty lists
    is_compression_list = []
    is_tension_list = []
    is_rotational_list = []
    stiffness_list = []
    is_pulley_list = []

    for i, conn in enumerate(kite_connectivity):
        # if wing elements
        if float(i) < len(wing_connectivity):
            is_tension_list.append(True)
            is_pulley_list.append(False)

            if i in config.kite.connectivity.te_line_indices:
                stiffness_list.append(config.kite.stiffness.trailing_edge)
                is_compression_list.append(False)
                is_rotational_list.append(False)
            elif i in config.kite.connectivity.tube_line_indices:
                stiffness_list.append(config.kite.stiffness.tube)
                is_compression_list.append(True)
                is_rotational_list.append(True)
            else:  # must be a canopy element
                stiffness_list.append(config.kite.stiffness.canopy)
                is_compression_list.append(False)
                is_rotational_list.append(False)

        # if bridle-lines
        else:
            is_compression_list.append(False)
            is_tension_list.append(True)
            is_rotational_list.append(True)
            stiffness_list.append(config.kite.stiffness.bridle)

            # TODO: might be better to use one index style/structure?
            # the (i - n_wing_elements) is to correct the indices
            if (i - n_wing_elements) in config.kite.pulley.line_indices:
                # print(f'pulley, i: {i}')
                # print(f'connection[i]: {kite_connectivity[i]}')
                is_pulley_list.append(True)
            else:
                is_pulley_list.append(False)

    params = {
        "c": config.solver.damping_constant,
        "dt": config.solver.dt,
        "t_steps": config.solver.n_time_steps,
        "abs_tol": config.solver.abs_tol,
        "rel_tol": config.solver.rel_tol,
        "max_iter": config.solver.max_iter,
        "pulley_other_line_pair": other_line_pair_corrected_dict,
        "k": np.array(stiffness_list),
        "is_compression": np.array(is_compression_list),
        "is_tension": np.array(is_tension_list),
        "is_pulley": np.array(is_pulley_list),
        "is_rotational": np.array(is_rotational_list),
        "n": int(len(config.kite.points_ini)),
        # "m_segment": 0.1,
        "aerostructural_tol": config.aero_structural.tol,  # N, < f_residual
        "l0": rest_lengths,
        "is_with_visc_damping": config.solver.is_with_visc_damping,
    }

    return params


## defining a distance to line function
def distance_point_to_line(point, line):
    """calculates distance between point and line
    return distance, vector and angle"""
    # defining p1, p2, p3
    p1, p2, p3 = line[0], line[1], point
    # defining the vector from p1 to p3
    p1_p3 = p3 - p1
    p1_p3_norm = np.linalg.norm(p1_p3)
    # defining the line vector from p1 to p2
    p1_p2 = p2 - p1
    p1_p2_norm = np.linalg.norm(p1_p2)
    p1_p2_unit = p1_p2 / p1_p2_norm

    # calculating the dot product and realize its equal to p1 to intersection point ps
    # dot(A,B) = |A|*|B|*cos(theta)
    # |B|*cos(theta) = dot(A,B) / |A|
    # take A as the line vector (p1_p2) and B as the diagonal (p1_p3)
    # |p1_p3|*cos(theta) = dot(p1_p2,p1_p3) / |p1_p2|

    # if we think trigonometry
    # cos(theta) = adjacent / diagonal
    # diagonal*cos(theta) = adjacent
    # |p1_p3| *cos(theta) = p1_ps_norm

    # ---> thus p1_ps_norm = dot(p1_p2,p1_p3) / |p1_p2|

    p1_ps_norm = np.dot(p1_p2, p1_p3) / p1_p2_norm
    # find the intersection point, by knowing it lies on the line
    ps = p1 + p1_p2_unit * p1_ps_norm
    # calculate the distance vector
    p3_ps = ps - p3
    # calculate the distance
    p3_ps_norm = np.linalg.norm(p3_ps)
    # calculate the angle between p1_p3 and p1_ps
    cos_angle_p1_p3_and_p1_ps = np.dot(p1_p3, p1_p2) / (p1_p2_norm * p1_p3_norm)
    # avoiding arccos input errors, it only takes -1< x <1
    if cos_angle_p1_p3_and_p1_ps > (1 - 1e-10):
        cos_angle_p1_p3_and_p1_ps = 1 - 1e-10
    elif cos_angle_p1_p3_and_p1_ps < (-1 + 1e-10):
        cos_angle_p1_p3_and_p1_ps = -1 + 1e-10
    angle_p1_p3_and_p1_ps = np.arccos(cos_angle_p1_p3_and_p1_ps)

    distance, vector, angle = p3_ps_norm, p3_ps, angle_p1_p3_and_p1_ps
    return distance, vector, angle


def extract_points_between_dict(config):
    # acquiring the attachment points that lie in between
    points_attachment_inbetween = []
    points_index_attachment_inbetween = []
    for idx in config.kite.connectivity.wing_ci[
        config.kite.connectivity.tube_line_indices
    ]:
        # if it is not a TE point
        if (
            idx
            not in config.kite.connectivity.wing_ci[
                config.kite.connectivity.te_line_indices
            ]
            and idx
            not in config.kite.connectivity.wing_cj[
                config.kite.connectivity.te_line_indices
            ]
            and idx not in config.kite.connectivity.plate_point_indices
        ):  # and not a plate point
            points_attachment_inbetween.append(config.kite.points_ini[idx])
            points_index_attachment_inbetween.append(idx)

    # need to sort these points per strut
    points_between_dict = {}
    tol_linepoint = 0.2
    for plate_index in config.kite.connectivity.plate_point_indices:
        left_le = config.kite.points_ini[plate_index[0]]
        left_te = config.kite.points_ini[plate_index[3]]
        for idx, point in enumerate(points_attachment_inbetween):
            if (
                abs(distance_point_to_line(point, [left_le, left_te])[0])
                < tol_linepoint
            ):
                points_between_dict[str(points_index_attachment_inbetween[idx])] = [
                    plate_index[0],
                    plate_index[3],
                ]

    # treating last strut separately
    right_le = config.kite.points_ini[
        config.kite.connectivity.plate_point_indices[-1][1]
    ]
    right_te = config.kite.points_ini[
        config.kite.connectivity.plate_point_indices[-1][2]
    ]
    for idx, point in enumerate(points_attachment_inbetween):
        if abs(distance_point_to_line(point, [right_le, right_te])[0]) < tol_linepoint:
            points_between_dict[str(points_index_attachment_inbetween[idx])] = [
                config.kite.connectivity.plate_point_indices[-1][1],
                config.kite.connectivity.plate_point_indices[-1][2],
            ]

    # the indices are sorted from left-to-right
    # the indices_mid are even sorted from TE to LE
    indices_le = [int(value[0]) for value in points_between_dict.values()]
    indices_mid = [int(key) for key in points_between_dict.keys()]
    indices_te = [int(value[1]) for value in points_between_dict.values()]

    return points_between_dict


def calculate_edge_lengths(ci, cj, pos):
    """returns the edge lengths between the nodes with index ci and cj
    for the given positions pos
    input : ci,cj,pos
    output: springL"""
    springL = np.zeros(ci.shape)
    for idx, (ci, cj) in enumerate(zip(ci, cj)):
        springL[idx] = np.linalg.norm(pos[cj, :] - pos[ci, :])
    return springL


def update_for_billowing(points_ini, springL_wing, u_p):
    ##TODO: hardcoded, ideally this would be resolved by the structural solver
    """adjusts the spring lengths of the wing for billowing (hardcoded)
    for the given u_p, for reference visit J.Poland MSc Thesis:
    http://resolver.tudelft.nl/uuid:39d67249-53c9-47b4-84c0-ddac948413a5"""

    def billowing_up(u_p):
        """defines ratio's to multiply the TE lengths by
        from J.Poland MSc Thesis - http://resolver.tudelft.nl/uuid:39d67249-53c9-47b4-84c0-ddac948413a5
        """
        bill_plate_1 = -0.056
        bill_plate_2 = -0.0216
        bill_plate_3 = +0.01
        bill_plate_4 = +0  # outlier, thus 0

        return np.array(
            [
                1 + bill_plate_1 * (1 - u_p),
                1 + bill_plate_2 * (1 - u_p),
                1 + bill_plate_3 * (1 - u_p),
                1 + bill_plate_4 * (1 - u_p),
            ]
        )

    ##TODO: fix hardcoding of indices
    ## Adjusts the spring lengths (multiply by ratios) of the wing for billowing, hardcoded indices
    springL_wing[8] = springL_wing[8] * billowing_up(u_p)[3]  # 17 -> 18
    springL_wing[14] = springL_wing[14] * billowing_up(u_p)[2]  # 16->17
    springL_wing[20] = springL_wing[20] * billowing_up(u_p)[1]  # 15- >16
    springL_wing[26] = springL_wing[26] * billowing_up(u_p)[0]  # 14- >15
    springL_wing[32] = springL_wing[32] * billowing_up(u_p)[1]  # 13 to 14
    springL_wing[38] = springL_wing[38] * billowing_up(u_p)[2]  # 12 to 13
    springL_wing[44] = springL_wing[44] * billowing_up(u_p)[3]  # 11 to 12

    # also adjusts the inner points to make things symmetrical again
    points_ini[14, 1] = points_ini[14, 1] * (
        billowing_up(u_p)[0]
    )  # Move points too to avoid asymmetry
    points_ini[15, 1] = points_ini[15, 1] * (
        billowing_up(u_p)[0]
    )  # Move points too to avoid asymmetry

    return points_ini, springL_wing


def extract_pulley_connectivity(points, bridle_ci, bridle_cj, pulley_data):
    PULLEY_point_indices = pulley_data["point_indices"]
    number_of_pulleys_in_back_lines = pulley_data["number_of_pulleys_in_back_lines"]

    # print(f'PULLEY_point_indices {len(PULLEY_point_indices)}, number_of_pulleys_in_back_lines: {number_of_pulleys_in_back_lines}')
    # Sorting pulleys based on z-coordinate
    sorted_pulley_point_indices = sorted(
        PULLEY_point_indices, key=lambda index: points[index][2]
    )

    # TODO: this is where you input the number of pulleys in the back-lines, that have a different feature
    # for pulleys where the line is below the pulley point
    pulley_point_indices_line_below = sorted_pulley_point_indices[
        :number_of_pulleys_in_back_lines
    ]
    if (
        len(PULLEY_point_indices) > number_of_pulleys_in_back_lines
    ):  # if there are pulleys in the front lines
        pulley_point_indices_line_above = sorted_pulley_point_indices[
            number_of_pulleys_in_back_lines:
        ]
    else:
        pulley_point_indices_line_above = []  # this should be the case for the V3.25

    all_possible_line_indices = []

    for idx, (idx_bridle_node_i, idx_bridle_node_j) in enumerate(
        zip(bridle_ci, bridle_cj)
    ):  # loop through each bridle line
        # if the current line its index i OR j  is a pulley point
        # AND the other line is LOWER than the pulley_point, i.e. pulley line is BELOW the pulley
        if (
            idx_bridle_node_i in pulley_point_indices_line_below
            and points[idx_bridle_node_i][2] > points[idx_bridle_node_j][2]
        ) or (
            idx_bridle_node_j in pulley_point_indices_line_below
            and points[idx_bridle_node_j][2] > points[idx_bridle_node_i][2]
        ):
            all_possible_line_indices.append(idx)

        # if the current line its index i OR j  is a pulley point
        # AND the other line is HIGHER than the pulley_point, i.e. pulley line is ABOVE the pulley
        if (
            idx_bridle_node_i in pulley_point_indices_line_above
            and points[idx_bridle_node_i][2] < points[idx_bridle_node_j][2]
        ) or (
            idx_bridle_node_j in pulley_point_indices_line_above
            and points[idx_bridle_node_j][2] < points[idx_bridle_node_i][2]
        ):
            all_possible_line_indices.append(idx)

    # loop through all_possible_line_indices, twice to try and find the matching line pair
    pulley_line_pair_indices = {}
    pulley_line_indices = []

    for line_index_1 in all_possible_line_indices:
        for line_index_2 in all_possible_line_indices:
            # break if the same line is compared
            if line_index_1 == line_index_2:
                break

            # list the point_indices of the two lines
            point_indices = [
                bridle_ci[line_index_1],
                bridle_cj[line_index_1],
                bridle_ci[line_index_2],
                bridle_cj[line_index_2],
            ]
            # sort the point_indices based on the z-coordinate of the points[index]
            sorted_indices = sorted(point_indices, key=lambda index: points[index][2])

            # IF the pulley line is BELOW the pulley point, i.e. IF LAST two indices correspond to the same point
            # AND the index is in the pulley_index list (i.e. it is actually a pulley)
            # AND the line is not already used, i.e. not in the used_pulley_line_indices set
            if (
                sorted_indices[2] == sorted_indices[3]
                and sorted_indices[3] in pulley_point_indices_line_below
                and line_index_1 not in pulley_line_indices
                and line_index_2 not in pulley_line_indices
            ):
                # Append the line indices, to the pulley_line_indices
                pulley_line_indices.append(line_index_1)
                pulley_line_indices.append(line_index_2)
                # Make new key indices for the pulley line pair indices
                pulley_line_pair_indices[str(line_index_1)] = line_index_2

            # IF the pulley line is ABOVE the pulley point, i.e. IF FIRST indices correspond to the same point
            # AND the index is in the pulley_index list (i.e. it is actually a pulley)
            # AND the line is not already used, i.e. not in the used_pulley_line_indices set
            elif (
                sorted_indices[0] == sorted_indices[1]
                and sorted_indices[0] in pulley_point_indices_line_above
                and line_index_1 not in pulley_line_indices
                and line_index_2 not in pulley_line_indices
            ):
                # Append the line indices, to the pulley_line_indices
                pulley_line_indices.append(line_index_1)
                pulley_line_indices.append(line_index_2)
                # Make new key indices for the pulley line pair indices
                pulley_line_pair_indices[str(line_index_1)] = line_index_2

    # i know that line key and line value make up a pulley
    # so i want to plot both line key and line value
    # extract them separately and append
    pulley_ci_key = [
        bridle_ci[int(key_index)] for key_index in pulley_line_pair_indices.keys()
    ]
    pulley_cj_key = [
        bridle_cj[int(key_index)] for key_index in pulley_line_pair_indices.keys()
    ]
    pulley_ci_value = [
        bridle_ci[int(value_index)] for value_index in pulley_line_pair_indices.values()
    ]
    pulley_cj_value = [
        bridle_cj[int(value_index)] for value_index in pulley_line_pair_indices.values()
    ]
    pulley_ci = np.concatenate((pulley_ci_key, pulley_ci_value))
    pulley_cj = np.concatenate((pulley_cj_key, pulley_cj_value))

    pulley_data["line_indices"] = np.ndarray.flatten(np.array(pulley_line_indices))

    # transform pulley_line_pair indices from a dict to a nested list
    # Why? Because attrs dataclasses can't handle numbers as attributes
    # Does not seem used anywhere
    pulley_line_pair_indices = [
        [key, value] for key, value in pulley_line_pair_indices.items()
    ]

    pulley_data["line_pair_indices"] = pulley_line_pair_indices
    pulley_data["ci"] = pulley_ci
    pulley_data["cj"] = pulley_cj

    additional_dict = {}
    for i in range(0, len(pulley_data["line_indices"]), 2):
        # line 1, key: line 1 and value: line 2 data
        line_1_key = str(pulley_data["line_indices"][i])
        line_1_idx_p3 = bridle_ci[pulley_data["line_indices"][i + 1]]
        line_1_idx_p4 = bridle_cj[pulley_data["line_indices"][i + 1]]
        line_1_rest_length_p3p4 = np.linalg.norm(
            points[line_1_idx_p3] - points[line_1_idx_p4]
        )
        additional_dict[line_1_key] = np.array(
            [line_1_idx_p3, line_1_idx_p4, line_1_rest_length_p3p4]
        )

        # line 2, key: line 2 and value: line 1 data
        line_2_key = str(pulley_data["line_indices"][i + 1])
        line_2_idx_p3 = bridle_ci[pulley_data["line_indices"][i]]
        line_2_idx_p4 = bridle_cj[pulley_data["line_indices"][i]]
        line_2_rest_length_p3p4 = np.linalg.norm(
            points[line_2_idx_p3] - points[line_2_idx_p4]
        )
        additional_dict[line_2_key] = np.array(
            [line_2_idx_p3, line_2_idx_p4, line_2_rest_length_p3p4]
        )

    # transform other_line_pair from a dict to a nested list
    # Why? Because attrs dataclasses can't handle numbers as attributes
    # USED in input_particleSystem for rotational springs
    additional_dict = [[key, value] for key, value in additional_dict.items()]
    pulley_data["other_line_pair"] = additional_dict

    return pulley_data


def new_calculate_mass_distribution(
    points_ini,
    bridle_ci,
    bridle_cj,
    wing_ci,
    wing_cj,
    WING_MASS,
    BRIDLE_RHO,
    BRIDLE_DIAMETER,
    KCU_MASS,
    KCU_INDEX,
    PULLEY_INDICES,
    PULLEY_MASS,
):
    # Define a spring force matrix of the right size
    node_masses = np.zeros(
        points_ini.shape[0]
    )  # Initialising with zero matrix in same shape as points

    ## Bridle lines
    for idx, (idx_bridle_node_i, idx_bridle_node_j) in enumerate(
        zip(bridle_ci, bridle_cj)
    ):  # loop through each bridle line
        # Calculate the length of the bridle line
        length_bridle = np.linalg.norm(
            points_ini[idx_bridle_node_i] - points_ini[idx_bridle_node_j]
        )
        # Calculate the mass of the bridle line
        mass_bridle = BRIDLE_RHO * np.pi * (BRIDLE_DIAMETER / 2) ** 2 * length_bridle
        # Add the mass of the bridle line to the nodes
        node_masses[idx_bridle_node_i] += mass_bridle / 2
        node_masses[idx_bridle_node_j] += mass_bridle / 2

    # print(f'bridle node_masses: {np.sum(node_masses)}')
    ## Pulleys
    for idx in PULLEY_INDICES:
        node_masses[idx] += PULLEY_MASS

    # print(f'pulley& bridle node_masses: {np.sum(node_masses)}')
    ## KCU
    node_masses[KCU_INDEX] += KCU_MASS

    # print(f'pulley& bridle & KCU node_masses: {np.sum(node_masses)}')

    ## Wing
    for idx, (idx_wing_node_i, idx_wing_node_i) in enumerate(
        zip(set(wing_ci), set(wing_cj))
    ):  # making them sets, to have it unique
        node_masses[idx] += WING_MASS / len(set(wing_ci))

    # print(f'pulley& bridle & KCU & wing node_masses: {np.sum(node_masses)}')

    return node_masses


def calculate_projected_area(points):
    # Project points onto the x,y plane
    xy_points = points[:, :2]

    # Find the convex hull
    hull = ConvexHull(xy_points)
    hull_points = xy_points[hull.vertices]

    # Using the shoelace formula
    x = hull_points[:, 0]
    y = hull_points[:, 1]

    return 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))


def read_yaml_file(
    file_path: str,
    yaml_file_name: str,
) -> dict:
    """Reads a yaml file and returns its content as a dictionary

    Args:
        file_path: Path to the yaml file
        yaml_file_name: Name of the yaml file

    Returns:
        dict: Content of the yaml file
    """
    if file_path is None:
        yaml_folder = os.path.dirname(__file__)
        file_path = os.path.join(yaml_folder, yaml_file_name)

    if not os.path.exists(file_path):
        raise FileNotFoundError(f"No such file: {file_path}")

    with open(file_path, "r") as file:
        dict_filled_with_yaml = yaml.load(file, Loader=yaml.SafeLoader)
    return dict_filled_with_yaml


def extract_points_and_connectivity(folder_path_kite_data, surfplan_file):
    """Extracts the points and connectivity of the V3
    This should be done from the surfplan_file, but is instead harcoded
    """

    file_path_points_npy = f"{folder_path_kite_data}/points.npy"
    if os.path.exists(file_path_points_npy):
        points_struc = np.load(file_path_points_npy)
    else:
        raise Exception(f"Error: no file found, with filepath: {file_path_points_npy}.")

    bridle_ci, bridle_cj = extract_bridle_connectivity()
    plate_point_indices = extract_plate_point_indices()
    wing_ci, wing_cj = extract_wing_connectivity(plate_point_indices)
    te_line_indices = extract_te_line_indices(plate_point_indices, wing_ci, wing_cj)
    tube_line_indices = extract_tube_line_indices()

    return (
        np.array(points_struc),
        bridle_ci,
        bridle_cj,
        plate_point_indices,
        wing_ci,
        wing_cj,
        te_line_indices,
        tube_line_indices,
    )


def extract_bridle_connectivity():
    ##TODO: fix hardcoding
    """hardcoded connectivity of the bridle lines"""
    bridle_ci_TE = [
        0,
        21,
        21,
        21,
        22,
        22,
        23,
        23,
        27,
        27,
        24,
        24,
        28,
        28,
        25,
        25,
        29,
        29,
        26,
        26,
        30,
        30,
    ]
    bridle_cj_TE = [
        21,
        22,
        23,
        27,
        24,
        28,
        24,
        1,
        28,
        10,
        25,
        26,
        29,
        30,
        18,
        17,
        11,
        12,
        16,
        15,
        13,
        14,
    ]
    bridle_ci_LE = [
        0,
        0,
        31,
        31,
        34,
        34,
        32,
        32,
        35,
        35,
        33,
        33,
        36,
        36,
        24,
        28,
        31,
        34,
    ]
    bridle_cj_LE = [31, 34, 32, 33, 35, 36, 2, 3, 9, 8, 4, 5, 7, 6, 19, 20, 19, 20]

    bridle_ci = np.append(bridle_ci_TE, bridle_ci_LE)
    bridle_cj = np.append(bridle_cj_TE, bridle_cj_LE)

    return bridle_ci, bridle_cj


def extract_te_line_indices(plate_point_indices, wing_ci, wing_cj):
    """Extracts a list of TE_indices from the plate_point_indices
    using its order: [left_LE, right_LE, right_TE,left_TE]"""
    TE_point_indices = []
    for plate in plate_point_indices:
        TE_point_indices.append(plate[2])
        TE_point_indices.append(plate[3])

    TE_line_indices, tube_line_indices = [], []
    for idx, (ci, cj) in enumerate(zip(wing_ci, wing_cj)):
        if ci in TE_point_indices and cj in TE_point_indices:
            TE_line_indices.append(idx)
    # return [2,8,14,20,26,32,38,44,50] #old hardcoded V3
    return np.array(list(set(TE_line_indices)))


def extract_tube_line_indices():
    ##TODO: fix hardcoding
    """hardcoded V3 indices of the inflatable tubes
    when looping over conn_i, the index that corresponds to a tube"""
    return np.array(
        [
            0,
            1,
            3,
            6,
            7,
            9,
            12,
            13,
            15,
            18,
            19,
            21,
            24,
            25,
            27,
            30,
            31,
            33,
            36,
            37,
            39,
            42,
            43,
            45,
            48,
            49,
            51,
        ]
    )


def extract_plate_point_indices():
    ##TODO: fix hardcoding
    """hardcoded indices of the kite plates"""
    ### Plate connectivity
    plate_1 = [19, 2, 18, 1]
    plate_2 = [2, 3, 17, 18]
    plate_3 = [3, 4, 16, 17]
    plate_4 = [4, 5, 15, 16]
    plate_5 = [5, 6, 14, 15]
    plate_6 = [6, 7, 13, 14]
    plate_7 = [7, 8, 12, 13]
    plate_8 = [8, 9, 11, 12]
    plate_9 = [9, 20, 10, 11]
    plate_point_indices = [
        plate_1,
        plate_2,
        plate_3,
        plate_4,
        plate_5,
        plate_6,
        plate_7,
        plate_8,
        plate_9,
    ]

    return np.array(plate_point_indices)


def extract_wing_connectivity(plate_point_indices):
    """hardcoded connectivity of the wing,
    based on the plate_point_indices"""

    wing_ci, wing_cj = [], []
    for i in np.arange(0, len(plate_point_indices)):
        # The 4 lines describing the tubular frame
        wing_ci.append(plate_point_indices[i][0])  # LE
        wing_cj.append(plate_point_indices[i][1])  # LE

        wing_ci.append(plate_point_indices[i][1])  # Strut right
        wing_cj.append(plate_point_indices[i][2])  # Strut right

        wing_ci.append(plate_point_indices[i][2])  # TE
        wing_cj.append(plate_point_indices[i][3])  # TE

        wing_ci.append(plate_point_indices[i][3])  # Strut left
        wing_cj.append(plate_point_indices[i][0])  # Strut left

        # Them diagonals
        wing_ci.append(plate_point_indices[i][0])
        wing_cj.append(plate_point_indices[i][2])

        wing_ci.append(plate_point_indices[i][1])
        wing_cj.append(plate_point_indices[i][3])

    wing_ci = np.reshape(wing_ci, len(wing_ci))
    wing_cj = np.reshape(wing_cj, len(wing_cj))

    return wing_ci, wing_cj


def update_dict_with_instantiated_classes(
    data_dict: dict, dict_nested_instances: dict
) -> dict:
    """Update the data_dict to include the instantiated classes.

    Args:
        data_dict (dict): The dictionary with the data to be used to create the class.
        dict_nested_instances (dict): The dictionary with the instantiated classes.

    Returns:
        The updated data_dict.
    """
    for key, value in data_dict.items():
        if key in dict_nested_instances:
            data_dict[key] = dict_nested_instances[key]
    return data_dict


def create_attr_class_from_dict(class_name: str, nest0_data_dict: dict):
    """Create an attrs class from a dictionary with up to three levels of nesting.

    Args:
        class_name (str): The name of the class to be created.
        nest0_data_dict (dict): The dictionary with the data to be used to create the class.

    Returns:
        The instantiated class.
    """
    nest0_attributes = {}
    nest1_attributes = {}
    nest2_attributes = {}
    nest3_attributes = {}
    nest4_attributes = {}

    nest1_instances = {}
    nest2_instances = {}
    nest3_instances = {}
    nest4_instances = {}

    # Level 0
    for nest0_key, nest0_value in nest0_data_dict.items():
        if isinstance(nest0_value, dict):
            nest1_data_dict = nest0_value

            # Level 1
            for nest1_key, nest1_value in nest1_data_dict.items():
                if isinstance(nest1_value, dict):
                    nest2_data_dict = nest1_value

                    # Level 2
                    for nest2_key, nest2_value in nest2_data_dict.items():
                        if isinstance(nest2_value, dict):
                            nest3_data_dict = nest2_value

                            # Level 3
                            for nest3_key, nest3_value in nest3_data_dict.items():
                                if isinstance(nest3_value, dict):
                                    nest4_data_dict = nest3_value

                                    # Level 4
                                    for (
                                        nest4_key,
                                        nest4_value,
                                    ) in nest4_data_dict.items():
                                        # find the last level attributes
                                        nest4_attributes[nest4_key] = attr.ib(
                                            default=nest4_value
                                        )

                                    #### LEVEL 4
                                    # 1. Create class,
                                    Nest4Class = attr.make_class(
                                        nest3_key, nest4_attributes, frozen=True
                                    )
                                    # 2. Update the data_dict to include the instantiated classes
                                    # 3. Instantiate the class
                                    nest4_instance = Nest4Class(**nest4_data_dict)
                                    # 4. Add the instantiated class to the dict
                                    nest4_instances[nest3_key] = nest4_instance
                                    # 5. Add to the attributes
                                    nest3_attributes[nest3_key] = attr.ib(
                                        default=nest4_instance
                                    )
                                else:
                                    nest3_attributes[nest3_key] = attr.ib(
                                        default=nest3_value
                                    )

                            ### LEVEL 3
                            # 1. Create class,
                            Nest3Class = attr.make_class(
                                nest2_key, nest3_attributes, frozen=True
                            )
                            # 2. Update the data_dict to include the instantiated classes
                            nest3_data_dict = update_dict_with_instantiated_classes(
                                nest3_data_dict, nest4_instances
                            )
                            # 3. Instantiate the class
                            nest3_instance = Nest3Class(**nest3_data_dict)
                            # 4. Add the instantiated class to the dict
                            nest3_instances[nest2_key] = nest3_instance
                            # 5. Add to the attributes
                            nest2_attributes[nest2_key] = attr.ib(
                                default=nest3_instance
                            )
                        else:
                            nest2_attributes[nest2_key] = attr.ib(default=nest2_value)

                    ### LEVEL 2
                    # 1. Create class,
                    Nest2Class = attr.make_class(
                        nest1_key, nest2_attributes, frozen=True
                    )
                    # 2. Update the data_dict to include the instantiated classes
                    nest2_data_dict = update_dict_with_instantiated_classes(
                        nest2_data_dict, nest3_instances
                    )
                    # 3. Instantiate the class
                    nest2_instance = Nest2Class(**nest2_data_dict)
                    # 4. Add the instantiated class to the dict
                    nest2_instances[nest1_key] = nest2_instance
                    # 5. Add to the attributes
                    nest1_attributes[nest1_key] = attr.ib(default=nest2_instance)
                else:
                    nest1_attributes[nest1_key] = attr.ib(default=nest1_value)

            ### LEVEL 1
            # 1. Create class,
            Nest1Class = attr.make_class(nest0_key, nest1_attributes, frozen=True)
            # 2. Update the data_dict to include the instantiated classes
            nest1_data_dict = update_dict_with_instantiated_classes(
                nest1_data_dict, nest2_instances
            )
            # 3. Instantiate the class
            nest1_instance = Nest1Class(**nest1_data_dict)
            # 4. Add the instantiated class to the dict
            nest1_instances[nest0_key] = nest1_instance
            # Add to the attributes
            nest0_attributes[nest0_key] = attr.ib(default=nest1_instance)
        else:
            nest0_attributes[nest0_key] = attr.ib(default=nest0_value)

    ### LEVEL 0
    # 1. Create class,
    Nest0Class = attr.make_class(class_name, nest0_attributes, frozen=True)
    # 2. Update the data_dict to include the instantiated classes
    nest0_data_dict = update_dict_with_instantiated_classes(
        nest0_data_dict, nest1_instances
    )
    # 3. Instantiate the class
    nest0_instance = Nest0Class(**nest0_data_dict)

    return nest0_instance


def setup_config(
    path_config,
    path_cases,
    case_dir,
    path_kite_config,
    path_processed_data_folder,
    data_kite_dir,
):
    """Setup the configuration for the simulation.

    Steps:
    1. Load the yaml defined settings into dicts
    2. Process kite-specific data and create kite configuration dict
    3. Create a comprehensive config data dictionary for simulation

    Args:
        path_config (str): The path to the configuration file.
        path_processed_data_folder (str): The path to the processed data folder.

    Returns:
        The complete configuration data dictionary.
    """
    # Load the yaml defined settings into dicts
    config_data = read_yaml_file(path_config, None)
    print(f"config_data: {config_data}")
    config_data.update(read_yaml_file(path_cases, "default_case.yaml"))
    print(f"config_data: {config_data}")
    config_data.update(
        read_yaml_file(
            Path(case_dir) / f"{config_data['sim_name']}.yaml",
            f"{config_data['sim_name']}.yaml",
        )
    )
    # path_kite_config = (
    #     Path(path_processed_data_folder)
    #     / config_data["kite_name"]
    #     / f"config_kite_{config_data['kite_name']}.yaml"
    # )
    path_kite_data = (
        Path(path_processed_data_folder)
        / str(config_data["kite_name"])
        / "processed_design_files"
    )
    config_data_kite = read_yaml_file(path_kite_config, None)

    # Process kite-specific YAML data into a structured dictionary
    kite_name = config_data["kite_name"]

    # Load kite-specific module for geometry extraction
    current_dir = os.path.dirname(__file__)
    folder_path_initialisation_kite_specific = os.path.join(current_dir, kite_name)

    # extract_points_and_connectivity = load_module_from_path(
    #     kite_name, f"{folder_path_initialisation_kite_specific}/analyze_surfplan.py"
    # ).extract_points_and_connectivity

    # Parse surfplan file (TODO: implement real surfplan-reader)
    surfplan_file = None

    # Extract geometric data from surfplan
    (
        points_struc,
        bridle_ci,
        bridle_cj,
        plate_point_indices,
        wing_ci,
        wing_cj,
        te_line_indices,
        tube_line_indices,
    ) = extract_points_and_connectivity(data_kite_dir, surfplan_file)

    # Scale points according to geometric scaling factor
    points_ini = np.array(points_struc) / config_data["geometric_scaling_factor"]

    # Calculate initial rest lengths
    bridle_rest_lengths_initial = np.array(
        calculate_edge_lengths(bridle_ci, bridle_cj, points_ini)
    )
    wing_rest_lengths_initial = np.array(
        calculate_edge_lengths(wing_ci, wing_cj, points_ini)
    )  # [m]

    # Apply billowing if enabled
    if config_data["is_billowing_on"]:
        points_ini, wing_rest_lengths = update_for_billowing(
            points_ini, wing_rest_lengths_initial, config_data["u_p"]
        )

    # Number of structural segments
    n_segments = int(len(plate_point_indices))

    # Using the V3_25 data (removed V9_60C as requested)
    # TODO: Most variables here should be generated/imported from the surfplan file instead
    BILLOWING_ANGLES = (
        np.array(config_data_kite["billowing_angles"])
        / config_data["geometric_scaling_factor"]
    )
    TUBE_DIAMETERS = (
        np.array(config_data_kite["tube_diameters"])
        / config_data["geometric_scaling_factor"]
    )
    CANOPY_MAX_HEIGHTS = (
        np.array(config_data_kite["canopy_max_heights"])
        / config_data["geometric_scaling_factor"]
    )

    # Process bridle data
    bridle_data = {}
    bridle_data["diameter"] = config_data_kite["bridle"]["diameter"]
    bridle_data["density"] = config_data_kite["bridle"]["density"]
    bridle_data["bridle_point_index"] = config_data_kite["bridle"]["bridle_point_index"]
    # Element indices
    bridle_data["depower_tape_index"] = config_data_kite["bridle"]["depower_tape_index"]
    bridle_data["left_steering_tape_index"] = config_data_kite["bridle"][
        "left_steering_tape_index"
    ]
    bridle_data["right_steering_tape_index"] = config_data_kite["bridle"][
        "right_steering_tape_index"
    ]

    # Process pulley data
    pulley_data = {}
    pulley_data["point_indices"] = config_data_kite["pulley"]["point_indices"]
    pulley_data["mass"] = config_data_kite["pulley"]["mass"]
    pulley_data["number_of_pulleys_in_back_lines"] = config_data_kite["pulley"][
        "number_of_pulleys_in_back_lines"
    ]

    # Process KCU data
    kcu_data = {}
    kcu_data["index"] = config_data_kite["kcu"]["index"]
    kcu_data["extra"] = {}

    # Handle rotational resistance
    if config_data_kite["is_with_rotational_resistance"]:
        stiffness_bend_strut = config_data_kite["stiffness_bend_strut"]
        stiffness_bend_leading_edge = config_data_kite["stiffness_bend_leading_edge"]
    else:
        stiffness_bend_strut = 0
        stiffness_bend_leading_edge = 0

    # Extract pulley connectivity
    pulley_data = extract_pulley_connectivity(
        points_ini, bridle_ci, bridle_cj, pulley_data
    )

    # Calculate mass distribution
    mass_points = new_calculate_mass_distribution(
        points_ini,
        bridle_ci,
        bridle_cj,
        wing_ci,
        wing_cj,
        config_data_kite["wing_mass"],
        config_data_kite["bridle"]["density"],
        config_data_kite["bridle"]["diameter"],
        config_data_kite["kcu"]["mass"],
        config_data_kite["kcu"]["index"],
        pulley_data["point_indices"],
        config_data_kite["pulley"]["mass"],
    )

    # Calculate reference distances
    # TODO: This is not the true chord, but only the chord taken from the structural discretization
    ref_chord_calculated = max(points_ini[:, 0]) - min(points_ini[:, 0])

    wing_connectivity = np.hstack((wing_ci, wing_cj))
    wing_connectivity = np.unique(wing_connectivity)
    wing_nodes = np.array([points_ini[i] for i in wing_connectivity])
    area_projected_calculated = calculate_projected_area(wing_nodes)
    area_surface_calculated = (
        config_data_kite["area_surface"] / config_data["geometric_scaling_factor"]
    )
    span_calculated = max(wing_nodes[:, 1]) - min(wing_nodes[:, 1])
    height_calculated = max(wing_nodes[:, 2]) - min(wing_nodes[:, 2])

    # Create the kite configuration dictionary
    dict_kite_config = {
        "kite": {
            "points_ini": points_ini,
            "n_points": int(len(points_ini)),
            "surfplan_filename": config_data_kite["surfplan_filename"],
            "area_projected": area_projected_calculated,
            "area_surface": area_surface_calculated,
            "ref_chord": ref_chord_calculated,
            "span": span_calculated,
            "height": height_calculated,
            "wing_mass": config_data_kite["wing_mass"],
            "is_with_elongation_limit": config_data_kite["is_with_elongation_limit"],
            "elongation_limit": config_data_kite["elongation_limit"],
            "is_with_compression_limit": config_data_kite["is_with_compression_limit"],
            "compression_limit": config_data_kite["compression_limit"],
            "limit_stiffness_factor": config_data_kite["limit_stiffness_factor"],
            "billowing_angles": BILLOWING_ANGLES,
            "n_segments": n_segments,
            "wing_rest_lengths_initial": wing_rest_lengths_initial,
            "bridle_rest_lengths_initial": bridle_rest_lengths_initial,
            "mass_points": mass_points,
            "bridle": {
                "diameter": bridle_data["diameter"],
                "density": bridle_data["density"],
                "bridle_point_index": bridle_data["bridle_point_index"],
                "depower_tape_index": bridle_data["depower_tape_index"],
                "left_steering_tape_index": bridle_data["left_steering_tape_index"],
                "right_steering_tape_index": bridle_data["right_steering_tape_index"],
            },
            "pulley": {
                "point_indices": np.array(pulley_data["point_indices"]),
                "mass": np.array(pulley_data["mass"]),
                "number_of_pulleys_in_back_lines": np.array(
                    pulley_data["number_of_pulleys_in_back_lines"]
                ),
                "line_indices": np.array(pulley_data["line_indices"]),
                "line_pair_indices": np.array(pulley_data["line_pair_indices"]),
                "ci": np.array(pulley_data["ci"]),
                "cj": np.array(pulley_data["cj"]),
                "other_line_pair": pulley_data["other_line_pair"],
            },
            "kcu": {
                "drag_coefficient": config_data_kite["kcu"]["drag_coefficient"],
                "diameter": config_data_kite["kcu"]["diameter"],
                "index": config_data_kite["kcu"]["index"],
                "mass": config_data_kite["kcu"]["mass"],
            },
            "connectivity": {
                "bridle_ci": bridle_ci,
                "bridle_cj": bridle_cj,
                "plate_point_indices": plate_point_indices,
                "wing_ci": wing_ci,
                "wing_cj": wing_cj,
                "te_line_indices": te_line_indices,
                "tube_line_indices": tube_line_indices,
            },
            "airfoil": {
                "tube_diameters": TUBE_DIAMETERS,
                "is_tube_diameter_dimensionless": config_data_kite[
                    "is_tube_diameter_dimensionless"
                ],
                "canopy_max_heights": CANOPY_MAX_HEIGHTS,
                "is_canopy_max_height_dimensionless": config_data_kite[
                    "is_canopy_max_height_dimensionless"
                ],
            },
            "stiffness": {
                "bridle": config_data_kite["stiffness_bridle"],
                "tube": config_data_kite["stiffness_tube"],
                "trailing_edge": config_data_kite["stiffness_trailing_edge"],
                "canopy": config_data_kite["stiffness_canopy"],
                # rotational
                "k_bend_strut": stiffness_bend_strut,
                "k_bend_leading_edge": stiffness_bend_leading_edge,
            },
        }
    }

    # Update the config data with the kite configuration
    config_data.update(dict_kite_config)

    # Correcting certain dict items, ensuring the float and np.array types
    config_data["vel_wind"] = np.array([float(i) for i in config_data["vel_wind"]])
    config_data["vel_kite"] = np.array([float(i) for i in config_data["vel_kite"]])
    config_data["acc_kite"] = np.array([float(i) for i in config_data["acc_kite"]])
    config_data["plot_elev"] = np.array([float(i) for i in config_data["plot_elev"]])
    config_data["plot_azim"] = np.array([float(i) for i in config_data["plot_azim"]])

    return create_attr_class_from_dict("config", config_data)


def instantiate_psystem(PROJECT_DIR, kite_name, grav_constant, config, config_kite):

    path_to_config = Path(PROJECT_DIR) / "data" / "config.yaml"
    case_dir = Path(PROJECT_DIR) / "src" / "kitesim" / "cases"
    path_cases = Path(PROJECT_DIR) / "src" / "kitesim" / "cases" / "default_case.yaml"
    path_processed_data_folder = Path(PROJECT_DIR) / "data" / kite_name
    path_kite_config = Path(PROJECT_DIR) / "data" / kite_name / f"config_kite.yaml"
    data_kite_dir = Path(PROJECT_DIR) / "data" / kite_name
    config = setup_config(
        path_to_config,
        path_cases,
        case_dir,
        path_kite_config,
        path_processed_data_folder,
        data_kite_dir,
    )

    ## extracting the connectivity from config_kite
    wing_connectivity = np.column_stack(
        (config_kite["wing_connectivity"]["ci"], config_kite["wing_connectivity"]["cj"])
    )
    bridle_connectivity = np.column_stack(
        (
            config_kite["bridle_connectivity"]["ci"],
            config_kite["bridle_connectivity"]["cj"],
        )
    )
    kite_connectivity = np.vstack((wing_connectivity, bridle_connectivity))

    ## PARAMS
    pss_param_dict = define_param_dict_input_to_pss(
        config, wing_connectivity, kite_connectivity
    )
    pss_param_dict.update({"g": -grav_constant[2]})

    # restructuring connectivity matrix
    pss_kite_connectivity = []
    for idx, _ in enumerate(kite_connectivity):
        if pss_param_dict["is_compression"][idx] and pss_param_dict["is_tension"][idx]:
            linktype = "default"
        elif (
            pss_param_dict["is_compression"][idx]
            and not pss_param_dict["is_tension"][idx]
        ):
            linktype = "nontensile"
        elif pss_param_dict["is_pulley"][idx]:
            linktype = "pulley"
        elif (
            pss_param_dict["is_tension"][idx]
            and not pss_param_dict["is_compression"][idx]
        ):
            linktype = "noncompressive"

        logging.debug(f"idx: {idx}")
        logging.debug(f"kite_connectivity[idx]: {kite_connectivity[idx]}")
        logging.debug(f"pss_param_dict['k'][idx]: {pss_param_dict['k'][idx]}")
        logging.debug(f"pss_param_dict['c']: {pss_param_dict['c']}")
        logging.debug(f"linktype: {linktype}")

        pss_kite_connectivity.append(
            [
                int(kite_connectivity[idx][0]),
                int(kite_connectivity[idx][1]),
                float(pss_param_dict["k"][idx]),
                float(pss_param_dict["c"]),
                string_to_springdampertype(linktype),
            ]
        )

    ## INITIAL CONDITIONS
    points_ini = np.array(config.kite.points_ini)
    for i, point in enumerate(points_ini):
        print(f"[{i}, {point[0]}, {point[1]}, {point[2]}, '<enter-type>']")
    # print(f"points_ini: {points_ini}")

    breakpoint()
    if config.is_with_initial_point_velocity:
        print("Error: initial point velocity has never been defined")
    else:
        vel_ini = np.zeros(points_ini.shape)
    m_array = config.kite.mass_points
    fixed_nodes = np.array(config.kite.bridle.bridle_point_index)
    # fill with: position, initial velocity?, mass, fixed boolean
    pss_initial_conditions = []
    n = config.kite.n_points
    for i in range(n):
        if i in fixed_nodes:
            pss_initial_conditions.append([points_ini[i], vel_ini[i], m_array[i], True])
        else:
            pss_initial_conditions.append(
                [points_ini[i], vel_ini[i], m_array[i], False]
            )

    psystem = ParticleSystem(
        pss_kite_connectivity,
        pss_initial_conditions,
        pss_param_dict,
    )

    if config.is_with_initial_plot:
        from kitesim.post_processing import plotting

        points_between_dict = extract_points_between_dict(config)
        plotting.plot_initial_geometry(config, points_between_dict)

    return psystem, pss_param_dict, config


def run_pss(psystem, params, f_external):
    f_ext = f_external
    t_vector_internal = np.linspace(
        params["dt"], params["t_steps"] * params["dt"], params["t_steps"]
    )
    E_kin = []
    f_int = []
    E_kin_tol = 1e-3  # 1e-29

    logging.debug(f"Running PS simulation, f_int: {psystem.f_int}")

    # And run the simulation
    for step_internal in t_vector_internal:
        psystem.kin_damp_sim(f_ext)

        E_kin.append(np.linalg.norm(psystem.x_v_current[1] ** 2))
        f_int.append(np.linalg.norm(psystem.f_int))

        converged = False
        if step_internal > 10:
            if np.max(E_kin[-10:-1]) <= E_kin_tol:
                converged = True
        if converged and step_internal > 1:
            print("Kinetic damping PS converged", step_internal)
            break
    return psystem
