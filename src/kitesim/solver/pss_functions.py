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
import numpy as np
import sys


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


def define_param_dict_input_to_pss(
    points_ini,
    wing_connectivity,
    kite_connectivity,
    config_dict,
    config_kite_dict,
    rest_lengths,
):
    n_wing_elements = len(wing_connectivity)

    # Transform pulley.other_line_pair to a dict
    # data_struc is [["3",value],["5", value], ...]

    pulley_data_from_config_kite_dict = {
        "point_indices": config_kite_dict["pulley"]["point_indices"],
        "mass": config_kite_dict["pulley"]["mass"],
        "number_of_pulleys_in_back_lines": config_kite_dict["pulley"][
            "number_of_pulleys_in_back_lines"
        ],
    }

    pulley_data = extract_pulley_connectivity(
        points_ini,
        config_kite_dict["bridle_connectivity"]["ci"],
        config_kite_dict["bridle_connectivity"]["cj"],
        pulley_data_from_config_kite_dict,
    )

    other_line_pair_dict = {}
    for entry in pulley_data["other_line_pair"]:
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

            if i in config_kite_dict["te_line_indices"]:
                stiffness_list.append(config_kite_dict["stiffness_trailing_edge"])
                is_compression_list.append(False)
                is_rotational_list.append(False)
            elif i in config_kite_dict["tube_line_indices"]:
                stiffness_list.append(config_kite_dict["stiffness_tube"])
                is_compression_list.append(True)
                is_rotational_list.append(True)
            else:  # must be a canopy element
                stiffness_list.append(config_kite_dict["stiffness_canopy"])
                is_compression_list.append(False)
                is_rotational_list.append(False)

        # if bridle-lines
        else:
            is_compression_list.append(False)
            is_tension_list.append(True)
            is_rotational_list.append(True)
            stiffness_list.append(config_kite_dict["stiffness_bridle"])

            # TODO: might be better to use one index style/structure?
            # the (i - n_wing_elements) is to correct the indices
            if (i - n_wing_elements) in config_kite_dict["pulley_line_indices"]:
                # print(f'pulley, i: {i}')
                # print(f'connection[i]: {kite_connectivity[i]}')
                is_pulley_list.append(True)
            else:
                is_pulley_list.append(False)

    params = {
        "pulley_other_line_pair": other_line_pair_corrected_dict,
        "k": np.array(stiffness_list),
        "is_compression": np.array(is_compression_list),
        "is_tension": np.array(is_tension_list),
        "is_pulley": np.array(is_pulley_list),
        "is_rotational": np.array(is_rotational_list),
        "l0": rest_lengths,
    }
    params.update(
        {
            "c": config_dict["solver"]["damping_constant"],
            "dt": config_dict["solver"]["dt"],
            "t_steps": config_dict["solver"]["n_time_steps"],
            "abs_tol": config_dict["solver"]["abs_tol"],
            "rel_tol": config_dict["solver"]["rel_tol"],
            "max_iter": config_dict["solver"]["max_iter"],
            "n": len(config_kite_dict["points"]),
            "aerostructural_tol": config_dict["aero_structural"]["tol"],
            "is_with_visc_damping": config_dict["solver"]["is_with_visc_damping"],
        }
    )

    return params


def instantiate_psystem(
    config_dict,
    config_kite_dict,
    points_ini,
    wing_connectivity,
    kite_connectivity,
    rest_lengths,
    m_array,
):
    ## PARAMS
    pss_param_dict = define_param_dict_input_to_pss(
        points_ini,
        wing_connectivity,
        kite_connectivity,
        config_dict,
        config_kite_dict,
        rest_lengths,
    )
    pss_param_dict.update({"g": -config_dict["grav_constant"][2]})

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
                SpringDamperType(linktype.lower()),
            ]
        )

    ## INITIAL CONDITIONS
    if config_dict["is_with_initial_point_velocity"]:
        print("Error: initial point velocity has never been defined")
    else:
        vel_ini = np.zeros(points_ini.shape)

    fixed_nodes = np.array(config_kite_dict["bridle"]["bridle_point_index"])
    pss_initial_conditions = []
    n = len(points_ini)
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
    return psystem, pss_param_dict, pss_kite_connectivity


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
            # print("Kinetic damping PS converged", step_internal)
            break
    return psystem
