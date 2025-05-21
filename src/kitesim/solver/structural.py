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


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def plot_3d_kite_structure(points, connectivity, fixed_nodes=None, pulley_nodes=None):
    """
    Plot the 3D structure of a kite with enhanced visualization features.

    Parameters:
    - points: array of 3D coordinates for each node
    - connectivity: list of [i, j, k, c, type] where i, j are node indices,
                   k is stiffness, c is damping, and type is the spring-damper type
    - fixed_nodes: indices of fixed nodes (optional)
    - pulley_nodes: indices of pulley nodes (optional)
    """
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection="3d")

    # Create sets for fixed and pulley nodes if not provided
    if fixed_nodes is None:
        fixed_nodes = set()
    else:
        fixed_nodes = set(np.atleast_1d(fixed_nodes))

    if pulley_nodes is None:
        pulley_nodes = set()
    else:
        pulley_nodes = set(np.atleast_1d(pulley_nodes))

    # Extract node masses from connectivity (using the m_array would be better if available)
    node_masses = {}
    for conn in connectivity:
        i, j = int(conn[0]), int(conn[1])
        if hasattr(conn[4], "value"):
            link_type = conn[4].value
        else:
            link_type = conn[4]

        # Initialize masses if not already in dictionary
        if i not in node_masses:
            node_masses[i] = 0
        if j not in node_masses:
            node_masses[j] = 0

    # Create sets to track which elements are tubular frame, te_lines, or other noncompressive
    tubular_frame_nodes = set()
    te_line_nodes = set()
    bridle_line_nodes = set()
    pulley_line_nodes = set()

    # Line style mapping
    line_styles = {
        "default": {
            "color": "black",
            "linestyle": "-",
            "linewidth": 2.5,
            "label": "Tubular Frame",
        },
        "noncompressive": {"color": "green", "linestyle": "-", "linewidth": 1.5},
        "pulley": {
            "color": "purple",
            "linestyle": "-",
            "linewidth": 1.5,
            "label": "Pulley Lines",
        },
    }

    # Track which labels have been used
    used_labels = set()

    # First pass to identify TE lines and bridle lines
    # This is necessary because we need to know which noncompressive lines are TE lines before plotting
    for conn in connectivity:
        i, j = int(conn[0]), int(conn[1])

        if hasattr(conn[4], "value"):
            link_type = conn[4].value
        else:
            link_type = conn[4]

        # Mark te_line_idx_list nodes (this is a placeholder - in actual code,
        # we would use the te_line_idx_list parameter to identify TE lines)
        # For now, we're just propagating the te_line_nodes set
        if link_type.lower() == "noncompressive":
            if i in te_line_nodes or j in te_line_nodes:
                te_line_nodes.add(i)
                te_line_nodes.add(j)

    # Plot connections with appropriate styling
    for conn in connectivity:
        i, j = int(conn[0]), int(conn[1])
        k, c = float(conn[2]), float(conn[3])

        if hasattr(conn[4], "value"):
            link_type = conn[4].value
        else:
            link_type = conn[4]

        x_vals = [points[i][0], points[j][0]]
        y_vals = [points[i][1], points[j][1]]
        z_vals = [points[i][2], points[j][2]]

        # Default styling
        style = line_styles.get(
            link_type.lower(), {"color": "gray", "linestyle": "-", "linewidth": 1}
        )

        # Separate noncompressive elements into TE lines and bridle lines
        if link_type.lower() == "noncompressive":
            if i in te_line_nodes or j in te_line_nodes:
                style["color"] = "orange"
                if "Canopy TE" not in used_labels:
                    style["label"] = "Canopy TE"
                    used_labels.add("Canopy TE")
                else:
                    style.pop("label", None)

            else:
                style["color"] = "blue"
                if "Bridle Lines" not in used_labels:
                    style["label"] = "Bridle Lines"
                    used_labels.add("Bridle Lines")
                else:
                    style.pop("label", None)

        # Track nodes for tubular frame and pulley lines
        if link_type.lower() == "default":
            tubular_frame_nodes.add(i)
            tubular_frame_nodes.add(j)
            if "Tubular Frame" not in used_labels:
                used_labels.add("Tubular Frame")
            else:
                style.pop("label", None)

        if link_type.lower() == "pulley":
            pulley_line_nodes.add(i)
            pulley_line_nodes.add(j)
            if "Pulley Lines" not in used_labels:
                used_labels.add("Pulley Lines")
            else:
                style.pop("label", None)

        # Include damping in the label if requested
        if "label" in style and "damping" not in style["label"]:
            style["label"] += f" (k={k:.1f}, c={c:.2f})"

        # Plot the line
        ax.plot(x_vals, y_vals, z_vals, **style)

    # Create legend labels for nodes
    node_handles = []
    node_labels = []

    # Plot nodes - separate loop to ensure nodes are drawn on top of lines
    for i, point in enumerate(points):
        # Plot the index of the node
        ax.text(
            point[0] + 0.02,
            point[1] + 0.02,
            point[2] + 0.02,
            str(i),
            color="black",
            fontsize=6,
        )
        if i in fixed_nodes:
            marker = ax.scatter(
                point[0],
                point[1],
                point[2],
                color="red",
                s=5,
                label="",  # We'll add to legend separately
            )
            if "Fixed Node" not in used_labels:
                node_handles.append(marker)
                node_labels.append("Fixed Node")
                used_labels.add("Fixed Node")
        elif i in pulley_nodes:
            marker = ax.scatter(
                point[0],
                point[1],
                point[2],
                color="purple",
                s=8,
                label="",  # We'll add to legend separately
            )
            if "Pulley Node" not in used_labels:
                node_handles.append(marker)
                node_labels.append("Pulley Node")
                used_labels.add("Pulley Node")
        else:
            marker = ax.scatter(
                point[0],
                point[1],
                point[2],
                color="black",
                s=8,
                label="",  # We'll add to legend separately
            )
            if "Free Node" not in used_labels:
                node_handles.append(marker)
                node_labels.append("Free Node")
                used_labels.add("Free Node")

    for idx, (i, j, k, _, line_type) in enumerate(connectivity):
        # Get coordinates
        p1 = np.array(points[i])
        p2 = np.array(points[j])

        # Midpoint for label
        midpoint = (p1 + p2) / 2
        # label = f"{line_type.name}\nk={k:.1e}"
        label = f"{idx}"

        # Add label slightly offset from midpoint
        offset = 0.02 * np.linalg.norm(p2 - p1)
        ax.text(
            midpoint[0] + offset,
            midpoint[1] + offset,
            midpoint[2] + offset,
            label,
            fontsize=6,
            color="blue",
        )

    # Set labels and title
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title("3D Kite Structure")

    # Equal aspect ratio
    if hasattr(points, "max") and hasattr(points, "min"):
        bb = points.max(axis=0) - points.min(axis=0)
        ax.set_box_aspect(bb)
    else:
        # If points is not a numpy array with max/min methods
        points_array = np.array(points)
        bb = points_array.max(axis=0) - points_array.min(axis=0)
        ax.set_box_aspect(bb)

    # Add legend - use a separate legend for nodes
    # Get existing handles and labels from the lines
    handles, labels = ax.get_legend_handles_labels()

    # Combine with node handles and labels
    all_handles = handles + node_handles
    all_labels = labels + node_labels

    # Create legend outside the plot area to ensure visibility
    plt.legend(
        all_handles,
        all_labels,
        loc="upper left",
        bbox_to_anchor=(1.05, 1),
        borderaxespad=0,
    )

    # Adjust layout to make room for the legend
    plt.tight_layout(rect=[0, 0, 0.85, 1])  # Leave space on the right for the legend

    plt.show()


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


def instantiate_psystem(
    config_dict,
    config_kite_dict,
    points_ini,
    wing_connectivity,
    kite_connectivity,
    rest_lengths,
    m_array,
    tubular_frame_line_idx_list,
    te_line_idx_list,
    pulley_point_indices,
):

    ## PARAMS
    n_wing_elements = len(wing_connectivity)

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

    pulley_data["line_indices"] = [
        idx + n_wing_elements for idx in pulley_data["line_indices"]
    ]

    pss_param_dict = {
        "pulley_other_line_pair": other_line_pair_corrected_dict,
        "l0": rest_lengths,
        "c": config_dict["solver"]["damping_constant"],
        "dt": config_dict["solver"]["dt"],
        "t_steps": config_dict["solver"]["n_time_steps"],
        "abs_tol": config_dict["solver"]["abs_tol"],
        "rel_tol": config_dict["solver"]["rel_tol"],
        "max_iter": config_dict["solver"]["max_iter"],
        "n": len(config_kite_dict["points"]),
        "aerostructural_tol": config_dict["aero_structural"]["tol"],
        "is_with_visc_damping": config_dict["solver"]["is_with_visc_damping"],
        "g": -config_dict["grav_constant"][2],
    }

    # creating PSS style connectivity matrix
    pss_kite_connectivity = []
    for idx, _ in enumerate(kite_connectivity):
        if idx in tubular_frame_line_idx_list:
            # compression and tension
            k = config_kite_dict["tubular_frame_stiffness"]
            c = config_kite_dict["tubular_frame_damping"]
            linktype = "default"
        elif idx in te_line_idx_list:
            # only tension
            k = config_kite_dict["trailing_edge_stiffness"]
            c = config_kite_dict["trailing_edge_damping"]
            linktype = "noncompressive"
        elif idx in pulley_data["line_indices"]:
            # pulley
            k = config_kite_dict["stiffness_bridle"]
            c = config_kite_dict["bridle_line_damping"]
            linktype = "pulley"
        else:
            # only compression
            k = config_kite_dict["stiffness_bridle"]
            c = config_kite_dict["canopy_damping"]
            linktype = "noncompressive"

        pss_kite_connectivity.append(
            [
                int(kite_connectivity[idx][0]),
                int(kite_connectivity[idx][1]),
                float(k),
                float(c),
                SpringDamperType(linktype.lower()),
            ]
        )

    ## INITIAL CONDITIONS
    if config_dict["is_with_initial_point_velocity"]:
        raise ValueError("Error: initial point velocity has never been defined")
    else:
        vel_ini = np.zeros((len(points_ini), 3))

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

    if config_dict["is_with_initial_structure_plot"]:
        plot_3d_kite_structure(
            points_ini,
            pss_kite_connectivity,
            fixed_nodes=fixed_nodes,
            pulley_nodes=pulley_point_indices,
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
