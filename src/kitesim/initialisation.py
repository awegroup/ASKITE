import numpy as np


def initialize_wing_structure(geometry_dict):
    """
    Create the structural nodes and connectivity for the kite structure.

    Args:
        geometry_dict (dict): Kite configuration dictionary.

    Returns:
        tuple: (struc_nodes, wing_ci, wing_cj, bridle_ci, bridle_cj, struc_node_le_indices,
                struc_node_te_indices, pulley_point_indices, tubular_frame_line_idx_list,
                te_line_idx_list, n_struc_ribs)
    """
    n_wing_nodes = len(geometry_dict["wing_nodes"]["data"])
    wing_ci = []
    wing_cj = []
    tubular_frame_line_idx_list = []
    te_line_idx_list = []
    struc_node_le_indices = []
    struc_node_te_indices = []
    struc_nodes_list = [[0, 0, 0]]
    pulley_point_indices = []
    node_masses_wing = np.zeros(n_wing_nodes + 1)
    for _, (node_id, x, y, z, node_type) in enumerate(
        geometry_dict["wing_nodes"]["data"]
    ):
        # TODO: Improve distribution to preserve inertia
        # Distributing wing mass uniformly
        node_masses_wing[node_id] = geometry_dict["wing_mass"] / n_wing_nodes

        struc_nodes_list.append([x, y, z])
        if node_type == "pulley":
            pulley_point_indices.append(node_id)

        # if uneven idx
        if node_id % 2 != 0:
            curr_le = node_id
            curr_te = node_id + 1
            nxt_le = node_id + 2
            nxt_te = node_id + 3
            struc_node_le_indices.append(curr_le)
            struc_node_te_indices.append(curr_te)
            # if not the last point
            if node_id < (n_wing_nodes - 1):
                # --- For first to up to second to last panel only connect to the next ---
                wing_ci.append(curr_le)
                wing_cj.append(curr_te)
                tubular_frame_line_idx_list.append(len(wing_ci) - 1)
                # connect to next panel’s LE (LE[i] -- LE[i+1])
                wing_ci.append(curr_le)
                wing_cj.append(nxt_le)
                tubular_frame_line_idx_list.append(len(wing_ci) - 1)
                # connect to next panel’s TE (TE[i] -- TE[i+1])
                wing_ci.append(curr_te)
                wing_cj.append(nxt_te)
                te_line_idx_list.append(len(wing_ci) - 1)
            else:
                # --- last panel only needs to connect to itself ---
                # connect LE->TE
                wing_ci.append(curr_le)
                wing_cj.append(curr_te)
                tubular_frame_line_idx_list.append(len(wing_ci) - 1)

    wing_connectivity = np.column_stack(
        (
            wing_ci,
            wing_cj,
        )
    )
    wing_rest_lengths_initial = np.zeros(np.array(wing_ci).shape)
    for idx, (ci, cj) in enumerate(zip(wing_ci, wing_cj)):
        wing_rest_lengths_initial[idx] = np.linalg.norm(
            np.array(struc_nodes_list)[cj, :] - np.array(struc_nodes_list)[ci, :]
        )

    return (
        struc_nodes_list,
        wing_ci,
        wing_cj,
        wing_connectivity,
        wing_rest_lengths_initial,
        struc_node_le_indices,
        struc_node_te_indices,
        tubular_frame_line_idx_list,
        te_line_idx_list,
        int(n_wing_nodes / 2),
        pulley_point_indices,
        node_masses_wing,
    )


def initialize_bridle_line_system(
    struc_nodes_list,
    node_masses_wing,
    n_wing_lines,
    pulley_point_indices,
    geometry_dict,
):
    """
    Initialize the bridle line system for the kite.

    Returns:
        tuple: (bridle_ci, bridle_cj, pulley_point_indices)
    """

    bridle_nodes_data = geometry_dict["bridle_nodes"]["data"]
    for idx, (node_id, x, y, z, node_type) in enumerate(bridle_nodes_data):
        struc_nodes_list.append([x, y, z])
        if node_type == "pulley":
            pulley_point_indices.append(node_id)

    bridle_lines = geometry_dict["bridle_lines"]
    headers = bridle_lines["headers"]  # e.g. ["LE_x", "LE_y", ...]
    data = bridle_lines["data"]  # list of lists
    bridle_lines_dict = {}
    for row in data:
        # Use first index of each data row as key, and the rest as value
        bridle_lines_dict[row[0]] = np.array(row[1:])

    node_masses = np.concatenate((node_masses_wing, np.zeros(len(bridle_nodes_data))))
    bridle_ci = []
    bridle_cj = []
    pulley_line_indices = []
    pulley_line_to_other_node_pair_dict = {}
    bridle_rest_lengths_initial = []
    steering_tape_indices = []
    for idx, bridle_data_i in enumerate(geometry_dict["bridle_connections"]["data"]):
        line_name = bridle_data_i[0]
        ci = int(bridle_data_i[1])
        cj = int(bridle_data_i[2])
        bridle_ci.append(ci)
        bridle_cj.append(cj)

        # TODO: if rest lengths are off, then you are currently not adding the correct mass
        # add mass
        line_density = float(bridle_lines_dict[line_name][3])  # density in kg/m^3
        line_diameter = float(bridle_lines_dict[line_name][1])  # diameter in m
        line_length = np.linalg.norm(
            np.array(struc_nodes_list[ci]) - np.array(struc_nodes_list[cj])
        )
        line_mass = line_density * np.pi * (line_diameter / 2) ** 2 * line_length
        node_masses[ci] += line_mass / 2
        node_masses[cj] += line_mass / 2

        if line_name == "Power Tape":
            power_tape_index = idx + n_wing_lines

        if line_name == "Steering Tape":
            steering_tape_indices.append(idx + n_wing_lines)

        if len(bridle_data_i) > 3:
            # if the connection ck exists, i.e. if there is a third connection, this line is a pulley
            pulley_line_index = idx + n_wing_lines
            pulley_line_indices.append(pulley_line_index)

            ck = int(bridle_data_i[3])
            if ci in pulley_point_indices:
                pulley_index = ci
            else:
                pulley_index = cj

            # Computing line length using the current point-to-point distance
            # and using it to divide the total pulley_length between both line segments
            point_to_point_line_len_current = np.linalg.norm(
                np.array(struc_nodes_list[ci]) - np.array(struc_nodes_list[cj])
            )
            point_to_point_line_len_other = np.linalg.norm(
                np.array(struc_nodes_list[pulley_index])
                - np.array(struc_nodes_list[ck])
            )
            total_point_to_point_line_len = (
                point_to_point_line_len_current + point_to_point_line_len_other
            )
            total_pulley_len = float(bridle_lines_dict[line_name][0])
            line_len_current = (
                point_to_point_line_len_current / total_point_to_point_line_len
            ) * total_pulley_len
            line_len_other = (
                point_to_point_line_len_other / total_point_to_point_line_len
            ) * total_pulley_len

            # Create a special mapping for the Structural Particle System Solver
            # key: pulley_line_index
            # value: [pulley_index, ck, line_len_other]
            # This is used to connect the pulley line to the other node pair
            pulley_line_to_other_node_pair_dict[str(pulley_line_index)] = np.array(
                [
                    pulley_index,
                    ck,
                    line_len_other,
                ]
            )
            bridle_rest_lengths_initial.append(line_len_current)

        else:
            bridle_rest_lengths_initial.append(float(bridle_lines_dict[line_name][0]))

    for idx in pulley_point_indices:
        # Add pulley mass to the node mass
        node_masses[idx] += geometry_dict["pulley_mass"]

    bridle_connectivity = np.column_stack(
        (
            bridle_ci,
            bridle_cj,
        )
    )

    return (
        np.array(struc_nodes_list),
        bridle_ci,
        bridle_cj,
        bridle_connectivity,
        bridle_rest_lengths_initial,
        pulley_point_indices,
        pulley_line_indices,
        pulley_line_to_other_node_pair_dict,
        power_tape_index,
        steering_tape_indices,
        node_masses,
    )


def main(geometry_dict):
    """
    Main entry point for kite structural initialisation.

    Args:
        geometry_dict (dict): Kite configuration dictionary.

    Returns:
        tuple: All structural and connectivity arrays, mass arrays, rest lengths, and pulley data.
    """

    (
        struc_nodes,
        wing_ci,
        wing_cj,
        wing_connectivity,
        wing_rest_lengths_initial,
        struc_node_le_indices,
        struc_node_te_indices,
        tubular_frame_line_idx_list,
        te_line_idx_list,
        n_struc_ribs,
        pulley_point_indices,
        node_masses_wing,
    ) = initialize_wing_structure(geometry_dict)

    (
        struc_nodes,
        bridle_ci,
        bridle_cj,
        bridle_connectivity,
        bridle_rest_lengths_initial,
        pulley_point_indices,
        pulley_line_indices,
        pulley_line_to_other_node_pair_dict,
        power_tape_index,
        steering_tape_indices,
        node_masses,
    ) = initialize_bridle_line_system(
        struc_nodes,
        node_masses_wing,
        len(wing_connectivity),
        pulley_point_indices,
        geometry_dict,
    )

    kite_connectivity = np.vstack((wing_connectivity, bridle_connectivity))
    rest_lengths = np.concatenate(
        (wing_rest_lengths_initial, bridle_rest_lengths_initial)
    )

    return (
        struc_nodes,
        wing_ci,
        wing_cj,
        bridle_ci,
        bridle_cj,
        struc_node_le_indices,
        struc_node_te_indices,
        pulley_point_indices,
        tubular_frame_line_idx_list,
        te_line_idx_list,
        n_struc_ribs,
        wing_connectivity,
        bridle_connectivity,
        kite_connectivity,
        node_masses,
        bridle_rest_lengths_initial,
        wing_rest_lengths_initial,
        rest_lengths,
        pulley_point_indices,
        pulley_line_indices,
        pulley_line_to_other_node_pair_dict,
        power_tape_index,
        steering_tape_indices,
    )
