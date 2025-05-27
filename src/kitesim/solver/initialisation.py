import numpy as np
import pandas as pd
from collections import defaultdict


def calculate_edge_lengths(ci, cj, pos):
    """
    Calculate edge lengths between nodes for given connectivity.

    Args:
        ci (np.ndarray): Start node indices.
        cj (np.ndarray): End node indices.
        pos (np.ndarray): Node positions (n_nodes,3).

    Returns:
        np.ndarray: Edge lengths for each connection.
    """
    springL = np.zeros(ci.shape)
    for idx, (ci, cj) in enumerate(zip(ci, cj)):
        springL[idx] = np.linalg.norm(pos[cj, :] - pos[ci, :])
    return springL


def initialize_wing_structure(config_kite_dict):
    """
    Create the structural nodes and connectivity for the kite structure.

    Args:
        config_kite_dict (dict): Kite configuration dictionary.

    Returns:
        tuple: (struc_nodes, wing_ci, wing_cj, bridle_ci, bridle_cj, struc_le_idx_list,
                struc_te_idx_list, pulley_point_indices, tubular_frame_line_idx_list,
                te_line_idx_list, n_struc_ribs)
    """

    # # Extract airfoil
    # af = config_kite_dict["airfoils"]
    # headers = af["headers"]  # e.g. ["LE_x", "LE_y", ...]
    # data = af["data"]  # list of lists
    # airfoil_data = config_kite_dict["airfoils"]["data"]
    # df_airfoil = pd.DataFrame(data, columns=headers)

    # wing_ci = []
    # wing_cj = []

    # # 1) Pick out exactly the panel‐indices you’ll represent in the structure:
    # n = len(df_airfoil)
    # selected = [
    #     i for i in range(n) if i == 0 or df_airfoil.loc[i, "is_strut"] or i == n - 1
    # ]

    # # 2) Build your structural‐node list and remember for each panel which
    # struc_nodes_list = [[0, 0, 0]]
    # #    structural index is its LE and which is its TE
    # le_map = {}
    # te_map = {}
    # struc_le_idx_list = []
    # struc_te_idx_list = []
    # for i in selected:
    #     row = df_airfoil.loc[i]
    #     # append LE
    #     le_idx = len(struc_nodes_list)
    #     struc_nodes_list.append(row[["LE_x", "LE_y", "LE_z"]].tolist())
    #     # append TE
    #     te_idx = len(struc_nodes_list)
    #     struc_nodes_list.append(row[["TE_x", "TE_y", "TE_z"]].tolist())
    #     le_map[i] = le_idx
    #     te_map[i] = te_idx
    #     # append to the list of indices
    #     struc_le_idx_list.append(le_idx)
    #     struc_te_idx_list.append(te_idx)

    # tubular_frame_line_idx_list = []
    # te_line_idx_list = []
    # # 3) Now build wing_ci/wing_cj with explicit if‐branches
    # for pos, i in enumerate(selected):
    #     curr_le = le_map[i]
    #     curr_te = te_map[i]

    #     if pos != len(selected) - 1:
    #         # --- For first to up to second to last panel only connect to the next ---
    #         wing_ci.append(curr_le)
    #         wing_cj.append(curr_te)
    #         tubular_frame_line_idx_list.append(len(wing_ci) - 1)
    #         # connect to next panel’s LE (LE[i] -- LE[i+1])
    #         nxt = selected[pos + 1]
    #         wing_ci.append(curr_le)
    #         wing_cj.append(le_map[nxt])
    #         tubular_frame_line_idx_list.append(len(wing_ci) - 1)
    #         # connect to next panel’s TE (TE[i] -- TE[i+1])
    #         wing_ci.append(curr_te)
    #         wing_cj.append(te_map[nxt])
    #         te_line_idx_list.append(len(wing_ci) - 1)
    #     elif pos == len(selected) - 1:
    #         # --- last panel only needs to connect to itself ---
    #         # connect LE->TE
    #         wing_ci.append(curr_le)
    #         wing_cj.append(curr_te)
    #         tubular_frame_line_idx_list.append(len(wing_ci) - 1)

    # for i, struc_node in enumerate(struc_nodes_list):
    #     if i == 0:
    #         continue
    #     struc_nodes_list[i][0] = struc_node[0] + 0.7
    #     struc_nodes_list[i][2] = struc_node[2] + 7.3

    ## new attempt
    n_wing_nodes = len(config_kite_dict["wing_nodes"]["data"])
    wing_ci = []
    wing_cj = []
    tubular_frame_line_idx_list = []
    te_line_idx_list = []
    struc_le_idx_list = []
    struc_te_idx_list = []
    struc_nodes_list = [[0, 0, 0]]
    pulley_point_indices = []
    node_masses_wing = np.zeros(n_wing_nodes + 1)
    for _, (node_id, x, y, z, node_type) in enumerate(
        config_kite_dict["wing_nodes"]["data"]
    ):
        # TODO: Improve distribution to preserve inertia
        # Distributing wing mass uniformly
        node_masses_wing[node_id] = config_kite_dict["wing_mass"] / n_wing_nodes

        struc_nodes_list.append([x, y, z])
        if node_type == "pulley":
            pulley_point_indices.append(node_id)

        # if uneven idx
        if node_id % 2 != 0:
            curr_le = node_id
            curr_te = node_id + 1
            nxt_le = node_id + 2
            nxt_te = node_id + 3
            struc_le_idx_list.append(curr_le)
            struc_te_idx_list.append(curr_te)
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

    # if config_kite_dict["is_with_bridle_line_system_surfplan"]:
    #     df_bridle_line_system_surfplan = pd.DataFrame(
    #         config_kite_dict["df_bridle_line_system_surfplan"]["data"],
    #         columns=config_kite_dict["df_bridle_line_system_surfplan"]["headers"],
    #     )

    #     # Initialize connectivity lists
    #     bridle_ci = []
    #     bridle_cj = []

    #     # Helper function to compare 3D points (rounded to tolerance)
    #     def point_key(p, tol=1e-8):
    #         return tuple(
    #             np.round(p, decimals=8)
    #         )  # prevent floating point uniqueness issues

    #     # Build a dictionary of existing nodes for fast lookup
    #     existing_node_map = {
    #         point_key(node): idx for idx, node in enumerate(struc_nodes_list)
    #     }

    #     # Iterate over each bridle line
    #     for _, row in df_bridle_line_system_surfplan.iterrows():
    #         # Extract p1 and p2
    #         p1 = [row["p1_x"], row["p1_y"], row["p1_z"]]
    #         p2 = [row["p2_x"], row["p2_y"], row["p2_z"]]

    #         # Deduplicate and register p1
    #         key1 = point_key(p1)
    #         if key1 in existing_node_map:
    #             idx1 = existing_node_map[key1]
    #         else:
    #             idx1 = len(struc_nodes_list)
    #             struc_nodes_list.append(p1)
    #             existing_node_map[key1] = idx1

    #         # Deduplicate and register p2
    #         key2 = point_key(p2)
    #         if key2 in existing_node_map:
    #             idx2 = existing_node_map[key2]
    #         else:
    #             idx2 = len(struc_nodes_list)
    #             struc_nodes_list.append(p2)
    #             existing_node_map[key2] = idx2

    #         # Store connectivity
    #         bridle_ci.append(idx1)
    #         bridle_cj.append(idx2)
    # else:
    #     # Extract brilde_nodes
    #     bridle_nodes = config_kite_dict["bridle_nodes"]
    #     headers = bridle_nodes["headers"]  # e.g. ["LE_x", "LE_y", ...]
    #     data = bridle_nodes["data"]  # list of lists
    #     bridle_nodes_data = config_kite_dict["bridle_nodes"]["data"]
    #     df_bridle_nodes = pd.DataFrame(data, columns=headers)

    #     pulley_point_indices = []
    #     for _, row in df_bridle_nodes.iterrows():
    #         struc_nodes_list.append(row[["x", "y", "z"]].tolist())
    #         if row["type"] == "pulley":
    #             pulley_point_indices.append(len(struc_nodes_list) - 1)

    #     bridle_ci = config_kite_dict["bridle_lines"]["ci"]
    #     bridle_cj = config_kite_dict["bridle_lines"]["cj"]

    #     # Given that we know
    wing_connectivity = np.column_stack(
        (
            wing_ci,
            wing_cj,
        )
    )
    wing_rest_lengths_initial = np.array(
        calculate_edge_lengths(
            np.array(wing_ci),
            np.array(wing_cj),
            np.array(struc_nodes_list),
        )
    )
    return (
        struc_nodes_list,
        wing_ci,
        wing_cj,
        wing_connectivity,
        wing_rest_lengths_initial,
        struc_le_idx_list,
        struc_te_idx_list,
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
    config_kite_dict,
):
    """
    Initialize the bridle line system for the kite.

    Returns:
        tuple: (bridle_ci, bridle_cj, pulley_point_indices)
    """
    # if config_kite_dict["is_with_bridle_line_system_surfplan"]:
    #     df_bridle_line_system_surfplan = pd.DataFrame(
    #         config_kite_dict["df_bridle_line_system_surfplan"]["data"],
    #         columns=config_kite_dict["df_bridle_line_system_surfplan"]["headers"],
    #     )

    #     # Initialize connectivity lists
    #     bridle_ci = []
    #     bridle_cj = []

    #     # Helper function to compare 3D points (rounded to tolerance)
    #     def point_key(p, tol=1e-8):
    #         return tuple(
    #             np.round(p, decimals=8)
    #         )  # prevent floating point uniqueness issues

    #     # Build a dictionary of existing nodes for fast lookup
    #     existing_node_map = {
    #         point_key(node): idx for idx, node in enumerate(struc_nodes_list)
    #     }

    #     # Iterate over each bridle line
    #     for _, row in df_bridle_line_system_surfplan.iterrows():
    #         # Extract p1 and p2
    #         p1 = [row["p1_x"], row["p1_y"], row["p1_z"]]
    #         p2 = [row["p2_x"], row["p2_y"], row["p2_z"]]

    #         # Deduplicate and register p1
    #         key1 = point_key(p1)
    #         if key1 in existing_node_map:
    #             idx1 = existing_node_map[key1]
    #         else:
    #             idx1 = len(struc_nodes_list)
    #             struc_nodes_list.append(p1)
    #             existing_node_map[key1] = idx1

    #         # Deduplicate and register p2
    #         key2 = point_key(p2)
    #         if key2 in existing_node_map:
    #             idx2 = existing_node_map[key2]
    #         else:
    #             idx2 = len(struc_nodes_list)
    #             struc_nodes_list.append(p2)
    #             existing_node_map[key2] = idx2

    #         # Store connectivity
    #         bridle_ci.append(idx1)
    #         bridle_cj.append(idx2)
    # else:
    #     # Extract brilde_nodes
    #     bridle_nodes = config_kite_dict["bridle_nodes"]
    #     headers = bridle_nodes["headers"]  # e.g. ["LE_x", "LE_y", ...]
    #     data = bridle_nodes["data"]  # list of lists
    #     bridle_nodes_data = config_kite_dict["bridle_nodes"]["data"]
    #     df_bridle_nodes = pd.DataFrame(data, columns=headers)

    #     pulley_point_indices = []
    #     for _, row in df_bridle_nodes.iterrows():
    #         struc_nodes_list.append(row[["x", "y", "z"]].tolist())
    #         if row["type"] == "pulley":
    #             pulley_point_indices.append(len(struc_nodes_list) - 1)

    # bridle_ci = config_kite_dict["bridle_lines"]["ci"]
    # bridle_cj = config_kite_dict["bridle_lines"]["cj"]

    # Given that we know

    # pulley_data: {
    #     "point_indices": [24, 28],
    #     "mass": 0.1,
    #     "number_of_pulleys_in_back_lines": 2,
    #     "line_indices": [np.int64(34), np.int64(32), np.int64(36), np.int64(33)],
    #     "line_pair_indices": [["6", 4], ["8", 5]],
    #     "ci": array([23, 27, 22, 22]),
    #     "cj": array([24, 28, 24, 28]),
    #     "other_line_pair": [
    #         ["6", array([22.0, 24.0, 2.86387114])],
    #         ["4", array([23.0, 24.0, 2.34058161])],
    #         ["8", array([22.0, 28.0, 2.86387114])],
    #         ["5", array([27.0, 28.0, 2.34058161])],
    #     ],
    #     "other_line_pair_corrected_dict": {
    #         "34": array([22.0, 24.0, 2.86387114]),
    #         "32": array([23.0, 24.0, 2.34058161]),
    #         "36": array([22.0, 28.0, 2.86387114]),
    #         "33": array([27.0, 28.0, 2.34058161]),
    #     },
    # }

    # # Extract brilde_nodes
    # bridle_nodes = config_kite_dict["bridle_nodes"]
    # headers = bridle_nodes["headers"]  # e.g. ["LE_x", "LE_y", ...]
    # data = bridle_nodes["data"]  # list of lists
    # bridle_nodes_data = config_kite_dict["bridle_nodes"]["data"]
    # df_bridle_nodes = pd.DataFrame(data, columns=headers)

    # pulley_point_indices = []
    # for _, row in df_bridle_nodes.iterrows():
    #     struc_nodes_list.append(row[["x", "y", "z"]].tolist())
    #     if row["type"] == "pulley":
    #         pulley_point_indices.append(len(struc_nodes_list) - 1)

    bridle_nodes_data = config_kite_dict["bridle_nodes"]["data"]
    for idx, (node_id, x, y, z, node_type) in enumerate(bridle_nodes_data):
        struc_nodes_list.append([x, y, z])
        if node_type == "pulley":
            pulley_point_indices.append(node_id)

    bridle_lines = config_kite_dict["bridle_lines"]
    headers = bridle_lines["headers"]  # e.g. ["LE_x", "LE_y", ...]
    data = bridle_lines["data"]  # list of lists
    bridle_lines_dict = {}
    for row in data:
        # Use first index of each data row as key, and the rest as value
        bridle_lines_dict[row[0]] = np.array(row[1:])

    node_masses = np.concatenate((node_masses_wing, np.zeros(len(bridle_nodes_data))))
    print(f"node_masses: {node_masses.shape}, {node_masses}")
    bridle_ci = []
    bridle_cj = []
    pulley_line_indices = []
    pulley_line_to_other_node_pair_dict = {}
    bridle_rest_lengths_initial = []
    for idx, bridle_data_i in enumerate(config_kite_dict["bridle_connections"]["data"]):
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
        print(
            f"line_name: {line_name}, line_density: {line_density}, line_diameter: {line_diameter}, line_length: {line_length}, line_mass: {line_mass}, node_masses[ci]: {node_masses[ci]}, node_masses[cj]: {node_masses[cj]}"
        )

        if line_name == "M":
            depower_tape_index = idx + n_wing_lines

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
        node_masses[idx] += config_kite_dict["pulley_mass"]

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
        depower_tape_index,
        node_masses,
    )


# def distribute_mass(
#     struc_nodes,
#     bridle_ci,
#     bridle_cj,
#     wing_ci,
#     wing_cj,
#     pulley_point_indices,
#     config_kite_dict,
# ):
#     """
#     Distribute mass to each node based on structural and bridle configuration.

#     Args:
#         struc_nodes (np.ndarray): Structural node positions (n_nodes,3).
#         bridle_ci, bridle_cj (list): Bridle connectivity indices.
#         wing_ci, wing_cj (list): Wing connectivity indices.
#         pulley_point_indices (list): Indices of pulley nodes.
#         config_kite_dict (dict): Kite configuration dictionary.

#     Returns:
#         np.ndarray: Mass array for each node (n_nodes,).
#     """

#     WING_MASS = config_kite_dict["wing_mass"]
#     BRIDLE_RHO = config_kite_dict["bridle"]["density"]
#     BRIDLE_DIAMETER = config_kite_dict["bridle"]["diameter"]
#     KCU_MASS = config_kite_dict["kcu"]["mass"]
#     KCU_INDEX = config_kite_dict["kcu"]["index"]
#     PULLEY_MASS = config_kite_dict["pulley"]["mass"]

#     # Define a spring force matrix of the right size
#     node_masses = np.zeros(
#         struc_nodes.shape[0]
#     )  # Initialising with zero matrix in same shape as points

#     ## Bridle lines
#     for idx, (idx_bridle_node_i, idx_bridle_node_j) in enumerate(
#         zip(bridle_ci, bridle_cj)
#     ):  # loop through each bridle line
#         # Calculate the len of the bridle line
#         len_bridle = np.linalg.norm(
#             struc_nodes[idx_bridle_node_i] - struc_nodes[idx_bridle_node_j]
#         )
#         # Calculate the mass of the bridle line
#         mass_bridle = BRIDLE_RHO * np.pi * (BRIDLE_DIAMETER / 2) ** 2 * len_bridle
#         # Add the mass of the bridle line to the nodes
#         node_masses[idx_bridle_node_i] += mass_bridle / 2
#         node_masses[idx_bridle_node_j] += mass_bridle / 2

#     # print(f'bridle node_masses: {np.sum(node_masses)}')
#     ## Pulleys
#     for idx in pulley_point_indices:
#         node_masses[idx] += PULLEY_MASS

#     # print(f'pulley& bridle node_masses: {np.sum(node_masses)}')
#     ## KCU
#     node_masses[KCU_INDEX] += KCU_MASS

#     # print(f'pulley& bridle & KCU node_masses: {np.sum(node_masses)}')

#     ## Wing
#     for idx, (idx_wing_node_i, idx_wing_node_i) in enumerate(
#         zip(set(wing_ci), set(wing_cj))
#     ):  # making them sets, to have it unique
#         node_masses[idx] += WING_MASS / len(set(wing_ci))

#     # print(f'pulley& bridle & KCU & wing node_masses: {np.sum(node_masses)}')

#     return node_masses


# # TODO: this function is way to bulky, it should be possible to make this simpler
# def extract_pulley_connectivity(
#     struc_nodes, bridle_ci, bridle_cj, wing_connectivity, config_kite_dict
# ):
#     """
#     Extract pulley connectivity and mapping information.

#     Args:
#         struc_nodes (np.ndarray): Structural node positions (n_nodes,3).
#         bridle_ci, bridle_cj (list): Bridle connectivity indices.
#         wing_connectivity (np.ndarray): Wing connectivity matrix.
#         config_kite_dict (dict): Kite configuration dictionary.

#     Returns:
#         dict: Pulley data dictionary with connectivity and mapping.
#     """

#     pulley_data = {
#         "point_indices": config_kite_dict["pulley"]["point_indices"],
#         "mass": config_kite_dict["pulley"]["mass"],
#         "number_of_pulleys_in_back_lines": config_kite_dict["pulley"][
#             "number_of_pulleys_in_back_lines"
#         ],
#     }

#     PULLEY_point_indices = pulley_data["point_indices"]
#     number_of_pulleys_in_back_lines = pulley_data["number_of_pulleys_in_back_lines"]

#     # print(f'PULLEY_point_indices {len(PULLEY_point_indices)}, number_of_pulleys_in_back_lines: {number_of_pulleys_in_back_lines}')
#     # Sorting pulleys based on z-coordinate
#     sorted_pulley_point_indices = sorted(
#         PULLEY_point_indices, key=lambda index: struc_nodes[index][2]
#     )

#     # TODO: this is where you input the number of pulleys in the back-lines, that have a different feature
#     # for pulleys where the line is below the pulley point
#     pulley_point_indices_line_below = sorted_pulley_point_indices[
#         :number_of_pulleys_in_back_lines
#     ]
#     if (
#         len(PULLEY_point_indices) > number_of_pulleys_in_back_lines
#     ):  # if there are pulleys in the front lines
#         pulley_point_indices_line_above = sorted_pulley_point_indices[
#             number_of_pulleys_in_back_lines:
#         ]
#     else:
#         pulley_point_indices_line_above = []  # this should be the case for the V3.25

#     all_possible_line_indices = []

#     for idx, (idx_bridle_node_i, idx_bridle_node_j) in enumerate(
#         zip(bridle_ci, bridle_cj)
#     ):  # loop through each bridle line
#         # if the current line its index i OR j  is a pulley point
#         # AND the other line is LOWER than the pulley_point, i.e. pulley line is BELOW the pulley
#         if (
#             idx_bridle_node_i in pulley_point_indices_line_below
#             and struc_nodes[idx_bridle_node_i][2] > struc_nodes[idx_bridle_node_j][2]
#         ) or (
#             idx_bridle_node_j in pulley_point_indices_line_below
#             and struc_nodes[idx_bridle_node_j][2] > struc_nodes[idx_bridle_node_i][2]
#         ):
#             all_possible_line_indices.append(idx)

#         # if the current line its index i OR j  is a pulley point
#         # AND the other line is HIGHER than the pulley_point, i.e. pulley line is ABOVE the pulley
#         if (
#             idx_bridle_node_i in pulley_point_indices_line_above
#             and struc_nodes[idx_bridle_node_i][2] < struc_nodes[idx_bridle_node_j][2]
#         ) or (
#             idx_bridle_node_j in pulley_point_indices_line_above
#             and struc_nodes[idx_bridle_node_j][2] < struc_nodes[idx_bridle_node_i][2]
#         ):
#             all_possible_line_indices.append(idx)

#     # loop through all_possible_line_indices, twice to try and find the matching line pair
#     pulley_line_pair_indices = {}
#     pulley_line_indices = []

#     for line_index_1 in all_possible_line_indices:
#         for line_index_2 in all_possible_line_indices:
#             # break if the same line is compared
#             if line_index_1 == line_index_2:
#                 break

#             # list the point_indices of the two lines
#             point_indices = [
#                 bridle_ci[line_index_1],
#                 bridle_cj[line_index_1],
#                 bridle_ci[line_index_2],
#                 bridle_cj[line_index_2],
#             ]
#             # sort the point_indices based on the z-coordinate of the points[index]
#             sorted_indices = sorted(
#                 point_indices, key=lambda index: struc_nodes[index][2]
#             )

#             # IF the pulley line is BELOW the pulley point, i.e. IF LAST two indices correspond to the same point
#             # AND the index is in the pulley_index list (i.e. it is actually a pulley)
#             # AND the line is not already used, i.e. not in the used_pulley_line_indices set
#             if (
#                 sorted_indices[2] == sorted_indices[3]
#                 and sorted_indices[3] in pulley_point_indices_line_below
#                 and line_index_1 not in pulley_line_indices
#                 and line_index_2 not in pulley_line_indices
#             ):
#                 # Append the line indices, to the pulley_line_indices
#                 pulley_line_indices.append(line_index_1)
#                 pulley_line_indices.append(line_index_2)
#                 # Make new key indices for the pulley line pair indices
#                 pulley_line_pair_indices[str(line_index_1)] = line_index_2

#             # IF the pulley line is ABOVE the pulley point, i.e. IF FIRST indices correspond to the same point
#             # AND the index is in the pulley_index list (i.e. it is actually a pulley)
#             # AND the line is not already used, i.e. not in the used_pulley_line_indices set
#             elif (
#                 sorted_indices[0] == sorted_indices[1]
#                 and sorted_indices[0] in pulley_point_indices_line_above
#                 and line_index_1 not in pulley_line_indices
#                 and line_index_2 not in pulley_line_indices
#             ):
#                 # Append the line indices, to the pulley_line_indices
#                 pulley_line_indices.append(line_index_1)
#                 pulley_line_indices.append(line_index_2)
#                 # Make new key indices for the pulley line pair indices
#                 pulley_line_pair_indices[str(line_index_1)] = line_index_2

#     # i know that line key and line value make up a pulley
#     # so i want to plot both line key and line value
#     # extract them separately and append
#     pulley_ci_key = [
#         bridle_ci[int(key_index)] for key_index in pulley_line_pair_indices.keys()
#     ]
#     pulley_cj_key = [
#         bridle_cj[int(key_index)] for key_index in pulley_line_pair_indices.keys()
#     ]
#     pulley_ci_value = [
#         bridle_ci[int(value_index)] for value_index in pulley_line_pair_indices.values()
#     ]
#     pulley_cj_value = [
#         bridle_cj[int(value_index)] for value_index in pulley_line_pair_indices.values()
#     ]
#     pulley_ci = np.concatenate((pulley_ci_key, pulley_ci_value))
#     pulley_cj = np.concatenate((pulley_cj_key, pulley_cj_value))

#     pulley_data["line_indices"] = np.ndarray.flatten(np.array(pulley_line_indices))

#     # transform pulley_line_pair indices from a dict to a nested list
#     # Why? Because attrs dataclasses can't handle numbers as attributes
#     # Does not seem used anywhere
#     pulley_line_pair_indices = [
#         [key, value] for key, value in pulley_line_pair_indices.items()
#     ]

#     pulley_data["line_pair_indices"] = pulley_line_pair_indices
#     pulley_data["ci"] = pulley_ci
#     pulley_data["cj"] = pulley_cj

#     additional_dict = {}
#     for i in range(0, len(pulley_data["line_indices"]), 2):
#         # line 1, key: line 1 and value: line 2 data
#         line_1_key = str(pulley_data["line_indices"][i])
#         line_1_idx_p3 = bridle_ci[pulley_data["line_indices"][i + 1]]
#         line_1_idx_p4 = bridle_cj[pulley_data["line_indices"][i + 1]]
#         line_1_rest_len_p3p4 = np.linalg.norm(
#             struc_nodes[line_1_idx_p3] - struc_nodes[line_1_idx_p4]
#         )
#         additional_dict[line_1_key] = np.array(
#             [line_1_idx_p3, line_1_idx_p4, line_1_rest_len_p3p4]
#         )

#         # line 2, key: line 2 and value: line 1 data
#         line_2_key = str(pulley_data["line_indices"][i + 1])
#         line_2_idx_p3 = bridle_ci[pulley_data["line_indices"][i]]
#         line_2_idx_p4 = bridle_cj[pulley_data["line_indices"][i]]
#         line_2_rest_len_p3p4 = np.linalg.norm(
#             struc_nodes[line_2_idx_p3] - struc_nodes[line_2_idx_p4]
#         )
#         additional_dict[line_2_key] = np.array(
#             [line_2_idx_p3, line_2_idx_p4, line_2_rest_len_p3p4]
#         )

#     # transform other_line_pair from a dict to a nested list
#     # Why? Because attrs dataclasses can't handle numbers as attributes
#     # USED in input_particleSystem for rotational springs
#     additional_dict = [[key, value] for key, value in additional_dict.items()]
#     pulley_data["other_line_pair"] = additional_dict

#     ## PARAMS
#     n_wing_elements = len(wing_connectivity)

#     other_line_pair_dict = {}
#     for entry in pulley_data["other_line_pair"]:
#         other_line_pair_dict[entry[0]] = entry[1]

#     # Correct the key
#     other_line_pair_corrected_dict = {}
#     for key_i in other_line_pair_dict.keys():
#         corrected_key = str(int(key_i) + n_wing_elements)
#         other_line_pair_corrected_dict[corrected_key] = other_line_pair_dict[key_i]

#     pulley_data["other_line_pair_corrected_dict"] = other_line_pair_corrected_dict

#     pulley_data["line_indices"] = [
#         idx + n_wing_elements for idx in pulley_data["line_indices"]
#     ]
#     return pulley_data


def main(config_kite_dict):
    """
    Main entry point for kite structural initialisation.

    Args:
        config_kite_dict (dict): Kite configuration dictionary.

    Returns:
        tuple: All structural and connectivity arrays, mass arrays, rest lengths, and pulley data.
    """

    (
        struc_nodes,
        wing_ci,
        wing_cj,
        wing_connectivity,
        wing_rest_lengths_initial,
        struc_le_idx_list,
        struc_te_idx_list,
        tubular_frame_line_idx_list,
        te_line_idx_list,
        n_struc_ribs,
        pulley_point_indices,
        node_masses_wing,
    ) = initialize_wing_structure(config_kite_dict)

    (
        struc_nodes,
        bridle_ci,
        bridle_cj,
        bridle_connectivity,
        bridle_rest_lengths_initial,
        pulley_point_indices,
        pulley_line_indices,
        pulley_line_to_other_node_pair_dict,
        depower_tape_index,
        node_masses,
    ) = initialize_bridle_line_system(
        struc_nodes,
        node_masses_wing,
        len(wing_connectivity),
        pulley_point_indices,
        config_kite_dict,
    )

    kite_connectivity = np.vstack((wing_connectivity, bridle_connectivity))
    rest_lengths = np.concatenate(
        (wing_rest_lengths_initial, bridle_rest_lengths_initial)
    )
    # m_array = distribute_mass(
    #     struc_nodes,
    #     bridle_ci,
    #     bridle_cj,
    #     wing_ci,
    #     wing_cj,
    #     pulley_point_indices,
    #     config_kite_dict,
    # )

    return (
        struc_nodes,
        wing_ci,
        wing_cj,
        bridle_ci,
        bridle_cj,
        struc_le_idx_list,
        struc_te_idx_list,
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
        depower_tape_index,
    )
