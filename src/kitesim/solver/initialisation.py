import numpy as np
import pandas as pd


def create_struc_nodes(config_kite_dict):
    """Creates the nodes for the structural system"""

    # Extract airfoil
    af = config_kite_dict["airfoils"]
    headers = af["headers"]  # e.g. ["LE_x", "LE_y", ...]
    data = af["data"]  # list of lists
    airfoil_data = config_kite_dict["airfoils"]["data"]
    df_airfoil = pd.DataFrame(data, columns=headers)

    wing_ci = []
    wing_cj = []

    # 1) Pick out exactly the panel‐indices you’ll represent in the structure:
    n = len(df_airfoil)
    selected = [
        i for i in range(n) if i == 0 or df_airfoil.loc[i, "is_strut"] or i == n - 1
    ]

    # 2) Build your structural‐node list and remember for each panel which
    struc_nodes_list = [[0, 0, 0]]
    #    structural index is its LE and which is its TE
    le_map = {}
    te_map = {}
    struc_le_idx_list = []
    struc_te_idx_list = []
    for i in selected:
        row = df_airfoil.loc[i]
        # append LE
        le_idx = len(struc_nodes_list)
        struc_nodes_list.append(row[["LE_x", "LE_y", "LE_z"]].tolist())
        # append TE
        te_idx = len(struc_nodes_list)
        struc_nodes_list.append(row[["TE_x", "TE_y", "TE_z"]].tolist())
        le_map[i] = le_idx
        te_map[i] = te_idx
        # append to the list of indices
        struc_le_idx_list.append(le_idx)
        struc_te_idx_list.append(te_idx)

    tubular_frame_line_idx_list = []
    te_line_idx_list = []
    # 3) Now build wing_ci/wing_cj with explicit if‐branches
    for pos, i in enumerate(selected):
        curr_le = le_map[i]
        curr_te = te_map[i]

        if pos != len(selected) - 1:
            # --- For first to up to second to last panel only connect to the next ---
            wing_ci.append(curr_le)
            wing_cj.append(curr_te)
            tubular_frame_line_idx_list.append(len(wing_ci) - 1)
            # connect to next panel’s LE (LE[i] -- LE[i+1])
            nxt = selected[pos + 1]
            wing_ci.append(curr_le)
            wing_cj.append(le_map[nxt])
            tubular_frame_line_idx_list.append(len(wing_ci) - 1)
            # connect to next panel’s TE (TE[i] -- TE[i+1])
            wing_ci.append(curr_te)
            wing_cj.append(te_map[nxt])
            te_line_idx_list.append(len(wing_ci) - 1)
        elif pos == len(selected) - 1:
            # --- last panel only needs to connect to itself ---
            # connect LE->TE
            wing_ci.append(curr_le)
            wing_cj.append(curr_te)
            tubular_frame_line_idx_list.append(len(wing_ci) - 1)

    for i, struc_node in enumerate(struc_nodes_list):
        if i == 0:
            continue
        struc_nodes_list[i][0] = struc_node[0] + 0.7
        struc_nodes_list[i][2] = struc_node[2] + 7.3

    if config_kite_dict["is_with_bridle_line_system_surfplan"]:
        df_bridle_line_system_surfplan = pd.DataFrame(
            config_kite_dict["df_bridle_line_system_surfplan"]["data"],
            columns=config_kite_dict["df_bridle_line_system_surfplan"]["headers"],
        )

        # Initialize connectivity lists
        bridle_ci = []
        bridle_cj = []

        # Helper function to compare 3D points (rounded to tolerance)
        def point_key(p, tol=1e-8):
            return tuple(
                np.round(p, decimals=8)
            )  # prevent floating point uniqueness issues

        # Build a dictionary of existing nodes for fast lookup
        existing_node_map = {
            point_key(node): idx for idx, node in enumerate(struc_nodes_list)
        }

        # Iterate over each bridle line
        for _, row in df_bridle_line_system_surfplan.iterrows():
            # Extract p1 and p2
            p1 = [row["p1_x"], row["p1_y"], row["p1_z"]]
            p2 = [row["p2_x"], row["p2_y"], row["p2_z"]]

            # Deduplicate and register p1
            key1 = point_key(p1)
            if key1 in existing_node_map:
                idx1 = existing_node_map[key1]
            else:
                idx1 = len(struc_nodes_list)
                struc_nodes_list.append(p1)
                existing_node_map[key1] = idx1

            # Deduplicate and register p2
            key2 = point_key(p2)
            if key2 in existing_node_map:
                idx2 = existing_node_map[key2]
            else:
                idx2 = len(struc_nodes_list)
                struc_nodes_list.append(p2)
                existing_node_map[key2] = idx2

            # Store connectivity
            bridle_ci.append(idx1)
            bridle_cj.append(idx2)
    else:
        # Extract brilde_nodes
        bridle_nodes = config_kite_dict["bridle_nodes"]
        headers = bridle_nodes["headers"]  # e.g. ["LE_x", "LE_y", ...]
        data = bridle_nodes["data"]  # list of lists
        bridle_nodes_data = config_kite_dict["bridle_nodes"]["data"]
        df_bridle_nodes = pd.DataFrame(data, columns=headers)

        pulley_point_indices = []
        for _, row in df_bridle_nodes.iterrows():
            struc_nodes_list.append(row[["x", "y", "z"]].tolist())
            if row["type"] == "pulley":
                pulley_point_indices.append(len(struc_nodes_list) - 1)

        bridle_ci = config_kite_dict["bridle_lines"]["ci"]
        bridle_cj = config_kite_dict["bridle_lines"]["cj"]

        # Given that we know

    n_struc_ribs = len(selected)

    return (
        np.array(struc_nodes_list),
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
    )


def distribute_mass(
    points_ini,
    bridle_ci,
    bridle_cj,
    wing_ci,
    wing_cj,
    pulley_point_indices,
    config_kite_dict,
):

    WING_MASS = config_kite_dict["wing_mass"]
    BRIDLE_RHO = config_kite_dict["bridle"]["density"]
    BRIDLE_DIAMETER = config_kite_dict["bridle"]["diameter"]
    KCU_MASS = config_kite_dict["kcu"]["mass"]
    KCU_INDEX = config_kite_dict["kcu"]["index"]
    PULLEY_MASS = config_kite_dict["pulley"]["mass"]

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
    for idx in pulley_point_indices:
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


def calculate_edge_lengths(ci, cj, pos):
    """returns the edge lengths between the nodes with index ci and cj
    for the given positions pos
    input : ci,cj,pos
    output: springL"""
    springL = np.zeros(ci.shape)
    for idx, (ci, cj) in enumerate(zip(ci, cj)):
        springL[idx] = np.linalg.norm(pos[cj, :] - pos[ci, :])
    return springL


def main(config_kite_dict):
    (
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
    ) = create_struc_nodes(config_kite_dict)
    wing_connectivity = np.column_stack(
        (
            wing_ci,
            wing_cj,
        )
    )
    bridle_connectivity = np.column_stack(
        (
            bridle_ci,
            bridle_cj,
        )
    )
    kite_connectivity = np.vstack((wing_connectivity, bridle_connectivity))
    m_array = distribute_mass(
        struc_nodes,
        bridle_ci,
        bridle_cj,
        wing_ci,
        wing_cj,
        pulley_point_indices,
        config_kite_dict,
    )
    bridle_rest_lengths_initial = np.array(
        calculate_edge_lengths(
            np.array(bridle_ci),
            np.array(bridle_cj),
            struc_nodes,
        )
    )
    wing_rest_lengths_initial = np.array(
        calculate_edge_lengths(
            np.array(wing_ci),
            np.array(wing_cj),
            struc_nodes,
        )
    )
    rest_lengths = np.concatenate(
        (wing_rest_lengths_initial, bridle_rest_lengths_initial)
    )
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
        m_array,
        bridle_rest_lengths_initial,
        wing_rest_lengths_initial,
        rest_lengths,
    )
