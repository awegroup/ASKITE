import numpy as np


def initialize_wing_structure(
    struc_geometry,
    struc_nodes,
    m_arr,
    conn_arr,
    l0_arr,
    k_arr,
    c_arr,
    linktype_arr,
):
    """
    Create the structural nodes and connectivity for the kite structure.

    Args:
        geometry_dict (dict): Kite configuration dictionary.

    Returns:
        tuple: (struc_nodes, wing_ci, wing_cj, bridle_ci, bridle_cj, struc_node_le_indices,
                struc_node_te_indices, pulley_point_indices, tubular_frame_line_idx_list,
                te_line_idx_list, n_struc_ribs)
    """
    ### node level ###
    # node_indices = []
    struc_node_le_indices = []
    struc_node_te_indices = []

    # Then append rest of the defined wing_nodes
    for node_idx, x, y, z in struc_geometry["wing_nodes"]["data"]:
        # node_indices.append(node_idx)
        struc_nodes.append(np.array([x, y, z]))
        m_arr.append(0)

        # if uneven --> this is a leading-edge node
        if node_idx % 2 != 0:
            struc_node_le_indices.append(node_idx)
        else:
            struc_node_te_indices.append(node_idx)

    ### element level ###
    # Create an element dict of dicts: { name → {rest_length:..., diameter:..., ...} }
    wing_elements_dict = {
        row[0]: dict(
            zip(struc_geometry["wing_spring_damper_elements"]["headers"][1:], row[1:])
        )
        for row in struc_geometry["wing_spring_damper_elements"]["data"]
    }
    tubular_frame_line_idx_list = []
    te_line_idx_list = []
    for conn_idx, (conn_name, ci, cj) in enumerate(
        struc_geometry["wing_connections"]["data"]
    ):

        m_element = wing_elements_dict[conn_name]["m"]
        m_arr[ci] += m_element / 2
        m_arr[cj] += m_element / 2

        conn_arr.append([ci, cj])
        l0_arr.append(wing_elements_dict[conn_name]["l0"])
        k_arr.append(wing_elements_dict[conn_name]["k"])
        c_arr.append(wing_elements_dict[conn_name]["c"])
        linktype_arr.append(wing_elements_dict[conn_name]["linktype"])

        if "le" in conn_name.lower() or "strut" in conn_name.lower():
            tubular_frame_line_idx_list.append(conn_idx)
        elif "te" in conn_name.lower():
            te_line_idx_list.append(conn_idx)

    return (
        # node level
        struc_nodes,
        m_arr,
        struc_node_le_indices,
        struc_node_te_indices,
        # element level
        conn_arr,
        l0_arr,
        k_arr,
        c_arr,
        linktype_arr,
        tubular_frame_line_idx_list,
        te_line_idx_list,
    )


def initialize_bridle_line_system(
    struc_geometry,
    struc_nodes,
    m_arr,
    conn_arr,
    l0_arr,
    k_arr,
    c_arr,
    linktype_arr,
):
    """
    Initialize the bridle line system for the kite.

    Returns:
        tuple: (bridle_ci, bridle_cj, pulley_point_indices)
    """

    ### node level ###

    # First append the bridle_point_node, as this node (KCU) should have index 0
    # Then append rest of the defined bridle_nodes

    for node_idx, x, y, z in struc_geometry["bridle_nodes"]["data"]:
        struc_nodes.append(np.array([x, y, z]))
        m_arr.append(0.0)

    ### element level ###

    # Create an element dict of dicts: { name → {l0:..., d:..., ...} }
    bridle_lines_dict = {
        row[0]: dict(zip(struc_geometry["bridle_lines"]["headers"][1:], row[1:]))
        for row in struc_geometry["bridle_lines"]["data"]
    }

    n_wing_conn = len(conn_arr)
    pulley_node_indices = []
    pulley_line_indices = []
    pulley_line_to_other_node_pair_dict = {}
    steering_tape_indices = []
    for conn_idx, conn_data in enumerate(struc_geometry["bridle_connections"]["data"]):

        conn_name = conn_data[0]
        ci = int(conn_data[1])
        cj = int(conn_data[2])
        l0 = bridle_lines_dict[conn_name]["l0"]
        material = bridle_lines_dict[conn_name]["material"]
        cross_sectional_area = np.pi * (bridle_lines_dict[conn_name]["d"] / 2) ** 2
        m_line = struc_geometry[material]["density"] * cross_sectional_area * l0

        conn_arr.append([ci, cj])
        m_arr[ci] += m_line / 2
        m_arr[cj] += m_line / 2

        if len(conn_data[1:]) == 2:
            linktype_arr.append("noncompressive")

            k = (struc_geometry[material]["youngs_modulus"] * cross_sectional_area) / l0
            c = (
                struc_geometry[material]["damping_per_stiffness"] * k
            )  # Rayleigh damping
            l0_arr.append(l0)
            k_arr.append(k)
            c_arr.append(c)

        elif len(conn_data[1:]) == 3:
            linktype_arr.append("pulley")

            # if the connection ck exists, i.e. if there is a third connection, this line is a pulley
            pulley_line_index = conn_idx + n_wing_conn
            pulley_line_indices.append(pulley_line_index)

            ck = int(conn_data[3])
            pulley_index = cj
            pulley_node_indices.append(pulley_index)
            m_arr[pulley_index] += struc_geometry["pulley_mass"]

            # Compute the straight-line distances between the nodes involved in the pulley:
            # - point_to_point_line_len_current: distance between ci and cj (the current segment)
            # - point_to_point_line_len_other: distance between cj and ck (the other segment)
            point_to_point_line_len_current = np.linalg.norm(
                np.array(struc_nodes[ci]) - np.array(struc_nodes[cj])
            )
            point_to_point_line_len_other = np.linalg.norm(
                np.array(struc_nodes[cj]) - np.array(struc_nodes[ck])
            )

            # The total straight-line length is the sum of both segments
            total_point_to_point_line_len = (
                point_to_point_line_len_current + point_to_point_line_len_other
            )

            # Divide the rest length proportionally between the two segments,
            # based on their straight-line distances ratios to the total length
            line_len_current = (
                point_to_point_line_len_current / total_point_to_point_line_len
            ) * l0
            line_len_other = (
                point_to_point_line_len_other / total_point_to_point_line_len
            ) * l0

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
            k_eff = (
                struc_geometry[material]["youngs_modulus"] * cross_sectional_area
            ) / (line_len_current + line_len_other)
            c_eff = (
                struc_geometry[material]["damping_per_stiffness"] * k_eff
            )  # Rayleigh damping
            l0_arr.append(line_len_current)
            k_arr.append(k_eff)
            c_arr.append(c_eff)

        else:
            raise ValueError(
                "bridle_connections should have 2 or 3 connections (ci,cj,ck), not more or less"
            )

        if conn_name == "Power Tape":
            power_tape_index = conn_idx + n_wing_conn

        if conn_name == "Steering Tape":
            steering_tape_indices.append(conn_idx + n_wing_conn)

    return (
        # node level
        struc_nodes,
        m_arr,
        power_tape_index,
        steering_tape_indices,
        pulley_node_indices,
        # element level
        conn_arr,
        l0_arr,
        k_arr,
        c_arr,
        linktype_arr,
        pulley_line_to_other_node_pair_dict,
    )


def main(struc_geometry):

    ### First append the bridle_point_node, as this node (KCU) should have index 0
    struc_nodes = []
    m_arr = []
    struc_nodes.append(np.array(struc_geometry["bridle_point_node"]))
    m_arr.append(struc_geometry["kcu_mass"])

    # initialize element level lists
    conn_arr = []
    l0_arr = []
    k_arr = []
    c_arr = []
    linktype_arr = []

    ### Analyze Wing Structure
    (
        # node level
        struc_nodes,
        m_arr,
        struc_node_le_indices,
        struc_node_te_indices,
        # element level
        conn_arr,
        l0_arr,
        k_arr,
        c_arr,
        linktype_arr,
        tubular_frame_line_idx_list,
        te_line_idx_list,
    ) = initialize_wing_structure(
        struc_geometry,
        struc_nodes,
        m_arr,
        conn_arr,
        l0_arr,
        k_arr,
        c_arr,
        linktype_arr,
    )

    ### Analyze Bridle Structure
    (
        # node level
        struc_nodes,
        m_arr,
        power_tape_index,
        steering_tape_indices,
        pulley_node_indices,
        # element level
        conn_arr,
        l0_arr,
        k_arr,
        c_arr,
        linktype_arr,
        pulley_line_to_other_node_pair_dict,
    ) = initialize_bridle_line_system(
        struc_geometry, struc_nodes, m_arr, conn_arr, l0_arr, k_arr, c_arr, linktype_arr
    )

    # explicit numpy arrays
    struc_nodes = np.array(struc_nodes)
    conn_arr = np.array(conn_arr)
    l0_arr = np.array(l0_arr)
    k_arr = np.array(k_arr)
    c_arr = np.array(c_arr)
    m_arr = np.array(m_arr)
    linktype_arr = np.array(linktype_arr)
    struc_node_le_indices = np.array(struc_node_le_indices)
    struc_node_te_indices = np.array(struc_node_te_indices)

    return (
        # node level
        struc_nodes,
        m_arr,
        struc_node_le_indices,
        struc_node_te_indices,
        power_tape_index,
        steering_tape_indices,
        pulley_node_indices,
        # element level
        conn_arr,
        l0_arr,
        k_arr,
        c_arr,
        linktype_arr,
        pulley_line_to_other_node_pair_dict,
    )
