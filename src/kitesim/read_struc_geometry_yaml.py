import numpy as np
import logging


def initialize_wing_structure(
    struc_geometry,
    struc_nodes,
    m_arr,
    kite_connectivity_arr,
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

    # Then append rest of the defined wing_particles
    for node_idx, x, y, z in struc_geometry["wing_particles"]["data"]:
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
        row[0]: dict(zip(struc_geometry["wing_elements"]["headers"][1:], row[1:]))
        for row in struc_geometry["wing_elements"]["data"]
    }
    tubular_frame_line_idx_list = []
    te_line_idx_list = []
    for conn_idx, (conn_name, ci, cj) in enumerate(
        struc_geometry["wing_connections"]["data"]
    ):

        m_element = wing_elements_dict[conn_name]["m"]
        m_arr[ci] += m_element / 2
        m_arr[cj] += m_element / 2

        kite_connectivity_arr.append([ci, cj])
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
        kite_connectivity_arr,
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
    kite_connectivity_arr,
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
    # Then append rest of the defined bridle_particles

    for node_idx, x, y, z in struc_geometry["bridle_particles"]["data"]:
        struc_nodes.append(np.array([x, y, z]))
        m_arr.append(0.0)

    ### element level ###

    # Create an element dict of dicts: { name → {l0:..., d:..., ...} }
    bridle_elements_dict = {
        row[0]: dict(zip(struc_geometry["bridle_elements"]["headers"][1:], row[1:]))
        for row in struc_geometry["bridle_elements"]["data"]
    }

    # initialize a connectivity counter, that starts with the number of wing_connections
    conn_idx_counter = len(kite_connectivity_arr)
    bridle_connectivity_arr = []
    bridle_diameter_arr = []
    pulley_node_indices = []
    pulley_line_indices = []
    pulley_line_to_other_node_pair_dict = {}
    steering_tape_indices = []
    for _, conn_data in enumerate(struc_geometry["bridle_connections"]["data"]):

        conn_name = conn_data[0]
        ci = int(conn_data[1])
        cj = int(conn_data[2])

        # computing the mass of the bridle line, and adding it 0.5 to each particle using m_arr
        l0 = bridle_elements_dict[conn_name]["l0"]
        material = bridle_elements_dict[conn_name]["material"]
        cross_sectional_area = np.pi * (bridle_elements_dict[conn_name]["d"] / 2) ** 2
        m_line = struc_geometry[material]["density"] * cross_sectional_area * l0
        m_arr[ci] += m_line / 2
        m_arr[cj] += m_line / 2

        # If there is third connections, this line is a pulley!
        # In here we will treat both ci-cj and cj-ck
        if len(conn_data[1:]) == 3:
            logging.debug(
                f"-- linktype should be pulley, linktype: {bridle_elements_dict[conn_name]["linktype"]}"
            )
            # adding pulley_node_indices
            pulley_node_indices.append(cj)

            # making the third node an integer
            ck = int(conn_data[3])

            # add the pulley mass, to the pulley_index cj
            m_arr[cj] += struc_geometry["pulley_mass"]

            #######################################
            # Computing k,c
            # Compute the straight-line distances between the nodes involved in the pulley:
            # - len_ci_cj: distance between ci and cj (the current segment)
            # - len_cj_ck: distance between cj and ck (the other segment)
            len_ci_cj = np.linalg.norm(
                np.array(struc_nodes[ci]) - np.array(struc_nodes[cj])
            )
            len_cj_ck = np.linalg.norm(
                np.array(struc_nodes[cj]) - np.array(struc_nodes[ck])
            )

            # The total straight-line length is the sum of both segments
            len_ci_cj_ck = len_ci_cj + len_cj_ck

            # Divide the rest length proportionally between the two segments,
            # based on their straight-line distances ratios to the total length
            l0_len_ci_cj = (len_ci_cj / len_ci_cj_ck) * l0
            l0_len_cj_ck = (len_cj_ck / len_ci_cj_ck) * l0

            k = (struc_geometry[material]["youngs_modulus"] * cross_sectional_area) / (
                l0
            )
            c = struc_geometry[material]["damping_per_stiffness"] * k

            #######################################
            # Dealing with ci-cj
            # add this new connection to the connectivity array, and also increase counter
            kite_connectivity_arr.append([ci, cj])
            bridle_connectivity_arr.append([ci, cj])
            bridle_diameter_arr.append(bridle_elements_dict[conn_name]["d"])
            l0_arr.append(l0)
            k_arr.append(k)
            c_arr.append(c)
            linktype_arr.append(bridle_elements_dict[conn_name]["linktype"])

            # Create a special mapping for the Structural Particle System Solver
            # key: pulley_line_index
            # value: [cj, ck, line_len_other]
            # This is used to connect the pulley line to the other node pair
            pulley_line_to_other_node_pair_dict[str(conn_idx_counter)] = np.array(
                [
                    cj,
                    ck,
                    l0_len_cj_ck,
                    l0_len_ci_cj,
                    ci,
                ]
            )

            # Mark the indices of connectivity
            pulley_line_index_ci_cj = conn_idx_counter
            pulley_line_indices.append(pulley_line_index_ci_cj)
            conn_idx_counter += 1

            #######################################
            # Dealing with cj-ck
            # add this new connection to the connectivity array, and also increase counter
            kite_connectivity_arr.append([cj, ck])
            bridle_connectivity_arr.append([cj, ck])
            bridle_diameter_arr.append(bridle_elements_dict[conn_name]["d"])
            l0_arr.append(l0)
            k_arr.append(k)
            c_arr.append(c)
            linktype_arr.append(bridle_elements_dict[conn_name]["linktype"])

            # Create a special mapping for the Structural Particle System Solver
            # key: pulley_line_index
            # value: [cj, ck, line_len_other]
            # This is used to connect the pulley line to the other node pair
            pulley_line_to_other_node_pair_dict[str(conn_idx_counter)] = np.array(
                [
                    cj,
                    ci,
                    l0_len_ci_cj,
                    l0_len_cj_ck,
                    ci,
                ]
            )

            # Mark the indices of connectivity
            pulley_line_index_cj_ck = conn_idx_counter
            pulley_line_indices.append(pulley_line_index_cj_ck)

        # if there is no third connections this line represents a knot-to-knot line, a regular spring damper
        elif len(conn_data[1:]) == 2:
            logging.debug(
                f"-- linktype should be noncompressive, linktype: {bridle_elements_dict[conn_name]["linktype"]}"
            )
            # add this new connection to the connectivity array, and also increase counter
            k = (struc_geometry[material]["youngs_modulus"] * cross_sectional_area) / l0
            c = (
                struc_geometry[material]["damping_per_stiffness"] * k
            )  # Rayleigh damping
            kite_connectivity_arr.append([ci, cj])
            bridle_connectivity_arr.append([ci, cj])
            bridle_diameter_arr.append(bridle_elements_dict[conn_name]["d"])
            l0_arr.append(l0)
            k_arr.append(k)
            c_arr.append(c)
            linktype_arr.append(bridle_elements_dict[conn_name]["linktype"])

        else:
            raise ValueError(
                "bridle_connections should have 2 or 3 connections (ci,cj,ck), not more or less"
            )

        if conn_name == "Power Tape":
            power_tape_index = conn_idx_counter

        if conn_name == "Steering Tape":
            steering_tape_indices.append(conn_idx_counter)

        ## increasing the counter
        conn_idx_counter += 1

    return (
        # node level
        struc_nodes,
        m_arr,
        power_tape_index,
        steering_tape_indices,
        pulley_node_indices,
        # element level
        kite_connectivity_arr,
        bridle_connectivity_arr,
        bridle_diameter_arr,
        l0_arr,
        k_arr,
        c_arr,
        linktype_arr,
        pulley_line_indices,
        pulley_line_to_other_node_pair_dict,
    )


def main(struc_geometry):

    ### First append the bridle_point_node, as this node (KCU) should have index 0
    struc_nodes = []
    m_arr = []
    struc_nodes.append(np.array(struc_geometry["bridle_point_node"]))
    m_arr.append(struc_geometry["kcu_mass"])

    # initialize element level lists
    kite_connectivity_arr = []
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
        kite_connectivity_arr,
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
        kite_connectivity_arr,
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
        kite_connectivity_arr,
        bridle_connectivity_arr,
        bridle_diameter_arr,
        l0_arr,
        k_arr,
        c_arr,
        linktype_arr,
        pulley_line_indices,
        pulley_line_to_other_node_pair_dict,
    ) = initialize_bridle_line_system(
        struc_geometry,
        struc_nodes,
        m_arr,
        kite_connectivity_arr,
        l0_arr,
        k_arr,
        c_arr,
        linktype_arr,
    )

    # explicit numpy arrays
    struc_nodes = np.array(struc_nodes)
    kite_connectivity_arr = np.array(kite_connectivity_arr)
    bridle_connectivity_arr = np.array(bridle_connectivity_arr)
    bridle_diameter_arr = np.array(bridle_diameter_arr)
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
        kite_connectivity_arr,
        bridle_connectivity_arr,
        bridle_diameter_arr,
        l0_arr,
        k_arr,
        c_arr,
        linktype_arr,
        pulley_line_indices,
        pulley_line_to_other_node_pair_dict,
    )
