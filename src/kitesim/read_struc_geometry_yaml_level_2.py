import numpy as np
import logging


def initialize_particles(struc_geometry, struc_nodes, m_arr):
    """
    Initialize particles for the kite structure.

    This function adds wing and bridle particles to the structural nodes list
    and initializes their masses to zero. It also identifies leading edge and
    trailing edge node indices based on even/odd node indexing, and generates
    additional canopy section nodes between struts.

    Args:
        struc_geometry (dict): Dictionary containing structural geometry data
        struc_nodes (list): List to append particle positions to
        m_arr (list): List to append particle masses to

    Returns:
        tuple: (struc_nodes, m_arr, struc_node_le_indices, struc_node_te_indices,
                strut_node_le_indices, strut_node_te_indices, canopy_sections)
            Updated lists with particle positions, masses, edge node indices,
            strut node indices, and canopy section connectivity
    """

    struc_node_le_indices = []
    struc_node_te_indices = []

    for node_idx, x, y, z in struc_geometry["wing_particles"]["data"]:
        # node_indices.append(node_idx)
        struc_nodes.append(np.array([x, y, z]))
        m_arr.append(0)
        # if uneven --> this is a leading-edge node
        if node_idx % 2 != 0:
            struc_node_le_indices.append(node_idx)
        else:
            struc_node_te_indices.append(node_idx)

    for node_idx, x, y, z in struc_geometry["bridle_particles"]["data"]:
        # node_indices.append(node_idx)
        struc_nodes.append(np.array([x, y, z]))
        m_arr.append(0)

    # add extra nodes for canopy
    strut_node_le_indices = []
    strut_node_te_indices = []
    strut_indices = []
    nodes_per_strut = int(0)
    for (
        name,
        ci,
        cj,
        strut_diam_le,
        strut_diam_te,
        le_diameter,
        node_indices,
    ) in struc_geometry["strut_tubes"]["data"]:
        strut_node_le_indices.append(ci)
        strut_node_te_indices.append(cj)
        nodes_per_strut = max(len(node_indices), nodes_per_strut)
        strut_indices.append(node_indices)

    canopy_section_le_indices = [
        idx for idx in struc_node_le_indices if idx not in strut_node_le_indices
    ]
    canopy_section_te_indices = [
        idx for idx in struc_node_te_indices if idx not in strut_node_te_indices
    ]
    canopy_sections = []

    # TODO: Add extra nodes along chord here, make input through configuration file?
    nodes_per_strut += 1

    # add extra nodes along struts such that the amount per strut is the same
    for i, indices in enumerate(strut_indices):
        nodes = len(indices)
        missing_nodes = nodes_per_strut - nodes
        for i in range(missing_nodes):
            coords_front = struc_nodes[indices[-3]]
            coords_back = struc_nodes[indices[-2]]
            ratio = (i + 1) / (missing_nodes + 1)
            x = coords_front[0] + ratio * (coords_back[0] - coords_front[0])
            y = coords_front[1] + ratio * (coords_back[1] - coords_front[1])
            z = coords_front[2] + ratio * (coords_back[2] - coords_front[2])
            struc_nodes.append(np.array([x, y, z]))
            m_arr.append(0)
            node_idx += 1
            indices.insert(-2, node_idx)

    simplified_bridle_points = []
    for i, indices in enumerate(strut_indices):
        # Project all intermediate nodes onto the line between first and last node
        start_pos = struc_nodes[indices[0]]
        end_pos = struc_nodes[indices[-1]]
        line_vector = end_pos - start_pos
        for j in range(1, len(indices) - 1):
            idx = indices[j]
            current_pos = struc_nodes[idx]
            # Project current node onto the line between start and end
            projection_scalar = np.dot(current_pos - start_pos, line_vector) / np.dot(
                line_vector, line_vector
            )
            projected_pos = start_pos + projection_scalar * line_vector
            struc_nodes[idx] = projected_pos
            simplified_bridle_points.append([idx, current_pos, projected_pos])
    simplified_bridle_points = np.array(simplified_bridle_points, dtype=object)

    for n1, n2 in zip(canopy_section_le_indices, canopy_section_te_indices):
        # Find closest indices to n1 in strut_node_le_indices
        neg_distances = [(idx - n1) for idx in strut_node_le_indices if idx < n1]
        strut_right_le = (
            strut_node_le_indices[
                strut_node_le_indices.index(
                    min(neg_distances, key=abs, default=n1) + n1
                )
            ]
            if neg_distances
            else n1
        )
        # Find closest index with positive offset
        pos_distances = [(idx - n1) for idx in strut_node_le_indices if idx > n1]
        strut_left_le = (
            strut_node_le_indices[
                strut_node_le_indices.index(
                    min(pos_distances, key=abs, default=n1) + n1
                )
            ]
            if pos_distances
            else n1
        )
        # Find closest indices to n2 in strut_node_te_indices
        neg_distances = [(idx - n1) for idx in strut_node_te_indices if idx < n1]
        strut_right_te = (
            strut_node_te_indices[
                strut_node_te_indices.index(
                    min(neg_distances, key=abs, default=n1) + n1
                )
            ]
            if neg_distances
            else n1
        )
        # Find closest index with positive offset
        pos_distances = [(idx - n2) for idx in strut_node_te_indices if idx > n2]
        strut_left_te = (
            strut_node_te_indices[
                strut_node_te_indices.index(
                    min(pos_distances, key=abs, default=n2) + n2
                )
            ]
            if pos_distances
            else n2
        )

        leading_edge_tube_indices = np.arange(strut_right_le, strut_left_le + 2, 2)
        trailing_edge_tube_indics = np.arange(strut_right_te, strut_left_te + 2, 2)
        leading_edge_tube_length = sum(
            np.linalg.norm(
                struc_nodes[leading_edge_tube_indices[i]]
                - struc_nodes[leading_edge_tube_indices[i + 1]]
            )
            for i in range(len(leading_edge_tube_indices) - 1)
        )
        trailing_edge_tube_length = sum(
            np.linalg.norm(
                struc_nodes[trailing_edge_tube_indics[i]]
                - struc_nodes[trailing_edge_tube_indics[i + 1]]
            )
            for i in range(len(trailing_edge_tube_indics) - 1)
        )
        ratio_le = (
            np.linalg.norm(struc_nodes[n1] - struc_nodes[strut_right_le])
            / leading_edge_tube_length
        )
        ratio_te = (
            np.linalg.norm(struc_nodes[n1] - struc_nodes[strut_right_le])
            / trailing_edge_tube_length
        )
        ratio_canopy = (ratio_le + ratio_te) / 2

        length_right = np.linalg.norm(
            struc_nodes[strut_right_le] - struc_nodes[strut_right_te]
        )
        length_left = np.linalg.norm(
            struc_nodes[strut_left_le] - struc_nodes[strut_left_te]
        )

        strut_indices_right = strut_indices[strut_node_le_indices.index(strut_right_le)]
        strut_indices_left = strut_indices[strut_node_le_indices.index(strut_left_le)]

        canopy_section_length = np.linalg.norm(struc_nodes[n1] - struc_nodes[n2])
        direction = struc_nodes[n2] - struc_nodes[n1]
        direction_normalized = direction / np.linalg.norm(direction)
        canopy_section_indices = [n1]
        for n in range(1, nodes_per_strut - 1):
            l_ratio_right = (
                np.linalg.norm(
                    struc_nodes[strut_indices_right[n]] - struc_nodes[strut_right_le]
                )
                / length_right
            )
            l_ratio_left = (
                np.linalg.norm(
                    struc_nodes[strut_indices_left[n]] - struc_nodes[strut_left_le]
                )
                / length_left
            )
            l_ratio = l_ratio_right * (1 - ratio_canopy) + l_ratio_left * ratio_canopy
            coordinates = (
                struc_nodes[n1] + canopy_section_length * l_ratio * direction_normalized
            )
            struc_nodes.append(coordinates)
            m_arr.append(0)
            node_idx += 1
            canopy_section_indices.append(node_idx)
        canopy_section_indices.append(n2)
        canopy_sections.append(canopy_section_indices)
    strut_sections = strut_indices
    return (
        struc_nodes,
        m_arr,
        struc_node_le_indices,
        struc_node_te_indices,
        strut_node_le_indices,
        strut_node_te_indices,
        canopy_sections,
        strut_sections,
        simplified_bridle_points,
    )


def initialize_wing_structure(
    struc_geometry,
    struc_nodes,
    m_arr,
    kite_connectivity_arr,
    l0_arr,
    k_arr,
    c_arr,
    linktype_arr,
    canopy_sections,
    strut_sections,
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

    # Struts
    for (
        name,
        ci,
        cj,
        strut_diam_le,
        strut_diam_te,
        le_diameter,
        node_indices,
    ) in struc_geometry["strut_tubes"]["data"]:
        c1s = node_indices[0:-1]
        c2s = node_indices[1:]
        length = np.linalg.norm(
            struc_nodes[node_indices[0]] - struc_nodes[node_indices[-1]]
        )
        for c1, c2 in zip(c1s, c2s):
            rest_length = np.linalg.norm(struc_nodes[c1] - struc_nodes[c2])
            # determine diameter of element, scale linearly from le to te
            l1 = np.linalg.norm(struc_nodes[ci] - struc_nodes[c1])
            l2 = np.linalg.norm(struc_nodes[ci] - struc_nodes[c2])
            diameter_n1 = strut_diam_le - (strut_diam_le - strut_diam_te) / length * l1
            diameter_n2 = strut_diam_le - (strut_diam_le - strut_diam_te) / length * l2
            diameter = (diameter_n1 + diameter_n2) / 2
            mass = 2 * np.pi * (diameter / 2) * rest_length * 170 / 1000
            # m_arr[c1] += mass/2
            # m_arr[c2] += mass/2
            kite_connectivity_arr.append([c1, c2])
            l0_arr.append(rest_length)
            k_arr.append(diameter)  # use k array to store diameter
            c_arr.append(struc_geometry["pressure"])  # use c array to store pressure
            linktype_arr.append("inflatable_beam")
    # for struct_sect in strut_sections:
    #     print(struct_sect)

    # leading edge tube
    for name, ci, cj, diameter in struc_geometry["leading_edge_tubes"]["data"]:
        rest_length = np.linalg.norm(struc_nodes[ci] - struc_nodes[cj])
        mass = 2 * np.pi * (diameter / 2) * rest_length * 170 / 1000
        # m_arr[ci] += mass/2
        # m_arr[cj] += mass/2
        kite_connectivity_arr.append([ci, cj])
        l0_arr.append(rest_length)
        k_arr.append(diameter)  # use k array to store diameter
        c_arr.append(struc_geometry["pressure"])  # use c array to store pressure
        linktype_arr.append("inflatable_beam")

    # combine and order canopy_and strut sections by first indices
    all_sections = canopy_sections + strut_sections
    all_sections.sort(key=lambda section: section[0])

    # Connect canopy sections along chord
    for canopy_section in canopy_sections:
        c1s = canopy_section[0:-1]
        c2s = canopy_section[1:]
        for c1, c2 in zip(c1s, c2s):
            rest_length = np.linalg.norm(struc_nodes[c1] - struc_nodes[c2])
            kite_connectivity_arr.append([c1, c2])
            l0_arr.append(rest_length)
            k_arr.append(5000)
            c_arr.append(0)
            linktype_arr.append("noncompressive")

    # canopy connections (crosses and rectangles)
    all_sections_1 = all_sections[0:-1]
    all_sections_2 = all_sections[1:]
    for section1, section2 in zip(all_sections_1, all_sections_2):
        # Connect corresponding nodes between adjacent struts with squares and diagonals
        # Skip the first connection (section1[0] to section2[0])
        for i in range(1, len(section1)):
            if i < len(section2):
                # Square connections: section1[i] to section2[i]
                rest_length = np.linalg.norm(
                    struc_nodes[section1[i]] - struc_nodes[section2[i]]
                )
                kite_connectivity_arr.append([section1[i], section2[i]])
                l0_arr.append(rest_length)
                k_arr.append(5000)
                c_arr.append(0)
                linktype_arr.append("noncompressive")

                # Diagonal connections: section1[i] to section2[i-1] (if i > 0)
                if i > 0:
                    rest_length = np.linalg.norm(
                        struc_nodes[section1[i]] - struc_nodes[section2[i - 1]]
                    )
                    kite_connectivity_arr.append([section1[i], section2[i - 1]])
                    l0_arr.append(rest_length)
                    k_arr.append(5000)
                    c_arr.append(0)
                    linktype_arr.append("noncompressive")

                # Diagonal connections: section1[i-1] to section2[i] (if i > 0)
                if i > 0:
                    rest_length = np.linalg.norm(
                        struc_nodes[section1[i - 1]] - struc_nodes[section2[i]]
                    )
                    kite_connectivity_arr.append([section1[i - 1], section2[i]])
                    l0_arr.append(rest_length)
                    k_arr.append(5000)
                    c_arr.append(0)
                    linktype_arr.append("noncompressive")

    # assign masses canopy
    def triangle_area(p1, p2, p3):
        # Vectors for two sides of the triangle
        v1 = p2 - p1
        v2 = p3 - p1
        # Cross product magnitude gives 2 * triangle area
        return 0.5 * np.linalg.norm(np.cross(v1, v2))

    def quad_area(A, B, C, D):
        A = np.array(A)
        B = np.array(B)
        C = np.array(C)
        D = np.array(D)

        # Split quad into triangles ABC and ACD
        area1 = triangle_area(A, B, C)
        area2 = triangle_area(A, C, D)
        return area1 + area2

    for i in range(len(canopy_sections) - 1):
        section_a = canopy_sections[i]
        section_b = canopy_sections[i + 1]
        # Create quads by connecting adjacent nodes in consecutive sections
        for j in range(len(section_a) - 1):
            quad = [section_a[j], section_a[j + 1], section_b[j + 1], section_b[j]]
            # Get coordinates and calculate area
            corners = [struc_nodes[node] for node in quad]
            area = quad_area(corners[0], corners[1], corners[2], corners[3])
            # Distribute quad mass to its 4 nodes (each gets 1/4 of quad mass)
            quad_mass = (
                area * struc_geometry["canopy_density"] / 1000
            )  # Convert g/m^2 to kg/m^2
            for node in quad:
                m_arr[node] += quad_mass / 4

    mass_canopy = np.sum(m_arr)

    le_indices = np.array(all_sections)[:, 0]

    # Calculate total length of inflatable tubes
    total_inflatable_length = 0

    # Add up lengths in leading edge
    for i in range(len(le_indices) - 1):
        total_inflatable_length += np.linalg.norm(
            struc_nodes[le_indices[i]] - struc_nodes[le_indices[i + 1]]
        )

    # Add up lengths in strut sections
    for section in strut_sections:
        for j in range(len(section) - 1):
            total_inflatable_length += np.linalg.norm(
                struc_nodes[section[j]] - struc_nodes[section[j + 1]]
            )

    # Calculate target mass for inflatable tubes
    target_mass_inflatable = struc_geometry["mass_without_bridles"] - mass_canopy

    # Distribute mass along leading edge based on segment lengths
    for i in range(len(le_indices) - 1):
        segment_length = np.linalg.norm(
            struc_nodes[le_indices[i]] - struc_nodes[le_indices[i + 1]]
        )
        segment_mass = target_mass_inflatable * (
            segment_length / total_inflatable_length
        )
        m_arr[le_indices[i]] += segment_mass / 2
        m_arr[le_indices[i + 1]] += segment_mass / 2

    # Distribute mass along strut sections based on segment lengths
    for section in strut_sections:
        for j in range(len(section) - 1):
            segment_length = np.linalg.norm(
                struc_nodes[section[j]] - struc_nodes[section[j + 1]]
            )
            segment_mass = target_mass_inflatable * (
                segment_length / total_inflatable_length
            )
            m_arr[section[j]] += segment_mass / 2
            m_arr[section[j + 1]] += segment_mass / 2

    return (
        # node level
        struc_nodes,
        m_arr,
        # element level
        kite_connectivity_arr,
        l0_arr,  # l
        k_arr,  # d
        c_arr,  # p
        linktype_arr,
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
    simplified_bridle_points,
):
    """
    Initialize the bridle line system for the kite.

    Returns:
        tuple: (bridle_ci, bridle_cj, pulley_point_indices)
    """

    ### node level ###

    # First append the bridle_point_node, as this node (KCU) should have index 0
    # Then append rest of the defined bridle_particles

    # for node_idx, x, y, z in struc_geometry["bridle_particles"]["data"]:
    #     struc_nodes.append(np.array([x, y, z]))
    #     m_arr.append(0.0)

    ### element level ###

    # Create an element dict of dicts: { name → {l0:..., d:..., ...} }
    bridle_lines_dict = {
        row[0]: dict(zip(struc_geometry["bridle_lines"]["headers"][1:], row[1:]))
        for row in struc_geometry["bridle_lines"]["data"]
    }

    # initialize a connectivity counter, that starts with the number of wing_connections
    conn_idx_counter = len(kite_connectivity_arr)
    bridle_connectivity_arr = []
    bridle_diameter_arr = []
    pulley_node_indices = []
    pulley_line_indices = []
    pulley_line_to_other_node_pair_dict = {}
    steering_tape_indices = []
    power_tape_index = 0
    for _, conn_data in enumerate(struc_geometry["bridle_connections"]["data"]):

        conn_name = conn_data[0]
        ci = int(conn_data[1])
        cj = int(conn_data[2])

        # computing the mass of the bridle line, and adding it 0.5 to each particle using m_arr
        l0 = bridle_lines_dict[conn_name]["rest_length"]
        material = bridle_lines_dict[conn_name]["material"]
        cross_sectional_area = (
            np.pi * (bridle_lines_dict[conn_name]["diameter"] / 2) ** 2
        )
        m_line = struc_geometry[material]["density"] * cross_sectional_area * l0

        # If there is third connections, this line is a pulley!
        # In here we will treat both ci-cj and cj-ck
        if len(conn_data[1:]) == 3:
            logging.debug(
                f"-- linktype should be pulley, linktype: {bridle_lines_dict[conn_name]["linktype"]}"
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
            # updating l0 to include simplification in bridle attachment point

            if ci in simplified_bridle_points[:, 0]:
                i = np.where(simplified_bridle_points[:, 0] == ci)[0][0]
                original_pos = simplified_bridle_points[i, 1]
                projected_pos = simplified_bridle_points[i, 2]
                original_lenth = np.linalg.norm(original_pos - struc_nodes[ck])
                new_length = np.linalg.norm(projected_pos - struc_nodes[ck])
                delta = original_lenth - new_length
                l0 -= delta

            if ck in simplified_bridle_points[:, 0]:
                i = np.where(simplified_bridle_points[:, 0] == ck)[0][0]
                original_pos = simplified_bridle_points[i, 1]
                projected_pos = simplified_bridle_points[i, 2]
                original_lenth = np.linalg.norm(original_pos - struc_nodes[ck])
                new_length = np.linalg.norm(projected_pos - struc_nodes[ck])
                delta = original_lenth - new_length
                l0 -= delta

            #######################################
            # Dealing with ci-cj
            # add this new connection to the connectivity array, and also increase counter
            kite_connectivity_arr.append([ci, cj])
            bridle_connectivity_arr.append([ci, cj])
            bridle_diameter_arr.append(bridle_lines_dict[conn_name]["diameter"])
            l0_arr.append(l0)
            k_arr.append(5000)
            c_arr.append(c)
            linktype_arr.append(bridle_lines_dict[conn_name]["linktype"])

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
            bridle_diameter_arr.append(bridle_lines_dict[conn_name]["diameter"])
            l0_arr.append(l0)
            k_arr.append(5000)
            c_arr.append(c)
            linktype_arr.append(bridle_lines_dict[conn_name]["linktype"])
            # add mass
            m_arr[ci] += m_line / 4
            m_arr[cj] += m_line / 2
            m_arr[ck] += m_line / 4

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
                f"-- linktype should be noncompressive, linktype: {bridle_lines_dict[conn_name]["linktype"]}"
            )
            # add this new connection to the connectivity array, and also increase counter
            k = (struc_geometry[material]["youngs_modulus"] * cross_sectional_area) / l0
            c = (
                struc_geometry[material]["damping_per_stiffness"] * k
            )  # Rayleigh damping
            #######################################
            # updating l0 to include simplification in bridle attachment point
            if ci in simplified_bridle_points[:, 0]:
                i = np.where(simplified_bridle_points[:, 0] == ci)[0][0]
                original_pos = simplified_bridle_points[i, 1]
                projected_pos = simplified_bridle_points[i, 2]
                original_lenth = np.linalg.norm(original_pos - struc_nodes[cj])
                new_length = np.linalg.norm(projected_pos - struc_nodes[cj])
                delta = original_lenth - new_length
                l0 -= delta
            elif cj in simplified_bridle_points[:, 0]:
                i = np.where(simplified_bridle_points[:, 0] == cj)[0][0]
                original_pos = simplified_bridle_points[i, 1]
                projected_pos = simplified_bridle_points[i, 2]
                original_lenth = np.linalg.norm(original_pos - struc_nodes[ci])
                new_length = np.linalg.norm(projected_pos - struc_nodes[ci])
                delta = original_lenth - new_length
                l0 -= delta

            m_arr[ci] += m_line / 2
            m_arr[cj] += m_line / 2
            kite_connectivity_arr.append([ci, cj])
            bridle_connectivity_arr.append([ci, cj])
            bridle_diameter_arr.append(bridle_lines_dict[conn_name]["diameter"])
            l0_arr.append(l0)
            k_arr.append(5000)
            c_arr.append(c)
            linktype_arr.append(bridle_lines_dict[conn_name]["linktype"])

        else:
            raise ValueError(
                "bridle_connections should have 2 or 3 connections (ci,cj,ck), not more or less"
            )
        try:
            if conn_name == "Power Tape":
                power_tape_index = conn_idx_counter

            if conn_name == "Steering Tape":
                steering_tape_indices.append(conn_idx_counter)
        except:
            power_tape_index = 0
            steering_tape_indices.append(0)

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
    m_arr.append(0)

    # initialize element level lists
    kite_connectivity_arr = []
    l0_arr = []
    k_arr = []
    c_arr = []
    linktype_arr = []

    # initialize particles
    (
        struc_nodes,
        m_arr,
        struc_node_le_indices,
        struc_node_te_indices,
        strut_node_le_indices,
        strut_node_te_indices,
        canopy_sections,
        strut_sections,
        simplified_bridle_points,
    ) = initialize_particles(struc_geometry, struc_nodes, m_arr)

    ### Analyze Wing Structure
    (
        # node level
        struc_nodes,
        m_arr,
        # element level
        kite_connectivity_arr,
        l0_arr,
        k_arr,
        c_arr,
        linktype_arr,
    ) = initialize_wing_structure(
        struc_geometry,
        struc_nodes,
        m_arr,
        kite_connectivity_arr,
        l0_arr,
        k_arr,
        c_arr,
        linktype_arr,
        canopy_sections,
        strut_sections,
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
        simplified_bridle_points,
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
        canopy_sections,
        strut_sections,
        simplified_bridle_points,
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
