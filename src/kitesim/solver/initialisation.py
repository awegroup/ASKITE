import numpy as np


def convert_points(points):
    """
    Convert a list of [idx, x, y, z, type] entries into an (N,3) NumPy array of [x, y, z].

    Parameters
    ----------
    points : list of lists
        Each element is [index, x, y, z, some_type].

    Returns
    -------
    numpy.ndarray
        An array of shape (len(points), 3) with just the x, y, z values.
    """
    # Using tuple unpacking to skip first and last items
    xyz = [[x, y, z] for _, x, y, z, _ in points]
    return np.array(xyz, dtype=float)


def calculate_edge_lengths(ci, cj, pos):
    """returns the edge lengths between the nodes with index ci and cj
    for the given positions pos
    input : ci,cj,pos
    output: springL"""
    springL = np.zeros(ci.shape)
    for idx, (ci, cj) in enumerate(zip(ci, cj)):
        springL[idx] = np.linalg.norm(pos[cj, :] - pos[ci, :])
    return springL


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


def initialising_solver(config_kite_dict):

    points = convert_points(config_kite_dict["points"])

    ## extracting the connectivity from config_kite
    wing_connectivity = np.column_stack(
        (
            config_kite_dict["wing_connectivity"]["ci"],
            config_kite_dict["wing_connectivity"]["cj"],
        )
    )
    bridle_connectivity = np.column_stack(
        (
            config_kite_dict["bridle_connectivity"]["ci"],
            config_kite_dict["bridle_connectivity"]["cj"],
        )
    )
    kite_connectivity = np.vstack((wing_connectivity, bridle_connectivity))

    ## computing rest lengths
    bridle_rest_lengths_initial = np.array(
        calculate_edge_lengths(
            np.array(config_kite_dict["bridle_connectivity"]["ci"]),
            np.array(config_kite_dict["bridle_connectivity"]["cj"]),
            points,
        )
    )
    wing_rest_lengths_initial = np.array(
        calculate_edge_lengths(
            np.array(config_kite_dict["wing_connectivity"]["ci"]),
            np.array(config_kite_dict["wing_connectivity"]["ci"]),
            points,
        )
    )
    rest_lengths = np.concatenate(
        (wing_rest_lengths_initial, bridle_rest_lengths_initial)
    )

    ## initializing mass
    m_array = new_calculate_mass_distribution(
        points,
        config_kite_dict["bridle_connectivity"]["ci"],
        config_kite_dict["bridle_connectivity"]["cj"],
        config_kite_dict["wing_connectivity"]["ci"],
        config_kite_dict["wing_connectivity"]["cj"],
        config_kite_dict["wing_mass"],
        config_kite_dict["bridle"]["density"],
        config_kite_dict["bridle"]["diameter"],
        config_kite_dict["kcu"]["mass"],
        config_kite_dict["kcu"]["index"],
        config_kite_dict["pulley"]["point_indices"],
        config_kite_dict["pulley"]["mass"],
    )

    return (
        points,
        wing_connectivity,
        bridle_connectivity,
        kite_connectivity,
        wing_rest_lengths_initial,
        rest_lengths,
        m_array,
    )
