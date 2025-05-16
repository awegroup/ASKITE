import numpy as np


# TODO: further extend this to also extract pulley_indices and so on.
def convert_points(points):
    """
    Convert a list of [idx, x, y, z, type, subtype] entries into:
      - coords:   (N,3) float array of [x,y,z]
      - ids:      (N,)   int array of the original point indices
      - le_idx:   (n_LE,) int array of point-indices whose subtype == "le"
      - te_idx:   (n_TE,) int array of point-indices whose subtype == "te"

    Parameters
    ----------
    points : list of lists
        Each element is [index, x, y, z, type, subtype].

    Returns
    -------
    coords : np.ndarray, shape (N,3)
    ids    : np.ndarray, shape (N,)
    le_idx : np.ndarray, shape (<=N,)
    te_idx : np.ndarray, shape (<=N,)
    """
    coords = []
    ids = []
    le_idx = []
    te_idx = []

    for idx, x, y, z, major, subtype in points:
        idx = int(idx)
        ids.append(idx)
        coords.append([x, y, z])

        # collect wing LE/TE indices
        if major.lower() == "wing":
            if subtype.lower() == "le":
                le_idx.append(idx)
            elif subtype.lower() == "te":
                te_idx.append(idx)

    return (
        np.array(coords, dtype=float),
        np.array(ids, dtype=int),
        np.array(le_idx, dtype=int),
        np.array(te_idx, dtype=int),
    )


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

    points, ids, le_idx, te_idx = convert_points(config_kite_dict["points"])

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
        le_idx,
        te_idx,
        wing_connectivity,
        bridle_connectivity,
        kite_connectivity,
        wing_rest_lengths_initial,
        rest_lengths,
        m_array,
    )
