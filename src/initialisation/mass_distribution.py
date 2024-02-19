## Mass-distribution
# import libraries
import numpy as np


def calculate_mass_distribution(
    points_ini,
    bridle_ci,
    bridle_cj,
    wing_ci,
    wing_cj,
    WING_MASS,
    BRIDLE_DATA,
    KCU_DATA,
    PULLEY_DATA,
):
    BRIDLE_RHO = BRIDLE_DATA["RHO"]
    BRIDLE_DIAMETER = BRIDLE_DATA["DIAMETER"]
    KCU_MASS = KCU_DATA["MASS"]
    KCU_INDEX = KCU_DATA["INDEX"]
    PULLEY_INDICES = PULLEY_DATA["POINT_INDICES"]
    PULLEY_MASS = PULLEY_DATA["MASS"]

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
