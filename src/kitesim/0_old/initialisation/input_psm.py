import numpy as np
from kitesim.initialisation import particles_with_rotational_resistance
from kitesim.post_processing import plotting


def define_connectivity_matrix(config):
    wing_connectivity = np.column_stack(
        (config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj)
    )
    bridle_connectivity = np.column_stack(
        (config.kite.connectivity.bridle_ci, config.kite.connectivity.bridle_cj)
    )
    connectivity_matrix = np.vstack((wing_connectivity, bridle_connectivity))

    return connectivity_matrix, wing_connectivity


def define_initial_conditions_kite(config):
    # defining parameters
    points_ini = np.array(config.kite.points_ini)
    if config.is_with_initial_point_velocity:
        print("Error: initial point velocity has never been defined")
    else:
        vel_ini = np.zeros(points_ini.shape)
    m_array = config.kite.mass_points
    fixed_nodes = np.array(config.kite.bridle.bridle_point_index)
    # fill with: position, initial velocity?, mass, fixed boolean
    conditions = []
    n = config.kite.n_points
    for i in range(n):
        if i in fixed_nodes:
            conditions.append([points_ini[i], vel_ini[i], m_array[i], True])
        else:
            conditions.append([points_ini[i], vel_ini[i], m_array[i], False])

    return conditions


def define_params(config, wing_connectivity, connectivity_matrix):
    bridle_rest_lengths = config.kite.bridle_rest_lengths_initial
    wing_rest_lengths = config.kite.wing_rest_lengths_initial
    rest_lengths = np.concatenate((wing_rest_lengths, bridle_rest_lengths))

    n_wing_elements = len(wing_connectivity)

    # Transform pulley.other_line_pair to a dict
    # data_struc is [["3",value],["5", value], ...]
    other_line_pair_dict = {}
    for entry in config.kite.pulley.other_line_pair:
        other_line_pair_dict[entry[0]] = entry[1]

    # Correct the key
    other_line_pair_corrected_dict = {}
    for key_i in other_line_pair_dict.keys():
        corrected_key = str(int(key_i) + n_wing_elements)
        other_line_pair_corrected_dict[corrected_key] = other_line_pair_dict[key_i]

    # initializing connectivity lists
    canopy_indices = []
    tube_indices = [i for i, conn in enumerate(wing_connectivity)]

    # initializing empty lists
    is_compression_list = []
    is_tension_list = []
    is_rotational_list = []
    stiffness_list = []
    is_pulley_list = []

    for i, conn in enumerate(connectivity_matrix):
        # if wing elements
        if float(i) < len(wing_connectivity):
            is_tension_list.append(True)
            is_pulley_list.append(False)

            if i in config.kite.connectivity.te_line_indices:
                stiffness_list.append(config.kite.stiffness.trailing_edge)
                is_compression_list.append(False)
                is_rotational_list.append(False)
            elif i in config.kite.connectivity.tube_line_indices:
                stiffness_list.append(config.kite.stiffness.tube)
                is_compression_list.append(True)
                is_rotational_list.append(True)
            else:  # must be a canopy element
                stiffness_list.append(config.kite.stiffness.canopy)
                is_compression_list.append(False)
                is_rotational_list.append(False)

        # if bridle-lines
        else:
            is_compression_list.append(False)
            is_tension_list.append(True)
            is_rotational_list.append(True)
            stiffness_list.append(config.kite.stiffness.bridle)

            # TODO: might be better to use one index style/structure?
            # the (i - n_wing_elements) is to correct the indices
            if (i - n_wing_elements) in config.kite.pulley.line_indices:
                # print(f'pulley, i: {i}')
                # print(f'connection[i]: {connectivity_matrix[i]}')
                is_pulley_list.append(True)
            else:
                is_pulley_list.append(False)

    params = {
        "c": config.solver.damping_constant,
        "dt": config.solver.dt,
        "t_steps": config.solver.n_time_steps,
        "abs_tol": config.solver.abs_tol,
        "rel_tol": config.solver.rel_tol,
        "max_iter": config.solver.max_iter,
        "pulley_other_line_pair": other_line_pair_corrected_dict,
        "k": np.array(stiffness_list),
        "is_compression": np.array(is_compression_list),
        "is_tension": np.array(is_tension_list),
        "is_pulley": np.array(is_pulley_list),
        "is_rotational": np.array(is_rotational_list),
        "n": int(len(config.kite.points_ini)),
        # "m_segment": 0.1,
        "aerostructural_tol": config.aero_structural.tol,  # N, < f_residual
        "l0": rest_lengths,
        "is_with_visc_damping": config.solver.is_with_visc_damping,
    }

    return params


def create_input_PSM(config):
    """Create the input PSM class.

    Args:
        config: The configuration class.

    Returns:
        The input PSM class."""
    # Should be the same for each kite
    (
        connectivity_matrix,
        wing_connectivity,
    ) = define_connectivity_matrix(config)
    params_dict = define_params(config, wing_connectivity, connectivity_matrix)
    initial_conditions = define_initial_conditions_kite(config)
    points_between_dict = (
        particles_with_rotational_resistance.extract_points_between_dict(config)
    )
    if config.is_with_initial_plot:
        plotting.plot_initial_geometry(config, points_between_dict)

    # TODO: this should be set in the kite config file
    is_with_rotational_resistance = False
    if config.kite_name == "V9_60C":
        is_with_rotational_resistance = True

    if is_with_rotational_resistance:
        (
            leadingedge_rotational_resistance_dict,
            strut_rotational_resistance_dict,
        ) = particles_with_rotational_resistance.extract_rotational_resistances_dicts(
            points_between_dict, config
        )
        # first do the struts
        params_dict = particles_with_rotational_resistance.initialize_bending_spring(
            config.kite.stiffness.k_bend_strut,
            initial_conditions,
            params_dict,
            connectivity_matrix,
            strut_rotational_resistance_dict,
        )
        # secondly do the leading-edge
        params_dict = particles_with_rotational_resistance.initialize_bending_spring(
            config.kite.stiffness.k_bend_leading_edge,
            initial_conditions,
            params_dict,
            connectivity_matrix,
            leadingedge_rotational_resistance_dict,
        )

    input_PSM = {
        "connectivity_matrix": connectivity_matrix,
        "initial_conditions": initial_conditions,
        "params_dict": params_dict,
    }
    return input_PSM
