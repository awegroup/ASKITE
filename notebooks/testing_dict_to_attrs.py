import attr
import logging
import numpy as np

config_data = {
    "kite_name": "V3_25",
    "sim_name": "straight",
    "u_p": 1.0,
    "depower_tape_extension_percentage": 8.0,
    "delta_ls_ini": 0.01,
    "depower_tape_extension_step": 0.01,
    "depower_tape_final_extension": 0.01,
    "steering_tape_extension_step": 0.05,
    "steering_tape_final_extension": -0.2,
    "is_with_printing": False,
    "is_print_mid_results": True,
    "is_with_initial_plot": False,
    "is_with_initial_point_velocity": False,
    "is_with_plotly_plot": False,
    "is_with_aero_geometry": True,
    "is_with_plotting": True,
    "plot_format": "png",
    "plot_elev": [9.14, 9.14],
    "plot_azim": [90.0, 230.0],
    "is_with_animation": False,
    "animation_fps": 10,
    "animation_dpi": 100,
    "animation_bitrate": 2000,
    "animation_elev": 9.14,
    "animation_azim": 230.0,
    "is_run_only_1_time_step": True,
    "bridle_initial_compression_factor": 1.0,
    "geometric_scaling_factor": 1.0,
    "is_with_velocity_initialization": False,
    "vel_app_initial": [10, 0, 0],
    "n_vel_initialisation_steps": 10,
    "is_billowing_on": True,
    "is_with_aero_bridle": True,
    "is_with_aero_tether": True,
    "coupling_method": "NN",
    "grav_constant": [0, 0, -9.81],
    "rho": 1.225,
    "mu": 1.789e-05,
}

config_data_2 = {
    "AeroConfig": {
        "ring_geometry": "5fil",
        "model": "VSM",
        "n_iter": 1500,
        "error_limit": 1e-05,
        "relax_factor": 0.03,
        "n_splits": 2,
        "n_chordwise_aero_nodes": 10,
        "cd_cylinder": 1.1,
        "cd_shear_cylinder": 0.02,
        "cd_multiplier": 1.0,
    },
    "SolverConfig": {
        "method": "root",
        "tol": 20.0,
        "max_fev": 2000,
        "integration_method": "df-sane",
        "newton_rtol": 10.0,
        "newton_max_iter": 2000,
        "newton_disp": False,
        "damping_constant": 0.9,
        "is_with_visc_damping": True,
        "dt": 0.01,
        "n_time_steps": 10000,
        "abs_tol": 1e-50,
        "rel_tol": 1e-05,
        "max_iter": 50000,
    },
    "AeroStructuralConfig": {
        "max_iter": 2500,
        "it_check": 5,
        "tol": 25.0,
        "crosswind_max_iter": 50,
        "crosswind_tol": 0.01,
        "crosswind_relax_factor": 0.5,
    },
    "TetherConfig": {"diameter": 0.04, "length": 100.0, "density": 1.0},
    "tol_fx_ratio_to_fz": 0.0,
    "tol_vk_optimization": 0.0,
    "vk_x_initial_guess_factor_of_vw": 0.0,
    "is_with_varying_va": False,
    "r_0_initial": 0.0,
    "is_with_vk_optimization": False,
    "is_circular_case": False,
    "is_print_intermediate_results": True,
    "is_with_gravity": False,
    "vel_wind": [12, 0, 2],
    "vel_kite": [0, 0, 0],
    "acc_kite": [0, 0, 0],
    "aero": {
        "ring_geometry": "5fil",
        "model": "VSM",
        "n_iter": 1500,
        "error_limit": 1e-05,
        "relax_factor": 0.03,
        "n_splits": 2,
        "n_chordwise_aero_nodes": 10,
        "cd_cylinder": 1.1,
        "cd_shear_cylinder": 0.02,
        "cd_multiplier": 1.0,
    },
    "solver": {
        "method": "root",
        "tol": 20.0,
        "max_fev": 2000,
        "integration_method": "df-sane",
        "newton_rtol": 10.0,
        "newton_max_iter": 2000,
        "newton_disp": False,
        "damping_constant": 0.9,
        "is_with_visc_damping": True,
        "dt": 0.01,
        "n_time_steps": 10000,
        "abs_tol": 1e-50,
        "rel_tol": 1e-05,
        "max_iter": 50000,
    },
    "aero_structural": {
        "max_iter": 2500,
        "it_check": 5,
        "tol": 25.0,
        "crosswind_max_iter": 50,
        "crosswind_tol": 0.01,
        "crosswind_relax_factor": 0.5,
    },
}

config_data_3 = {
    "tether": {"diameter": 0.04, "length": 100.0, "density": 1.0},
}
config_data_kite_1 = {
    "kite": {
        "points_ini": np.array(
            [
                [0.0, 0.0, 0.0],
                [1.54503268, 4.1300398, 7.26136495],
                [-0.01769469, 3.9841196, 8.49703339],
                [-0.23864989, 3.14708537, 9.81674623],
                [-0.38529401, 1.96770111, 10.60735655],
                [-0.45841928, 0.66695416, 10.94897378],
                [-0.45841928, -0.66695416, 10.94897378],
                [-0.38529401, -1.96770111, 10.60735655],
                [-0.23864989, -3.14708537, 9.81674623],
                [-0.01769469, -3.9841196, 8.49703339],
                [1.54503268, -4.1300398, 7.26136495],
                [1.71039665, -3.97159687, 8.49204017],
                [2.01066217, -3.12943185, 9.78719788],
                [2.118729, -1.95450565, 10.55854696],
                [2.16733995, -0.66340879, 10.88984164],
                [2.16733995, 0.66340879, 10.88984164],
                [2.118729, 1.95450565, 10.55854696],
                [2.01066217, 3.12943185, 9.78719788],
                [1.71039665, 3.97159687, 8.49204017],
                [0.86307674, 4.1565, 7.42381981],
                [0.86307674, -4.1565, 7.42381981],
                [0.17884119, 0.0, 0.93301438],
                [0.41766592, 0.0, 2.00472644],
                [0.47651992, 0.51570164, 2.41872723],
                [1.01186906, 0.85948465, 4.67117852],
                [1.44706777, 2.71993277, 7.35496513],
                [1.53424698, 1.34379336, 7.83340153],
                [0.47651992, -0.51570164, 2.41872723],
                [1.01186906, -0.85948465, 4.67117852],
                [1.44706777, -2.71993277, 7.35496513],
                [1.53424698, -1.34379336, 7.83340153],
                [-0.06579256, 1.3268467, 5.5315951],
                [-0.04266184, 2.05530303, 7.25510078],
                [-0.09211119, 1.26740889, 7.82771638],
                [-0.06579256, -1.3268467, 5.5315951],
                [-0.04266184, -2.05530303, 7.25510078],
                [-0.09211119, -1.26740889, 7.82771638],
            ]
        ),
        "n_points": 37,
        "surfplan_filename": "V3_25",
        "area_projected": 19.488047451623245,
        "area_surface": 25.0,
        "ref_chord": 2.625759225757671,
        "span": 8.313,
        "height": 3.687608834028282,
        "wing_mass": 11.0,
        "is_with_elongation_limit": True,
        "elongation_limit": 1.15,
        "is_with_compression_limit": True,
        "compression_limit": 0.85,
        "limit_stiffness_factor": 10.0,
        "billowing_angles": np.array(
            [1.0, 10.0, 15.0, 18.0, 20.0, 18.0, 15.0, 10.0, 1.0]
        ),
    }
}
config_data_kite_2 = {
    "kite": {
        "n_segments": 9,
        "wing_rest_lengths_initial": np.array(
            [
                1.39902134,
                1.72814392,
                1.25180306,
                0.70153804,
                1.37594874,
                1.99757003,
                1.57831856,
                2.2495754,
                1.5737963,
                1.72814392,
                2.55132262,
                2.49668735,
                1.42741599,
                2.50453344,
                1.40964868,
                2.2495754,
                2.74403898,
                2.78619536,
                1.34684531,
                2.62642736,
                1.33381044,
                2.50453344,
                2.88043685,
                2.90721778,
                1.33390831,
                2.62642736,
                1.32681758,
                2.62642736,
                2.94414227,
                2.94414227,
                1.34684531,
                2.50453344,
                1.33381044,
                2.62642736,
                2.90721778,
                2.88043685,
                1.42741599,
                2.2495754,
                1.40964868,
                2.50453344,
                2.78619536,
                2.74403898,
                1.57831856,
                1.72814392,
                1.5737963,
                2.2495754,
                2.49668735,
                2.55132262,
                1.39902134,
                0.70153804,
                1.25180306,
                1.72814392,
                1.99757003,
                1.37594874,
            ]
        ),
        "bridle_rest_lengths_initial": np.array(
            [
                0.95,
                1.098,
                1.60059473,
                1.60059473,
                2.86387114,
                2.86387114,
                2.34058161,
                6.13647292,
                2.34058161,
                6.13647292,
                3.29444619,
                3.24146386,
                3.29444619,
                3.24146386,
                1.71141601,
                2.53003642,
                1.71141601,
                2.53003642,
                2.85324488,
                3.19461359,
                2.85324488,
                3.19461359,
                5.68888347,
                5.68888347,
                1.87127109,
                2.29704124,
                1.87127109,
                2.29704124,
                2.29419993,
                2.79149197,
                2.29419993,
                2.79149197,
                2.88145187,
                3.19952736,
                2.88145187,
                3.19952736,
                4.29761369,
                4.29761369,
                3.52849123,
                3.52849123,
            ]
        ),
        "mass_points": np.array(
            [
                0.66134536,
                0.60542511,
                0.57072139,
                0.57521298,
                0.5760255,
                0.57889838,
                0.57889838,
                0.5760255,
                0.57521298,
                0.57072139,
                0.60542511,
                0.56545764,
                0.57285149,
                0.57577073,
                0.578854,
                0.578854,
                0.57577073,
                0.57285149,
                0.56545764,
                0.620686,
                0.070686,
                8.44741109,
                0.06165064,
                0.09102212,
                0.24485627,
                0.06806483,
                0.08390189,
                0.09102212,
                0.24485627,
                0.06806483,
                0.08390189,
                0.12090058,
                0.06283584,
                0.07567094,
                0.12090058,
                0.06283584,
                0.07567094,
            ]
        ),
        "bridle": {
            "diameter": 0.01,
            "density": 230.0,
            "bridle_point_index": 0,
            "depower_tape_index": 1,
            "left_steering_tape_index": 2,
            "right_steering_tape_index": 3,
        },
    }
}
config_data_kite_3 = {
    "kite": {
        "pulley": {
            "point_indices": np.array([24, 28]),
            "mass": 0.1,
            "number_of_pulleys_in_back_lines": np.array(2),
            "line_indices": [6, 4, 8, 5],
            "line_pair_indices": {"6": 4, "8": 5},
            "ci": np.array([23, 27, 22, 22]),
            "cj": np.array([24, 28, 24, 28]),
            "other_line_pair": {
                "6": np.array([22.0, 24.0, 2.86387114]),
                "4": np.array([23.0, 24.0, 2.34058161]),
                "8": np.array([22.0, 28.0, 2.86387114]),
                "5": np.array([27.0, 28.0, 2.34058161]),
            },
        }
    }
}

config_data_kite_4 = {
    "kite": {
        "kcu": {"drag_coefficient": 0.47, "diameter": 0.38, "index": 21, "mass": 8.4},
        "connectivity": {
            "bridle_ci": np.array(
                [
                    0,
                    21,
                    21,
                    21,
                    22,
                    22,
                    23,
                    23,
                    27,
                    27,
                    24,
                    24,
                    28,
                    28,
                    25,
                    25,
                    29,
                    29,
                    26,
                    26,
                    30,
                    30,
                    0,
                    0,
                    31,
                    31,
                    34,
                    34,
                    32,
                    32,
                    35,
                    35,
                    33,
                    33,
                    36,
                    36,
                    24,
                    28,
                    31,
                    34,
                ]
            ),
            "bridle_cj": np.array(
                [
                    21,
                    22,
                    23,
                    27,
                    24,
                    28,
                    24,
                    1,
                    28,
                    10,
                    25,
                    26,
                    29,
                    30,
                    18,
                    17,
                    11,
                    12,
                    16,
                    15,
                    13,
                    14,
                    31,
                    34,
                    32,
                    33,
                    35,
                    36,
                    2,
                    3,
                    9,
                    8,
                    4,
                    5,
                    7,
                    6,
                    19,
                    20,
                    19,
                    20,
                ]
            ),
            "plate_point_indices": np.array(
                [
                    [19, 2, 18, 1],
                    [2, 3, 17, 18],
                    [3, 4, 16, 17],
                    [4, 5, 15, 16],
                    [5, 6, 14, 15],
                    [6, 7, 13, 14],
                    [7, 8, 12, 13],
                    [8, 9, 11, 12],
                    [9, 20, 10, 11],
                ]
            ),
            "wing_ci": np.array(
                [
                    19,
                    2,
                    18,
                    1,
                    19,
                    2,
                    2,
                    3,
                    17,
                    18,
                    2,
                    3,
                    3,
                    4,
                    16,
                    17,
                    3,
                    4,
                    4,
                    5,
                    15,
                    16,
                    4,
                    5,
                    5,
                    6,
                    14,
                    15,
                    5,
                    6,
                    6,
                    7,
                    13,
                    14,
                    6,
                    7,
                    7,
                    8,
                    12,
                    13,
                    7,
                    8,
                    8,
                    9,
                    11,
                    12,
                    8,
                    9,
                    9,
                    20,
                    10,
                    11,
                    9,
                    20,
                ]
            ),
            "wing_cj": np.array(
                [
                    2,
                    18,
                    1,
                    19,
                    18,
                    1,
                    3,
                    17,
                    18,
                    2,
                    17,
                    18,
                    4,
                    16,
                    17,
                    3,
                    16,
                    17,
                    5,
                    15,
                    16,
                    4,
                    15,
                    16,
                    6,
                    14,
                    15,
                    5,
                    14,
                    15,
                    7,
                    13,
                    14,
                    6,
                    13,
                    14,
                    8,
                    12,
                    13,
                    7,
                    12,
                    13,
                    9,
                    11,
                    12,
                    8,
                    11,
                    12,
                    20,
                    10,
                    11,
                    9,
                    10,
                    11,
                ]
            ),
            "te_line_indices": np.array([32, 2, 38, 8, 44, 14, 50, 20, 26]),
            "tube_line_indices": np.array(
                [
                    0,
                    1,
                    3,
                    6,
                    7,
                    9,
                    12,
                    13,
                    15,
                    18,
                    19,
                    21,
                    24,
                    25,
                    27,
                    30,
                    31,
                    33,
                    36,
                    37,
                    39,
                    42,
                    43,
                    45,
                    48,
                    49,
                    51,
                ]
            ),
        },
        "airfoil": {
            "tube_diameters": np.array(
                [
                    0.1,
                    0.151561,
                    0.178254,
                    0.19406,
                    0.202418,
                    0.202418,
                    0.19406,
                    0.178254,
                    0.151561,
                    0.1,
                ]
            ),
            "is_tube_diameter_dimensionless": False,
            "canopy_max_heights": np.array(
                [0.02, 0.03, 0.04, 0.05, 0.06, 0.06, 0.05, 0.04, 0.03, 0.02]
            ),
            "is_canopy_max_height_dimensionless": True,
        },
        "stiffness": {
            "bridle": 50000.0,
            "tube": 50000.0,
            "trailing_edge": 50000.0,
            "canopy": 1000.0,
            "k_bend_strut": 0,
            "k_bend_leading_edge": 0,
        },
    },
}


def update_dict_with_instantiated_classes(
    data_dict: dict, dict_nested_instances: dict
) -> dict:
    """Update the data_dict to include the instantiated classes."""
    for key, value in data_dict.items():
        if key in dict_nested_instances:
            data_dict[key] = dict_nested_instances[key]
    return data_dict


def create_attr_class_from_dict_NEW(class_name: str, nest0_data_dict: dict) -> any:
    """Create an attrs class from a dictionary with up to three levels of nesting."""
    nest0_attributes = {}
    nest1_attributes = {}
    nest2_attributes = {}
    nest3_attributes = {}
    nest4_attributes = {}

    nest1_instances = {}
    nest2_instances = {}
    nest3_instances = {}
    nest4_instances = {}

    # Level 0
    for nest0_key, nest0_value in nest0_data_dict.items():
        print("nest0_key: ", nest0_key, type(nest0_key))
        nest0_key = str(nest0_key)
        if isinstance(nest0_value, dict):
            nest1_data_dict = nest0_value

            # Level 1
            for nest1_key, nest1_value in nest1_data_dict.items():
                if isinstance(nest1_value, dict):
                    nest2_data_dict = nest1_value

                    # Level 2
                    for nest2_key, nest2_value in nest2_data_dict.items():
                        if isinstance(nest2_value, dict):
                            nest3_data_dict = nest2_value

                            # Level 3
                            for nest3_key, nest3_value in nest3_data_dict.items():
                                if isinstance(nest3_value, dict):
                                    nest4_data_dict = nest3_value

                                    # Level 4
                                    for (
                                        nest4_key,
                                        nest4_value,
                                    ) in nest4_data_dict.items():
                                        # find the last level attributes
                                        nest4_attributes[nest4_key] = attr.ib(
                                            default=nest4_value
                                        )

                                    #### LEVEL 4
                                    # 1. Create class,
                                    Nest4Class = attr.make_class(
                                        nest3_key, nest4_attributes, frozen=True
                                    )
                                    # 2. Update the data_dict to include the instantiated classes
                                    # 3. Instantiate the class
                                    nest4_instance = Nest4Class(**nest4_data_dict)
                                    # 4. Add the instantiated class to the dict
                                    nest4_instances[nest3_key] = nest4_instance
                                    # 5. Add to the attributes
                                    nest3_attributes[nest3_key] = attr.ib(
                                        default=nest4_instance
                                    )
                                else:
                                    nest3_attributes[nest3_key] = attr.ib(
                                        default=nest3_value
                                    )

                            ### LEVEL 3
                            # 1. Create class,
                            Nest3Class = attr.make_class(
                                nest2_key, nest3_attributes, frozen=True
                            )
                            # 2. Update the data_dict to include the instantiated classes
                            nest3_data_dict = update_dict_with_instantiated_classes(
                                nest3_data_dict, nest4_instances
                            )
                            # 3. Instantiate the class
                            nest3_instance = Nest3Class(**nest3_data_dict)
                            # 4. Add the instantiated class to the dict
                            nest3_instances[nest2_key] = nest3_instance
                            # 5. Add to the attributes
                            nest2_attributes[nest2_key] = attr.ib(
                                default=nest3_instance
                            )
                        else:
                            nest2_attributes[nest2_key] = attr.ib(default=nest2_value)

                    ### LEVEL 2
                    # 1. Create class,
                    Nest2Class = attr.make_class(
                        nest1_key, nest2_attributes, frozen=True
                    )
                    # 2. Update the data_dict to include the instantiated classes
                    nest2_data_dict = update_dict_with_instantiated_classes(
                        nest2_data_dict, nest3_instances
                    )
                    # 3. Instantiate the class
                    nest2_instance = Nest2Class(**nest2_data_dict)
                    # 4. Add the instantiated class to the dict
                    nest2_instances[nest1_key] = nest2_instance
                    # 5. Add to the attributes
                    nest1_attributes[nest1_key] = attr.ib(default=nest2_instance)
                else:
                    nest1_attributes[nest1_key] = attr.ib(default=nest1_value)

            ### LEVEL 1
            # 1. Create class,
            Nest1Class = attr.make_class(nest0_key, nest1_attributes, frozen=True)
            # 2. Update the data_dict to include the instantiated classes
            nest1_data_dict = update_dict_with_instantiated_classes(
                nest1_data_dict, nest2_instances
            )
            # 3. Instantiate the class
            nest1_instance = Nest1Class(**nest1_data_dict)
            # 4. Add the instantiated class to the dict
            nest1_instances[nest0_key] = nest1_instance
            # Add to the attributes
            nest0_attributes[nest0_key] = attr.ib(default=nest1_instance)
        else:
            nest0_attributes[nest0_key] = attr.ib(default=nest0_value)

    ### LEVEL 0
    # 1. Create class,
    Nest0Class = attr.make_class(class_name, nest0_attributes, frozen=True)
    # 2. Update the data_dict to include the instantiated classes
    nest0_data_dict = update_dict_with_instantiated_classes(
        nest0_data_dict, nest1_instances
    )
    # 3. Instantiate the class
    nest0_instance = Nest0Class(**nest0_data_dict)

    return nest0_instance


# the below used to work
def create_attr_class_from_dict_infinity(class_name: str, data_dict: dict) -> any:
    """Recursively create an attrs class from a dictionary."""
    attributes = {}
    nested_attributes = {}
    nested_nested_attributes = {}
    dict_nested_instances = {}
    dict_nested_nested_instances = {}

    logging.debug("Inside the function")
    for key, value in data_dict.items():
        logging.debug(f"key: {key}, value: {value}")
        # if a nested dictionary is found, create a new class for it
        if isinstance(value, dict):
            logging.debug(" ")
            logging.debug("Nested dict is found")
            nested_data_dict = value
            nested_class_name = str(key)

            for nested_key, nested_value in nested_data_dict.items():
                logging.debug(f"nested_key: {nested_key}, nested_value: {nested_value}")

                # looping it through on itself
                nested_nested_class_instance = create_attr_class_from_dict_infinity(
                    nested_class_name, nested_data_dict
                )
                nested_attributes[nested_key] = attr.ib(
                    default=nested_nested_class_instance
                )
                dict_nested_nested_instances[nested_class_name] = (
                    nested_nested_class_instance
                )

            # (1) Create the nested class
            NestedClass = attr.make_class(
                nested_class_name, nested_attributes, frozen=True
            )
            # (2) Update the nested_data_dict to include the instantiated nested nested classes
            for nested_key, nested_value in nested_data_dict.items():
                if nested_key in dict_nested_nested_instances:
                    nested_data_dict[nested_key] = dict_nested_nested_instances[
                        nested_key
                    ]

            # (3) Instantiate the nested class
            nested_instance = NestedClass(**nested_data_dict)
            logging.debug(f"nested_instance: {nested_instance}")

            ## Needed for adding to the main-class
            dict_nested_instances[nested_class_name] = nested_instance
            attributes[key] = attr.ib(default=nested_instance)

        else:
            attributes[key] = attr.ib(default=value)

    # (1) create the class
    DictParentClass = attr.make_class(class_name, attributes, frozen=True)

    # (2) Update the data_dict to include the instantiated nested classes
    for key, value in data_dict.items():
        if key in dict_nested_instances:
            data_dict[key] = dict_nested_instances[key]

    # (3) Instantiate the class
    dict_parent_instance = DictParentClass(**data_dict)

    return dict_parent_instance


# example dict
example_dict = {
    "key1": 1,
    "key2": 2,
    "key3": 3,
    "key4": {"nested_key1": 4, "nested_key2": 5, "nested_key3": 6},
    "nest0": {
        "nest1": {
            "nested_nested_key1": 7,
            "nested_nested_key2": 8,
            "nest2": {"nest3": 9},
        }
    },
}


# ## finding level of depth
def max_nesting_level(data_dict):
    max_level = 0

    def _max_nesting_level(data_dict, current_level):
        nonlocal max_level
        if not isinstance(data_dict, dict):
            return
        max_level = max(max_level, current_level)
        for value in data_dict.values():
            _max_nesting_level(value, current_level + 1)

    _max_nesting_level(data_dict, 0)
    return max_level


def changing_dict_list_entries_to_arrays(dict_input: dict) -> dict:
    """This function loops through a dictionary and changes all the lists to np.arrays

    Args:
        dict_input (dict): dictionary with lists as values

    Returns:
        dict: dictionary with np.arrays as values
    """
    for key, value in dict_input.items():
        if isinstance(value, list):
            # and not isinstance(value, str)
            # and not isinstance(value, dict):
            dict_input[key] = np.array([float(value_i) for value_i in value])
        elif isinstance(value, dict):
            nested_dict = value
            dict_input[key] = changing_dict_list_entries_to_arrays(nested_dict)

    return dict_input


print(f"Depth of dictionairy: {max_nesting_level(config_data)}")
# config_data = changing_dict_list_entries_to_arrays(config_data)
# random_name = create_attr_class_from_dict_infinity("random_name", config_data)
random_name = create_attr_class_from_dict_NEW("random_name", config_data_kite_1)
# print(random_name)
