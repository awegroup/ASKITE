import numpy as np
import os
from pathlib import Path
import logging
import attr

from kitesim.structural import structural_mesher
from kitesim.initialisation.initialisation_utils import (
    load_module_from_path,
    calculate_projected_area,
)
from kitesim.initialisation import mass_distribution, pulley_connectivity
from kitesim.cases import cases_yaml_reader


## Initialize Kite-specific configuration
def ini_kite_config(config_data, config_data_kite, path_kite_data):

    ## Finding the Connectivity
    # Reading out surfplan
    # TODO: implement real surfplan-reader here surfplan_file = load_surfplan.fetch_surfplan_file(config.kite.surfplan_filename,sys.path) # reading out the surfplan file
    surfplan_file = None

    kite_name = config_data["kite_name"]

    current_dir = os.path.dirname(__file__)
    folder_path_initialisation_kite_specific = os.path.join(current_dir, kite_name)

    extract_points_and_connectivity = load_module_from_path(
        kite_name, f"{folder_path_initialisation_kite_specific}/analyze_surfplan.py"
    ).extract_points_and_connectivity

    (
        points_struc,
        bridle_ci,
        bridle_cj,
        plate_point_indices,
        wing_ci,
        wing_cj,
        te_line_indices,
        tube_line_indices,
    ) = extract_points_and_connectivity(path_kite_data, surfplan_file)

    # scale points if so desired
    points_ini = np.array(points_struc) / config_data["geometric_scaling_factor"]

    bridle_rest_lengths_initial = np.array(
        structural_mesher.calculate_edge_lengths(bridle_ci, bridle_cj, points_ini)
    )
    wing_rest_lengths_initial = np.array(
        structural_mesher.calculate_edge_lengths(wing_ci, wing_cj, points_ini)
    )  # [m]

    # Billowing (both TE and middle TE points to make it symmetrical), when config.is_billowing_on is True
    if config_data["is_billowing_on"]:
        points_ini, wing_rest_lengths = structural_mesher.update_for_billowing(
            points_ini, wing_rest_lengths_initial, config_data["u_p"]
        )

    # number of structural segments
    n_segments = int(len(plate_point_indices))

    ### Getting additional data
    if kite_name == "V3_25":
        # TODO: Most variables here should be generated/imported from the surfplan file instead
        BILLOWING_ANGLES = (
            np.array(config_data_kite["billowing_angles"])
            / config_data["geometric_scaling_factor"]
        )
        TUBE_DIAMETERS = (
            np.array(config_data_kite["tube_diameters"])
            / config_data["geometric_scaling_factor"]
        )
        CANOPY_MAX_HEIGHTS = (
            np.array(config_data_kite["canopy_max_heights"])
            / config_data["geometric_scaling_factor"]
        )

        ## aditional data for bridles
        bridle_data = {}
        bridle_data["diameter"] = config_data_kite["bridle"]["diameter"]
        bridle_data["density"] = config_data_kite["bridle"]["density"]
        bridle_data["bridle_point_index"] = config_data_kite["bridle"][
            "bridle_point_index"
        ]
        # element indices
        bridle_data["depower_tape_index"] = config_data_kite["bridle"][
            "depower_tape_index"
        ]
        bridle_data["left_steering_tape_index"] = config_data_kite["bridle"][
            "left_steering_tape_index"
        ]
        bridle_data["right_steering_tape_index"] = config_data_kite["bridle"][
            "right_steering_tape_index"
        ]

        ## additional data for pulleys
        pulley_data = {}
        pulley_data["point_indices"] = config_data_kite["pulley"]["point_indices"]
        pulley_data["mass"] = config_data_kite["pulley"]["mass"]
        pulley_data["number_of_pulleys_in_back_lines"] = config_data_kite["pulley"][
            "number_of_pulleys_in_back_lines"
        ]

        ## kcu data
        kcu_data = {}
        kcu_data["index"] = config_data_kite["kcu"]["index"]
        kcu_data["extra"] = {}

        ## stiffness rotational ##TODO: can this be done cleaner?
        stiffness_bend_strut = 0
        stiffness_bend_leading_edge = 0

    elif kite_name == "V9_60C":
        # TODO: Should be generated/imported directly from the surfplan file instead

        RIB_DB_WHOLE_STRUTS = np.load(
            f"{path_kite_data}/rib_db_whole_model_struts.npy",
            allow_pickle=True,
        )
        extract_airfoil_geometry_data_from_ribs = load_module_from_path(
            kite_name, f"{folder_path_initialisation_kite_specific}/analyze_surfplan.py"
        ).extract_airfoil_geometry_data_from_ribs
        print(extract_airfoil_geometry_data_from_ribs)
        (
            TUBE_DIAMETERS,
            CANOPY_MAX_HEIGHTS,
            BILLOWING_ANGLES,
        ) = extract_airfoil_geometry_data_from_ribs(
            RIB_DB_WHOLE_STRUTS,
            folder_path_initialisation_kite_specific,
            config_data["kite_name"],
        )

        ## bridle data
        bridle_data = {}
        bridle_data["diameter"] = config_data_kite["bridle"]["diameter"]
        bridle_data["density"] = config_data_kite["bridle"]["density"]
        bridle_data["bridle_point_index"] = np.load(
            f"{path_kite_data}/bridlepoint_index.npy", allow_pickle=True
        )
        # element indices
        bridle_data["depower_tape_index"] = config_data_kite["bridle"][
            "depower_tape_index"
        ]
        bridle_data["left_steering_tape_index"] = config_data_kite["bridle"][
            "left_steering_tape_index"
        ]
        bridle_data["right_steering_tape_index"] = config_data_kite["bridle"][
            "right_steering_tape_index"
        ]

        ## pulley_data
        pulley_data = {}
        pulley_data["point_indices"] = np.load(
            f"{path_kite_data}/pulley_point_indices.npy", allow_pickle=True
        )
        pulley_data["mass"] = config_data_kite["pulley"]["mass"]
        pulley_data["number_of_pulleys_in_back_lines"] = config_data_kite["pulley"][
            "number_of_pulleys_in_back_lines"
        ]

        ## kcu data
        kcu_data = {}
        kcu_data["extra"] = {
            "kcu_point_indices": np.load(
                f"{path_kite_data}/kcu_point_indices.npy", allow_pickle=True
            ),
            "kcu_line_indices": np.load(
                f"{path_kite_data}/kcu_line_indices.npy", allow_pickle=True
            ),
            "kcu_plate_indices": np.load(
                f"{path_kite_data}/kcu_plate_indices.npy", allow_pickle=True
            ),
        }

        ## rotational stiffness
        stiffness_bend_strut = config_data_kite["stiffness_bend_strut"]
        stiffness_bend_leading_edge = config_data_kite["stiffness_bend_leading_edge"]

    ## Extracting pulley connectivity
    pulley_data = pulley_connectivity.extract_pulley_connectivity(
        points_ini, bridle_ci, bridle_cj, pulley_data
    )
    mass_points = mass_distribution.new_calculate_mass_distribution(
        points_ini,
        bridle_ci,
        bridle_cj,
        wing_ci,
        wing_cj,
        config_data_kite["wing_mass"],
        config_data_kite["bridle"]["density"],
        config_data_kite["bridle"]["diameter"],
        config_data_kite["kcu"]["mass"],
        config_data_kite["kcu"]["index"],
        pulley_data["point_indices"],
        config_data_kite["pulley"]["mass"],
    )

    ##TODO: This is not the true chord, but only the chord taken from the structural discretization
    ## Calculating reference distance
    ref_chord_calculated = max(points_ini[:, 0]) - min(points_ini[:, 0])

    wing_connectivity = np.hstack((wing_ci, wing_cj))
    wing_connectivity = np.unique(wing_connectivity)
    wing_nodes = np.array([points_ini[i] for i in wing_connectivity])
    area_projected_calculated = calculate_projected_area(wing_nodes)
    area_surface_calculated = (
        config_data_kite["area_surface"] / config_data["geometric_scaling_factor"]
    )
    span_calculated = max(wing_nodes[:, 1]) - min(wing_nodes[:, 1])
    height_calculated = max(wing_nodes[:, 2]) - min(wing_nodes[:, 2])

    dict_kite_config = {
        "kite": {
            "points_ini": points_ini,
            "n_points": int(len(points_ini)),
            "surfplan_filename": config_data_kite["surfplan_filename"],
            "area_projected": area_projected_calculated,
            "area_surface": area_surface_calculated,
            "ref_chord": ref_chord_calculated,
            "span": span_calculated,
            "height": height_calculated,
            "wing_mass": config_data_kite["wing_mass"],
            "is_with_elongation_limit": config_data_kite["is_with_elongation_limit"],
            "elongation_limit": config_data_kite["elongation_limit"],
            "is_with_compression_limit": config_data_kite["is_with_compression_limit"],
            "compression_limit": config_data_kite["compression_limit"],
            "limit_stiffness_factor": config_data_kite["limit_stiffness_factor"],
            "billowing_angles": BILLOWING_ANGLES,
            "n_segments": n_segments,
            "wing_rest_lengths_initial": wing_rest_lengths_initial,
            "bridle_rest_lengths_initial": bridle_rest_lengths_initial,
            "mass_points": mass_points,
            "bridle": {
                "diameter": bridle_data["diameter"],
                "density": bridle_data["density"],
                "bridle_point_index": bridle_data["bridle_point_index"],
                "depower_tape_index": bridle_data["depower_tape_index"],
                "left_steering_tape_index": bridle_data["left_steering_tape_index"],
                "right_steering_tape_index": bridle_data["right_steering_tape_index"],
            },
            "pulley": {
                "point_indices": np.array(pulley_data["point_indices"]),
                "mass": np.array(pulley_data["mass"]),
                "number_of_pulleys_in_back_lines": np.array(
                    pulley_data["number_of_pulleys_in_back_lines"]
                ),
                "line_indices": np.array(pulley_data["line_indices"]),
                "line_pair_indices": np.array(pulley_data["line_pair_indices"]),
                "ci": np.array(pulley_data["ci"]),
                "cj": np.array(pulley_data["cj"]),
                "other_line_pair": pulley_data["other_line_pair"],
            },
            "kcu": {
                "drag_coefficient": config_data_kite["kcu"]["drag_coefficient"],
                "diameter": config_data_kite["kcu"]["diameter"],
                "index": config_data_kite["kcu"]["index"],
                "mass": config_data_kite["kcu"]["mass"],
            },
            "connectivity": {
                "bridle_ci": bridle_ci,
                "bridle_cj": bridle_cj,
                "plate_point_indices": plate_point_indices,
                "wing_ci": wing_ci,
                "wing_cj": wing_cj,
                "te_line_indices": te_line_indices,
                "tube_line_indices": tube_line_indices,
            },
            "airfoil": {
                "tube_diameters": TUBE_DIAMETERS,
                "is_tube_diameter_dimensionless": config_data_kite[
                    "is_tube_diameter_dimensionless"
                ],
                "canopy_max_heights": CANOPY_MAX_HEIGHTS,
                "is_canopy_max_height_dimensionless": config_data_kite[
                    "is_canopy_max_height_dimensionless"
                ],
            },
            "stiffness": {
                "bridle": config_data_kite["stiffness_bridle"],
                "tube": config_data_kite["stiffness_tube"],
                "trailing_edge": config_data_kite["stiffness_trailing_edge"],
                "canopy": config_data_kite["stiffness_canopy"],
                # rotational
                "k_bend_strut": stiffness_bend_strut,
                "k_bend_leading_edge": stiffness_bend_leading_edge,
            },
        }
    }
    return dict_kite_config


def update_dict_with_instantiated_classes(
    data_dict: dict, dict_nested_instances: dict
) -> dict:
    """Update the data_dict to include the instantiated classes.

    Args:
        data_dict (dict): The dictionary with the data to be used to create the class.
        dict_nested_instances (dict): The dictionary with the instantiated classes.

    Returns:
        The updated data_dict.
    """
    for key, value in data_dict.items():
        if key in dict_nested_instances:
            data_dict[key] = dict_nested_instances[key]
    return data_dict


def create_attr_class_from_dict(class_name: str, nest0_data_dict: dict):
    """Create an attrs class from a dictionary with up to three levels of nesting.

    Args:
        class_name (str): The name of the class to be created.
        nest0_data_dict (dict): The dictionary with the data to be used to create the class.

    Returns:
        The instantiated class.
    """
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


# TODO: somehow this should be possible, but how?!?
# def create_attr_class_from_dict_infinity(class_name: str, nest0_data_dict: dict):
#     """Create an attrs class from a dictionary with up to three levels of nesting."""
#     nest0_attributes = {}
#     nest1_attributes = {}

#     nest1_instances = {}
#     nest2_instances = {}

#     # Level 0
#     for nest0_key, nest0_value in nest0_data_dict.items():
#         if isinstance(nest0_value, dict):
#             nest1_data_dict = nest0_value

#             # Level 1
#             for nest1_key, nest1_value in nest1_data_dict.items():

#                 #### LEVEL 2, but infinite
#                 # 3. Instantiate the class by LOOPING IT ON ITSELF
#                 nest2_instance = create_attr_class_from_dict_NEW_infinity(
#                     nest0_key, nest1_data_dict
#                 )
#                 # 4. Add the instantiated class to the dict
#                 nest2_instances[nest0_key] = nest2_instance
#                 # 5. Add to the attributes
#                 nest1_attributes[nest1_key] = attr.ib(default=nest2_instance)

#             ### LEVEL 1
#             # 1. Create class,
#             Nest1Class = attr.make_class(nest0_key, nest1_attributes, frozen=True)
#             # 2. Update the data_dict to include the instantiated classes
#             nest1_data_dict = update_dict_with_instantiated_classes(
#                 nest1_data_dict, nest2_instances
#             )
#             # 3. Instantiate the class
#             nest1_instance = Nest1Class(**nest1_data_dict)
#             # 4. Add the instantiated class to the dict
#             nest1_instances[nest0_key] = nest1_instance
#             # Add to the attributes
#             nest0_attributes[nest0_key] = attr.ib(default=nest1_instance)
#         else:
#             nest0_attributes[nest0_key] = attr.ib(default=nest0_value)

#     ### LEVEL 0
#     # 1. Create class,
#     Nest0Class = attr.make_class(class_name, nest0_attributes, frozen=True)
#     # 2. Update the data_dict to include the instantiated classes
#     nest0_data_dict = update_dict_with_instantiated_classes(
#         nest0_data_dict, nest1_instances
#     )
#     # 3. Instantiate the class
#     nest0_instance = Nest0Class(**nest0_data_dict)

#     return nest0_instance


def setup_config(
    path_config,
    path_processed_data_folder,
):
    """Setup the configuration for the simulation.

        Steps:
        1. Load the yaml defined settings into dicts
        2. Create extra dicts for some of the child classes
        3. Create a big nested dict_config_data with all of these dicts
        4. Dynamically create a frozen attrs nested config class from this big dictionary
    Args:
        path_config (str): The path to the configuration file.
        path_processed_data_folder (str): The path to the processed data folder.

    Returns:
        The configuration class.
    """

    # Load the yaml defined settings into dicts
    config_data = cases_yaml_reader.read_yaml_file(path_config, None)
    config_data.update(cases_yaml_reader.read_yaml_file(None, "default_case.yaml"))
    config_data.update(
        cases_yaml_reader.read_yaml_file(None, f"{config_data['sim_name']}.yaml")
    )
    path_kite_config = (
        Path(path_processed_data_folder)
        / config_data["kite_name"]
        / f"config_kite_{config_data['kite_name']}.yaml"
    )
    path_kite_data = (
        Path(path_processed_data_folder)
        / str(config_data["kite_name"])
        / "processed_design_files"
    )
    config_data_kite = cases_yaml_reader.read_yaml_file(path_kite_config, None)
    dict_kite_config = ini_kite_config(
        config_data,
        config_data_kite,
        path_kite_data,
    )
    config_data.update(dict_kite_config)

    # Correcting certain dict items, ensuring the float and np.array types
    config_data["vel_wind"] = np.array([float(i) for i in config_data["vel_wind"]])
    config_data["vel_kite"] = np.array([float(i) for i in config_data["vel_kite"]])
    config_data["acc_kite"] = np.array([float(i) for i in config_data["acc_kite"]])
    config_data["plot_elev"] = np.array([float(i) for i in config_data["plot_elev"]])
    config_data["plot_azim"] = np.array([float(i) for i in config_data["plot_azim"]])

    return create_attr_class_from_dict("config", config_data)
