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

from .defining_config_class import (
    Config,
    AeroConfig,
    SolverConfig,
    AeroStructuralConfig,
    TetherConfig,
    KiteConfig,
    BridleConfig,
    PulleyConfig,
    KCUConfig,
    ConnectivityConfig,
    AirfoilGeometry,
    StiffnessConfig,
)


## Initialize kite config child-classes
def ini_bridle_config(BridleConfig, config_data_kite):
    bridle_config = BridleConfig(
        diameter=config_data_kite["bridle"]["diameter"],
        density=config_data_kite["bridle"]["density"],
        bridle_point_index=config_data_kite["bridle"]["bridle_point_index"],
        depower_tape_index=config_data_kite["bridle"]["depower_tape_index"],
        left_steering_tape_index=config_data_kite["bridle"]["left_steering_tape_index"],
        right_steering_tape_index=config_data_kite["bridle"][
            "right_steering_tape_index"
        ],
    )
    return bridle_config


def ini_pulley_config(PulleyConfig, config_data_kite):
    pulley_config = PulleyConfig(
        point_indices=config_data_kite["pulley"]["point_indices"],
        mass=config_data_kite["pulley"]["mass"],
        number_of_pulleys_in_back_lines=config_data_kite["pulley"][
            "number_of_pulleys_in_back_lines"
        ],
        line_indices=config_data_kite["pulley"]["line_indices"],
        line_pair_indices=config_data_kite["pulley"]["line_pair_indices"],
        ci=config_data_kite["pulley"]["ci"],
        cj=config_data_kite["pulley"]["cj"],
        other_line_pair=config_data_kite["pulley"]["other_line_pair"],
    )
    return pulley_config


def ini_kcu_config(KCUConfig, config_data_kite):
    kcu_config = KCUConfig(
        drag_coefficient=config_data_kite["kcu"]["drag_coefficient"],
        diameter=config_data_kite["kcu"]["diameter"],
        index=config_data_kite["kcu"]["index"],
        mass=config_data_kite["kcu"]["mass"],
    )
    return kcu_config


def ini_connectivity_config(ConnectivityConfig, config_data_kite):
    connectivity_config = ConnectivityConfig(
        bridle_ci=config_data_kite["connectivity"]["bridle_ci"],
        bridle_cj=config_data_kite["connectivity"]["bridle_cj"],
        plate_point_indices=config_data_kite["connectivity"]["plate_point_indices"],
        wing_ci=config_data_kite["connectivity"]["wing_ci"],
        wing_cj=config_data_kite["connectivity"]["wing_cj"],
        te_line_indices=config_data_kite["connectivity"]["te_line_indices"],
        tube_line_indices=config_data_kite["connectivity"]["tube_line_indices"],
    )
    return connectivity_config


def ini_airfoil_geometry_config(AirfoilGeometry, config_data_kite):
    airfoil_geometry_config = AirfoilGeometry(
        tube_diameters=config_data_kite["airfoil"]["tube_diameters"],
        is_tube_diameter_dimensionless=config_data_kite["airfoil"][
            "is_tube_diameter_dimensionless"
        ],
        canopy_max_heights=config_data_kite["airfoil"]["canopy_max_heights"],
        is_canopy_max_height_dimensionless=config_data_kite["airfoil"][
            "is_canopy_max_height_dimensionless"
        ],
    )
    return airfoil_geometry_config


def ini_stiffness_config(StiffnessConfig, config_data_kite):
    stiffness_config = StiffnessConfig(
        bridle=config_data_kite["stiffness"]["bridle"],
        tube=config_data_kite["stiffness"]["tube"],
        trailing_edge=config_data_kite["stiffness"]["trailing_edge"],
        canopy=config_data_kite["stiffness"]["canopy"],
        k_bend_strut=config_data_kite["stiffness"]["k_bend_strut"],
        k_bend_leading_edge=config_data_kite["stiffness"]["k_bend_leading_edge"],
    )
    return stiffness_config


## Initialize child-classes first
def ini_aero_config(AeroConfig, config_data):
    # Initialize Aero configuration
    aero_config = AeroConfig(
        ring_geometry=config_data["AeroConfig"]["ring_geometry"],
        model=config_data["AeroConfig"]["model"],
        n_iter=config_data["AeroConfig"]["n_iter"],
        error_limit=config_data["AeroConfig"]["error_limit"],
        relax_factor=config_data["AeroConfig"]["relax_factor"],
        n_splits=config_data["AeroConfig"]["n_splits"],
        n_chordwise_aero_nodes=config_data["AeroConfig"]["n_chordwise_aero_nodes"],
        cd_cylinder=config_data["AeroConfig"]["cd_cylinder"],
        cd_shear_cylinder=config_data["AeroConfig"]["cd_shear_cylinder"],
        cd_multiplier=config_data["AeroConfig"]["cd_multiplier"],
    )
    return aero_config


def ini_solver_config(SolverConfig, config_data):
    solver_config = SolverConfig(
        method=config_data["SolverConfig"]["method"],
        tol=config_data["SolverConfig"]["tol"],
        max_fev=config_data["SolverConfig"]["max_fev"],
        integration_method=config_data["SolverConfig"]["integration_method"],
        newton_rtol=config_data["SolverConfig"]["newton_rtol"],
        newton_max_iter=config_data["SolverConfig"]["newton_max_iter"],
        newton_disp=config_data["SolverConfig"]["newton_disp"],
        damping_constant=config_data["SolverConfig"]["damping_constant"],
        is_with_visc_damping=config_data["SolverConfig"]["is_with_visc_damping"],
        dt=config_data["SolverConfig"]["dt"],
        n_time_steps=config_data["SolverConfig"]["n_time_steps"],
        abs_tol=config_data["SolverConfig"]["abs_tol"],
        rel_tol=config_data["SolverConfig"]["rel_tol"],
        max_iter=config_data["SolverConfig"]["max_iter"],
    )
    return solver_config


def ini_aero_structural_config(AeroStructuralConfig, config_data):
    aero_structural_config = AeroStructuralConfig(
        max_iter=config_data["AeroStructuralConfig"]["max_iter"],
        it_check=config_data["AeroStructuralConfig"]["it_check"],
        tol=config_data["AeroStructuralConfig"]["tol"],
        crosswind_max_iter=config_data["AeroStructuralConfig"]["crosswind_max_iter"],
        crosswind_tol=config_data["AeroStructuralConfig"]["crosswind_tol"],
        crosswind_relax_factor=config_data["AeroStructuralConfig"][
            "crosswind_relax_factor"
        ],
    )
    return aero_structural_config


def ini_tether_config(TetherConfig, config_data):
    tether_config = TetherConfig(
        diameter=config_data["TetherConfig"]["diameter"],
        length=config_data["TetherConfig"]["length"],
        density=config_data["TetherConfig"]["density"],
    )
    return tether_config


## Initialize Kite-specific configuration
def ini_kite_config(KiteConfig, config_data, config_data_kite, path_kite_data):

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


def instantiate_kite_config(
    KiteConfig,
    BridleConfig,
    PulleyConfig,
    KCUConfig,
    ConnectivityConfig,
    AirfoilGeometry,
    StiffnessConfig,
    dict_kite_config,
):
    kite_config = KiteConfig(
        points_ini=dict_kite_config["points_ini"],
        n_points=dict_kite_config["n_points"],
        surfplan_filename=dict_kite_config["surfplan_filename"],
        area_projected=dict_kite_config["area_projected"],
        area_surface=dict_kite_config["area_surface"],
        ref_chord=dict_kite_config["ref_chord"],
        span=dict_kite_config["span"],
        height=dict_kite_config["height"],
        wing_mass=dict_kite_config["wing_mass"],
        is_with_elongation_limit=dict_kite_config["is_with_elongation_limit"],
        elongation_limit=dict_kite_config["elongation_limit"],
        is_with_compression_limit=dict_kite_config["is_with_compression_limit"],
        compression_limit=dict_kite_config["compression_limit"],
        limit_stiffness_factor=dict_kite_config["limit_stiffness_factor"],
        billowing_angles=dict_kite_config["billowing_angles"],
        n_segments=dict_kite_config["n_segments"],
        wing_rest_lengths_initial=dict_kite_config["wing_rest_lengths_initial"],
        bridle_rest_lengths_initial=dict_kite_config["bridle_rest_lengths_initial"],
        mass_points=dict_kite_config["mass_points"],
        bridle=ini_bridle_config(BridleConfig, dict_kite_config),
        pulley=ini_pulley_config(PulleyConfig, dict_kite_config),
        kcu=ini_kcu_config(KCUConfig, dict_kite_config),
        connectivity=ini_connectivity_config(ConnectivityConfig, dict_kite_config),
        airfoil=ini_airfoil_geometry_config(AirfoilGeometry, dict_kite_config),
        stiffness=ini_stiffness_config(StiffnessConfig, dict_kite_config),
    )
    return kite_config


# # Loading the yaml config file
# with open(("settings/config.yaml"), "r") as config_file:
#     config_data = yaml.load(config_file, Loader=yaml.SafeLoader)


def ini_config(
    Config,
    config_data,
    aero_config,
    solver_config,
    aero_structural_config,
    tether_config,
    kite_config,
):
    # Initialize the overarching Config object
    config = Config(
        # KITE NAME
        kite_name=config_data["kite_name"],
        # CASE SETTINGS
        sim_name=config_data["sim_name"],
        is_with_vk_optimization=config_data["is_with_vk_optimization"],
        is_circular_case=config_data["is_circular_case"],
        is_run_only_1_time_step=config_data["is_run_only_1_time_step"],
        is_print_intermediate_results=config_data["is_print_intermediate_results"],
        is_with_gravity=config_data["is_with_gravity"],
        is_with_velocity_initialization=config_data["is_with_velocity_initialization"],
        # INFLOW CONDITIONS
        vel_wind=np.array([float(i) for i in config_data["vel_wind"]]),
        vel_kite=np.array([float(i) for i in config_data["vel_kite"]]),
        acc_kite=np.array([float(i) for i in config_data["acc_kite"]]),
        # ACTUATION
        u_p=config_data["u_p"],
        depower_tape_extension_percentage=config_data[
            "depower_tape_extension_percentage"
        ],
        delta_ls_ini=config_data["delta_ls_ini"],
        depower_tape_extension_step=config_data["depower_tape_extension_step"],
        depower_tape_final_extension=config_data["depower_tape_final_extension"],
        steering_tape_extension_step=config_data["steering_tape_extension_step"],
        steering_tape_final_extension=config_data["steering_tape_final_extension"],
        # OUTPUT SETTINGS
        is_with_printing=config_data["is_with_printing"],
        is_print_mid_results=config_data["is_print_mid_results"],
        is_with_initial_plot=config_data["is_with_initial_plot"],
        is_with_initial_point_velocity=config_data["is_with_initial_point_velocity"],
        is_with_plotly_plot=config_data["is_with_plotly_plot"],
        is_with_aero_geometry=config_data["is_with_aero_geometry"],
        # PLOTTING SETTINGS
        is_with_plotting=config_data["is_with_plotting"],
        plot_format=config_data["plot_format"],
        plot_elev=np.array([float(i) for i in config_data["plot_elev"]]),
        plot_azim=np.array([float(i) for i in config_data["plot_azim"]]),
        # ANIMATION SETTINGS
        is_with_animation=config_data["is_with_animation"],
        animation_fps=config_data["animation_fps"],
        animation_dpi=config_data["animation_dpi"],
        animation_bitrate=config_data["animation_bitrate"],
        animation_elev=config_data["animation_elev"],
        animation_azim=config_data["animation_azim"],
        # SIMULATION SETTINGS
        ## initialisation
        bridle_initial_compression_factor=config_data[
            "bridle_initial_compression_factor"
        ],
        geometric_scaling_factor=config_data["geometric_scaling_factor"],
        n_vel_initialisation_steps=config_data["n_vel_initialisation_steps"],
        ## physics
        is_billowing_on=config_data["is_billowing_on"],
        is_with_aero_bridle=config_data["is_with_aero_bridle"],
        is_with_aero_tether=config_data["is_with_aero_tether"],
        coupling_method=config_data["coupling_method"],
        # CASES
        ## reel-in symmetric straight
        vel_app_initial=np.array(config_data["vel_app_initial"]),
        ## crosswind flight settings
        tol_fx_ratio_to_fz=config_data["tol_fx_ratio_to_fz"],
        tol_vk_optimization=config_data["tol_vk_optimization"],
        vk_x_initial_guess_factor_of_vw=config_data["vk_x_initial_guess_factor_of_vw"],
        ## circular flight settings
        is_with_varying_va=config_data["is_with_varying_va"],
        r_0_initial=config_data["r_0_initial"],
        # PHYSICAL CONSTANTS
        grav_constant=np.array(config_data["grav_constant"]),
        rho=config_data["rho"],
        mu=config_data["mu"],
        # CHILD CLASSES
        aero=aero_config,
        solver=solver_config,
        aero_structural=aero_structural_config,
        tether=tether_config,
        kite=kite_config,
    )
    return config


# def create_attr_class_from_dict(class_name: str, data_dict: dict):
#     """Recursively create an attrs class from a dictionary."""
#     attributes = {}
#     child_classses = []
#     for key, value in data_dict.items():
#         # if value is a dictionary, that means we need to create a nested attrs class
#         if isinstance(value, dict):
#             # Use the key, that indicates the nest, as the class name
#             child_class_name = f"{key}"
#             child_attributes = {}
#             # loop through all the nested dict items and create attributes
#             for key, value in value.items():
#                 child_attributes[key] = attr.ib(default=value)
#             # make the class
#             ChildClass = attr.make_class(
#                 child_class_name, child_attributes, frozen=True
#             )

#             # Create class object
#             child_class_object = ChildClass(**value)

#             # # Initialize classes with the dictionaries
#             general_configg = General_Config(**config_data)
#             # kite_configg = Kite_Config(**dict_kite_config)

#             # # Recursively create a nested attrs class for nested dictionaries
#             # nested_class = create_attr_class_from_dict(
#             #     f"{class_name}_{key.capitalize()}", value
#             # )
#             # attributes[key] = attr.ib(default=nested_class(**value))
#         else:
#             attributes[key] = attr.ib(default=value)

#     # append the child classes to the attributes
#     for child_class in child_classses:
#         attributes[child_class.__name__] = attr.ib(default=child_class)

#     return attr.make_class(class_name, attributes, frozen=True)


def create_attr_class_from_dict(class_name: str, data_dict: dict):
    """Recursively create an attrs class from a dictionary."""
    attributes = {}
    child_classes = []
    for key, value in data_dict.items():
        # if nested item in dictionary, create a nested class
        if isinstance(value, dict):
            # Use the key to name the child class
            child_class_name = f"{class_name}_{key.capitalize()}"
            child_attributes = {}
            # Loop through the nested dictionary items and create attributes
            for nested_key, nested_value in value.items():
                if isinstance(nested_value, dict):
                    # If the nested value is itself a dictionary, recursively create a nested class
                    nested_class = create_attr_class_from_dict(
                        f"{child_class_name}_{nested_key.capitalize()}", nested_value
                    )
                    child_attributes[nested_key] = attr.ib(
                        default=nested_class(**nested_value)
                    )
                else:
                    child_attributes[nested_key] = attr.ib(default=nested_value)

            # Make the child class
            ChildClass = attr.make_class(
                child_class_name, child_attributes, frozen=True
            )
            # Append the child class to the list of child classes
            child_classes.append(ChildClass)
        else:
            attributes[key] = attr.ib(default=value)

    # Append the child classes to the attributes
    for child_class in child_classes:
        attributes[child_class.__name__] = attr.ib(default=child_class)

    # Create the main class
    return attr.make_class(class_name, attributes, frozen=True)


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


def setup_config(
    path_config,
    path_processed_data_folder,
):
    # Define paths

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
        KiteConfig,
        config_data,
        config_data_kite,
        path_kite_data,
    )

    ###############
    #  OLD METHOD
    ###############
    aero_config = ini_aero_config(AeroConfig, config_data)
    solver_config = ini_solver_config(SolverConfig, config_data)
    aero_structural_config = ini_aero_structural_config(
        AeroStructuralConfig, config_data
    )
    tether_config = ini_tether_config(TetherConfig, config_data)

    kite_config = instantiate_kite_config(
        KiteConfig,
        BridleConfig,
        PulleyConfig,
        KCUConfig,
        ConnectivityConfig,
        AirfoilGeometry,
        StiffnessConfig,
        dict_kite_config["kite"],
    )

    # Create Config
    config = ini_config(
        Config,
        config_data,
        aero_config,
        solver_config,
        aero_structural_config,
        tether_config,
        kite_config,
    )

    ##############
    # NEW METHOD

    # (1) create a dict from the yaml files
    # (2) create extra dicts for some of the child classes
    # (3) create a big nested dict_config_data with all of these dicts
    # (4) dynamically create a frozen attrs nested config class from this big dictionary

    ##############

    # TODO: change yaml from AeroConfig to just aero
    ## adding the dicts to the config_data without changing the yaml for now
    config_data["aero"] = config_data["AeroConfig"]
    config_data["solver"] = config_data["SolverConfig"]
    config_data["aero_structural"] = config_data["AeroStructuralConfig"]
    config_data["tether"] = config_data["TetherConfig"]
    config_data["kite"] = dict_kite_config["kite"]
    # print(f'config_data["TetherConfig"]: {config_data["TetherConfig"]}')
    # print(f'dict_kite_config["kite"]: {dict_kite_config["kite"]}')

    print(config_data)
    general_configg = create_attr_class_from_dict_infinity(
        "general_configg", config_data
    )

    # Adding the kite
    # config_data.update(dict_kite_config)

    ## looping through config_data to find a list or array, and transform it to a np.array
    def changing_dict_list_entries_to_arrays(dict_input: dict) -> dict:
        """This function loops through a dictionary and changes all the lists to np.arrays

        Args:
            dict_input (dict): dictionary with lists as values

        Returns:
            dict: dictionary with np.arrays as values
        """
        for key, value in dict_input.items():
            if isinstance(value, list) 
                and not isinstance(value, str) 
                and not isinstance(value, dict):
                dict_input[key] = np.array([float(value_i) for value_i in value])
            elif isinstance(value, dict):
                nested_dict = value
                dict_input[key] = changing_dict_list_entries_to_arrays(nested_dict)

        return dict_input

    config_data = changing_dict_list_entries_to_arrays(config_data)

    ## finding level of depth
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

    print(f"Depth of dictionairy: {max_nesting_level(config_data)}")

    # general_configg = create_attr_class_from_dict_infinity("general_configg", config_data)

    # ##################
    # # PRINTING DIFFERENCES
    # ##################

    print("---------------------------------")
    print("Differences between General_Config and config:")
    print(set(attr.asdict(general_configg)) - set(attr.asdict(config)))
    # print(set(config_data) - set(attr.asdict(config)))
    print(f" ")
    print(set(attr.asdict(config)) - set(attr.asdict(general_configg)))
    # print(set(attr.asdict(config)) - set(config_data))
    print("---------------------------------")

    # print(f"config.solver: {(config.solver)}")
    # print(f"general_configg.SolverConfig: {(general_configg.solver)}")
    # print(f'config.solver - general_configg.solver: {set(attr.asdict(  config.solver)) - set(attr.asdict(general_configg.SolverConfig))}')

    return config
