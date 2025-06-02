"""
### Info

Author: Jelle Poland \
Citing: https://doi.org/10.3390/en16145264 \
License: ... \
Github: ...
"""

### Initialisation

# Making things autoreload - needed for Jupyter Kernel/Interactive env.
# %load_ext autoreload
# %autoreload 2
# %matplotlib widget

import os
from pathlib import Path
import time
import numpy as np
import pandas as pd
import logging
import json

from kitesim.initialisation import initialisation_main
from kitesim.logging_config import *
from kitesim.VortexStepMethod import VSM_old
from kitesim.aerodynamic import bridle_line_system_aero
from kitesim.coupling import coupling_struc2aero, coupling_aero2struc
from kitesim.solver import solver_utils
from kitesim.aerodynamic import tether_aero, bridle_line_system_aero
from kitesim.VortexStepMethod import VSM_old
from PSS.particleSystem.SpringDamper import SpringDamperType
import PSS.particleSystem as PSS


def generating_PSS_test_input(sim_input, path_results_folder: str):
    """Runs the aero-structural solver for the given input parameters"""

    # Unpacking input_dict
    points = sim_input["points"]
    vel_app = sim_input["vel_app"]
    config = sim_input["config"]
    input_bridle_aero = sim_input["input_bridle_aero"]
    input_tether_aero = sim_input["input_tether_aero"]
    input_VSM = sim_input["input_VSM"]
    input_PSM = sim_input["input_PSM"]

    # GENERAL INITIALIsATION
    ## case settings
    sim_name = config.sim_name
    is_with_vk_optimization = config.is_with_vk_optimization
    is_circular_case = config.is_circular_case
    is_run_only_1_time_step = config.is_run_only_1_time_step
    is_print_intermediate_results = config.is_print_intermediate_results
    is_with_gravity = config.is_with_gravity
    is_with_velocity_initialisation = config.is_with_velocity_initialisation

    logging.debug(
        f"Running aero-structural simulation for {sim_name}, shape: {np.shape(sim_name)}"
    )
    logging.debug(
        f" type of input_PSM: {type(input_PSM)}, shape: {np.shape(input_PSM)}"
    )
    logging.debug(
        f" type of input_VSM: {type(input_VSM)}, shape: {np.shape(input_VSM)}"
    )
    logging.debug(
        f" type of input_bridle_aero: {type(input_bridle_aero)}, shape: {np.shape(input_bridle_aero)}"
    )
    logging.debug(
        f" type of input_tether_aero: {type(input_tether_aero)}, shape: {np.shape(input_tether_aero)}"
    )
    logging.debug(f" type of points: {type(points)}, shape: {np.shape(points)}")
    logging.debug(f" type of vel_app: {type(vel_app)}, shape: {np.shape(vel_app)}")
    logging.debug(f" type of config: {type(config)}, shape: {np.shape(config)}")
    logging.debug(
        f" type of sim_input: {type(sim_input)}, shape: {np.shape(sim_input)}"
    )
    logging.info(
        f" type of connectivity_matrix: {type(input_PSM['connectivity_matrix'])}, shape: {np.shape(input_PSM['connectivity_matrix'])}"
    )
    logging.debug(
        f" type of initial_conditions: {type(input_PSM['initial_conditions'])}, len: {len(input_PSM['initial_conditions'])}"
    )
    logging.debug(
        f" type of params_dict: {type(input_PSM['params_dict'])}, len: {len(input_PSM['params_dict'])}"
    )
    logging.debug(f"params: {input_PSM['params_dict']}")
    logging.debug(f"initial_conditions: {input_PSM['initial_conditions']}")
    logging.debug(f"connectivity_matrix: {input_PSM['connectivity_matrix']}")

    # defining the new PSS inputs
    PSS_initial_conditions = input_PSM["initial_conditions"]
    logging.debug(f"PSS_initial_conditions: {PSS_initial_conditions}")

    # refining the params_dict
    PSS_params_dict = input_PSM["params_dict"]
    PSS_params_dict.update({"g": -config.grav_constant[2]})
    logging.debug(f"PSS_params_dict: {PSS_params_dict.keys()}")

    # restructuring connectivity matrix
    connectivity_matrix = input_PSM["connectivity_matrix"]

    ## converting the string spring types to the enum
    def string_to_springdampertype(link_type: str) -> SpringDamperType:
        """
        Convert a string representation of a link type to a SpringDamperType enum value.

        Args:
            link_type (str): String representation of the link type.

        Returns:
            SpringDamperType: Corresponding enum value.

        Raises:
            ValueError: If the input string doesn't match any SpringDamperType.
        """
        try:
            return SpringDamperType(link_type.lower())
        except ValueError:
            raise ValueError(f"Invalid link type: {link_type}")

    PSS_connectivity_matrix = []
    for idx, _ in enumerate(connectivity_matrix):
        if (
            PSS_params_dict["is_compression"][idx]
            and PSS_params_dict["is_tension"][idx]
        ):
            linktype = "default"
        elif (
            PSS_params_dict["is_compression"][idx]
            and not PSS_params_dict["is_tension"][idx]
        ):
            linktype = "nontensile"
        elif PSS_params_dict["is_pulley"][idx]:
            linktype = "pulley"
        elif (
            PSS_params_dict["is_tension"][idx]
            and not PSS_params_dict["is_compression"][idx]
        ):
            linktype = "noncompressive"

        logging.debug(f"idx: {idx}")
        logging.debug(f"connectivity_matrix[idx]: {connectivity_matrix[idx]}")
        logging.debug(f"PSS_params_dict['k'][idx]: {PSS_params_dict['k'][idx]}")
        logging.debug(f"PSS_params_dict['c']: {PSS_params_dict['c']}")
        logging.debug(f"linktype: {linktype}")

        PSS_connectivity_matrix.append(
            [
                int(connectivity_matrix[idx][0]),
                int(connectivity_matrix[idx][1]),
                float(PSS_params_dict["k"][idx]),
                float(PSS_params_dict["c"]),
                linktype,
            ]
        )

    logging.debug(f"PSS_connectivity_matrix: {PSS_connectivity_matrix}")

    ## Instating the PSM from LightSailSim, PSS
    psystem = PSS.ParticleSystem(
        PSS_connectivity_matrix,
        PSS_initial_conditions,
        PSS_params_dict,
    )

    # ## instantiating the PSM
    # psystem = ParticleSystem(
    #     input_PSM["connectivity_matrix"],
    #     input_PSM["initial_conditions"],
    #     input_PSM["params_dict"],
    # )

    # TODO: would be cleaner if this was extracted from input_PSM or config
    params = input_PSM["params_dict"]

    ## setting up the position-dataframe
    t_vector = np.linspace(
        params["dt"], params["t_steps"] * params["dt"], params["t_steps"]
    )
    x = {}
    for i in range(params["n"]):
        x[f"x{i + 1}"] = np.zeros(len(t_vector))
        x[f"y{i + 1}"] = np.zeros(len(t_vector))
        x[f"z{i + 1}"] = np.zeros(len(t_vector))

    position = pd.DataFrame(index=t_vector, columns=x)
    aero_structural_tol = params["aerostructural_tol"]

    # INITIALISATION OF VARIABLES
    ## parameters used in the loop
    start_time = time.time()
    is_convergence = False
    residual_f_list = []
    f_tether_drag = np.zeros(3)
    is_residual_below_tol = False

    ##TODO: should this be done differently?
    # Defining variables outside the loop, to speed up RUNTIME
    damping_ratio = config.solver.damping_constant
    connectivity = config.kite.connectivity
    n_segments = config.kite.n_segments
    is_with_varying_va = config.is_with_varying_va
    coupling_method = config.coupling_method
    n_chordwise_aero_nodes = config.aero.n_chordwise_aero_nodes
    aero_structural_max_iter = config.aero_structural.max_iter

    ## actuation
    ### depower tape
    index_fixed_node = config.kite.bridle.bridle_point_index
    len_wing_rest_length = len(config.kite.wing_rest_lengths_initial)
    index_depower_tape = len_wing_rest_length + config.kite.bridle.depower_tape_index
    initial_length_depower_tape = params["l0"][index_depower_tape]
    depower_tape_extension_step = config.depower_tape_extension_step
    depower_tape_final_extension = config.depower_tape_final_extension
    ### steering tape
    index_steering_tape_left = (
        len_wing_rest_length + config.kite.bridle.left_steering_tape_index
    )
    index_steering_tape_right = (
        len_wing_rest_length + config.kite.bridle.right_steering_tape_index
    )
    initial_length_steering_tape_right = params["l0"][index_steering_tape_right]
    steering_tape_extension_step = config.steering_tape_extension_step
    steering_tape_final_extension = config.steering_tape_final_extension
    ### printing initial lengths
    print(
        f"Initial depower tape length: {psystem.extract_rest_length[index_depower_tape]:.3f}m"
    )
    print(
        f"Desired depower tape length: {initial_length_depower_tape + depower_tape_final_extension:.3f}m"
    )
    print(
        f"Initial steering tape length right: {psystem.extract_rest_length[index_steering_tape_right]:.3f}m"
    )
    print(
        f"Desired steering tape length right: {initial_length_steering_tape_right + steering_tape_final_extension:.3f}m"
    )

    ## gravity
    mass_points = config.kite.mass_points
    if is_with_gravity:
        force_gravity = np.array(
            [np.array(config.grav_constant) * m_pt for m_pt in mass_points]
        )
    else:
        force_gravity = np.zeros(points.shape)

    ## velocity initialisation
    n_vel_initialisation_steps = config.n_vel_initialisation_steps
    if is_with_velocity_initialisation:
        # set the velocity to a lower value, to initialize the simulation and get rid of velocity induces shock
        vel_app_linspace = np.linspace(
            config.vel_app_initial, np.copy(vel_app), n_vel_initialisation_steps
        )
        vel_app = np.copy(config.vel_app_initial)

    if is_with_vk_optimization:
        # tolerance for fx, if fx>tol*fz, then update vk_x
        tol_fx_ratio_to_fz = config.tol_fx_ratio_to_fz
        tol_vk_optimization = config.tol_vk_optimization
        vk_x_initial_guess = config.vk_x_initial_guess_factor_of_vw * config.vel_wind[2]
        pre_simulation_of_vk = solver_utils.optimalisation_of_vk_for_fx_0(
            vk_x_initial_guess,
            vel_app,
            points,
            config,
            input_VSM,
            input_bridle_aero,
            tol_vk_optimization,
        )
        vel_app[0] = -pre_simulation_of_vk
        force_aero = np.zeros(points.shape)

    ## circular initialisation
    force_centrifugal = np.zeros(points.shape)
    vel_app_distributed = None
    if is_circular_case:
        r_0 = config.r_0_initial
        length_tether = config.tether.length

        # determining sign of the centrifugal force
        ## steering_tape_final_extension is the right-one
        ## contraction is right-turn
        if steering_tape_final_extension < 0:
            centrifugal_force_sign = 1
        ## extension is left-turn
        elif steering_tape_final_extension > 0:
            centrifugal_force_sign = -1
        else:
            raise ValueError(
                "steering_tape_final_extension is zero, this is not allowed in a circular case"
            )

    ## coupling initialisation
    (
        points_wing_segment_corners_aero_orderded,
        index_transformation_struc_to_aero,
    ) = coupling_struc2aero.extract_wingpanel_corners_aero_orderded(
        points, config.kite.connectivity.plate_point_indices
    )
    ### creating a dict with key value pairs, to transform from aero to struc
    index_transformation_aero_to_struc_dict = {}
    for i, value in enumerate(index_transformation_struc_to_aero):
        index_transformation_aero_to_struc_dict[value] = i

    print(f" ")
    print(f"Running aero-structural simulation")
    print(f"----------------------------------- ")

    ######################################################################
    # SIMULATION LOOP
    ######################################################################
    ## propagating the simulation for each timestep and saving results
    for i, step in enumerate(t_vector):
        if is_with_vk_optimization:
            # TODO: the func. below re-runs, the aero analysis. can be done more efficient
            # Updating the velocity of the kite in the x-direction, to make Fx go to zero
            # run if fx > tol_fx_ratio_to_fz*Fz
            force_resultant_x = solver_utils.calculate_fx(
                -vel_app[0],
                vel_app,
                points,
                connectivity,  #
                input_tether_aero,  #
                input_VSM,
                input_bridle_aero,
            )

            f_tether_drag_x = tether_aero.calculate_tether_drag(
                input_tether_aero,
                vel_app,
            )
            f_tether_drag = np.array([f_tether_drag_x, 0, 0])

            if np.abs(force_resultant_x) > tol_fx_ratio_to_fz * np.abs(
                np.sum(force_aero[:, 2])
            ):
                print(
                    f"--- updating vel_kite_x, fx is too large: {force_resultant_x:.3f}N"
                )
                vel_app[0] = -solver_utils.optimalisation_of_vk_for_fx_0(
                    -vel_app[0],
                    vel_app,
                    points,
                    n_segments,
                    connectivity,
                    input_VSM,
                    input_bridle_aero,
                    tol_vk_optimization,
                )

        ## external force
        begin_time_f_ext = time.time()

        # Centrifugal
        if is_circular_case:
            (force_centrifugal, r_0, v_k_i_array, center_of_gravity, r_k) = (
                solver_utils.calculate_force_centrifugal_distribution_with_tether_tension(
                    r_0,
                    force_aero,
                    mass_points,
                    vel_app,
                    points,
                    length_tether,
                    is_print_intermediate_results,
                )
            )
            force_centrifugal *= centrifugal_force_sign

            ##TODO: remove hard-coding ONLY apply non-zero Fc after 25 iterations
            n_turn_initialisation_iterations = 100
            if i < n_turn_initialisation_iterations:
                force_centrifugal = np.zeros(points.shape)

            # Computing the apparant wind speed distribution
            # Using the turning-radius ratios of the aerodynamic panels to kite
            # v_a_i = v_w - v_k_i | v_k_i = v_k * (r_0 + y_i) / r_k
            # The 1/4chord control point locations are used as y_i,
            # Obtained from the controlpoints dict, only defined after first iteration
            if i > n_turn_initialisation_iterations and is_with_varying_va:
                vel_kite = input_bridle_aero.vel_wind - vel_app
                vel_app_distributed = [
                    input_bridle_aero.vel_wind
                    - (vel_kite * (r_0 + cp["coordinates"][1]) / r_k)
                    for cp in controlpoints
                ]

        # Struc --> aero
        points_wing_segment_corners_aero_orderded = points[
            index_transformation_struc_to_aero
        ]
        # Wing Aerodynamic
        (
            force_aero_wing_VSM,
            moment_aero_wing_VSM,
            F_rel,
            ringvec,
            controlpoints,
            wingpanels,
            rings,
            coord_L,
            coord_refined,
        ) = VSM_old.calculate_force_aero_wing_VSM(
            points_wing_segment_corners_aero_orderded,
            vel_app,
            input_VSM,
            vel_app_distributed,
        )
        # Aero --> struc
        if coupling_method == "NN":
            force_aero_wing = coupling_aero2struc.aero2struc_NN(
                n_chordwise_aero_nodes,
                wingpanels,
                force_aero_wing_VSM,
                points_wing_segment_corners_aero_orderded,
                index_transformation_aero_to_struc_dict,
                points,
            )
        elif coupling_method == "MSc_Oriol":
            force_aero_wing = coupling_aero2struc.aero2struc(
                points,
                connectivity.wing_ci,
                connectivity.wing_cj,
                connectivity.plate_point_indices,
                force_aero_wing_VSM,
                moment_aero_wing_VSM,
                ringvec,
                controlpoints,
            )
        else:
            raise ValueError("Coupling method not recognized; wrong name or typo")

        # Bridle Aerodynamics
        if input_bridle_aero.is_with_aero_bridle:
            force_aero_bridle = (
                bridle_line_system_aero.calculate_force_aero_bridle_thedens2022(
                    points, vel_app, input_bridle_aero
                )
            )
        else:
            force_aero_bridle = [0]
        force_aero = force_aero_wing + force_aero_bridle

        ## summing up
        force_external = force_aero + force_gravity + force_centrifugal
        # force_external[0:3] += f_tether_drag  # TODO: remove hard-code?
        ### f_external is flat, and force_external is 2D, ##TODO: could add this to the name?
        f_external = force_external.flatten()
        end_time_f_ext = time.time()

        begin_time_f_int = time.time()

        logging.debug(
            f"f_external[2]: {f_external[2]} shape of f_external: {np.shape(np.array(f_external))}, len: {len(f_external)}"
        )
        ######### running PS simulation
        # # saving external force
        # import dill

        # path_saved_folder = "/home/jellepoland/ownCloud/phd/code/kitesim/results/V3_25/2024_08_13_17h/straight/input"
        # to_be_saved_dict = {"f_external": f_external}
        # for name, data in to_be_saved_dict.items():
        #     with open(f"{path_saved_folder}/{name}.pkl", "wb") as f:
        #         dill.dump(data, f)

        ### saving to json
        class NumpyEncoder(json.JSONEncoder):
            def default(self, obj):
                if isinstance(obj, np.ndarray):
                    return obj.tolist()
                if isinstance(obj, np.integer):
                    return int(obj)
                if isinstance(obj, np.floating):
                    return float(obj)
                return super(NumpyEncoder, self).default(obj)

        def save_data_to_json(data_dict: dict, folder_path: str):
            """Save dictionary data to JSON files.

            Args:
                data_dict (dict): Dictionary containing data to be saved.
                folder_path (Path): Path to the folder where JSON files will be saved.
            """
            folder_path.mkdir(parents=True, exist_ok=True)

            for name, data in data_dict.items():
                file_path = folder_path / f"{name}.json"
                with open(file_path, "w") as f:
                    json.dump(data, f, cls=NumpyEncoder)
                logging.info(f"Saved {name} to {file_path}")

        to_be_saved_dict = {
            "connectivity_matrix": PSS_connectivity_matrix,
            "initial_conditions": PSS_initial_conditions,
            "params": PSS_params_dict,
            "f_external": f_external,
        }
        save_data_to_json(to_be_saved_dict, path_results_folder)
        break
    return


# Import modules
def main():
    """Main function"""
    # Find the root directory of the repository
    root_dir = os.path.abspath(os.path.dirname(__file__))
    while not os.path.isfile(os.path.join(root_dir, ".gitignore")):
        root_dir = os.path.abspath(os.path.join(root_dir, ".."))
        if root_dir == "/":
            raise FileNotFoundError(
                "Could not find the root directory of the repository."
            )
    # defining paths
    # path_config = "../data/config.yaml"
    path_config = Path(root_dir) / "data" / "config.yaml"
    # underlying mechanism assumes specific folder structure inside processed_data
    ## kite config files in folder: processed_data/kite_name
    ## kite data files in folder: processed_data/kite_name/processed_design_files
    path_processed_data_folder = Path(root_dir) / "processed_data"
    path_results_folder = Path(root_dir) / "processed_data" / "V3_25" / "PSS_test_input"

    # Loading the input
    sim_input = initialisation_main.get_sim_input(
        path_config,
        path_processed_data_folder,
    )

    # Generating the PSS_test_results
    generating_PSS_test_input(sim_input, path_results_folder)


if __name__ == "__main__":
    main()
