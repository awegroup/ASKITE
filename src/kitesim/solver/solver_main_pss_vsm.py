# %% Defing the spring-damper system function
# defining function for the spring-damper system

import time
import numpy as np
import pandas as pd
import logging
import os
import sys
from pathlib import Path
from kitesim.coupling import coupling_struc2aero, coupling_aero2struc
from kitesim.solver import solver_utils
from kitesim.solver.vsm_functions import initialize_vsm, run_vsm_package
from kitesim.solver.pss_functions import instantiate_psystem, run_pss


def run_aerostructural_solver(config, config_kite, PROJECT_DIR):
    """Runs the aero-structural solver for the given input parameters"""

    # # Unpacking input_dict
    # points = sim_input["points"]

    # vel_app = sim_input["vel_app"]
    # config = sim_input["config"]
    # input_bridle_aero = sim_input["input_bridle_aero"]
    # input_tether_aero = sim_input["input_tether_aero"]
    # input_VSM = sim_input["input_VSM"]
    # input_PSM = sim_input["input_PSM"]

    # # GENERAL INITIALIsATION
    # ## case settings
    # sim_name = config.sim_name
    # is_with_vk_optimization = config.is_with_vk_optimization
    # is_circular_case = config.is_circular_case
    # is_run_only_1_time_step = config.is_run_only_1_time_step
    # is_print_intermediate_results = config.is_print_intermediate_results
    # is_with_gravity = config.is_with_gravity
    # is_with_velocity_initialisation = config.is_with_velocity_initialisation

    # # logging
    # logging.info(
    #     f"Running aero-structural simulation for {sim_name}, shape: {np.shape(sim_name)}"
    # )
    # logging.debug(
    #     f" type of input_PSM: {type(input_PSM)}, shape: {np.shape(input_PSM)}"
    # )
    # logging.debug(
    #     f" type of input_VSM: {type(input_VSM)}, shape: {np.shape(input_VSM)}"
    # )
    # logging.info(
    #     f" type of input_bridle_aero: {type(input_bridle_aero)}, shape: {np.shape(input_bridle_aero)}"
    # )
    # logging.info(
    #     f" type of input_tether_aero: {type(input_tether_aero)}, shape: {np.shape(input_tether_aero)}"
    # )
    # logging.info(f" type of points: {type(points)}, shape: {np.shape(points)}")
    # logging.info(f" type of vel_app: {type(vel_app)}, shape: {np.shape(vel_app)}")
    # logging.info(f" type of config: {type(config)}, shape: {np.shape(config)}")
    # logging.info(f" type of sim_input: {type(sim_input)}, shape: {np.shape(sim_input)}")
    # logging.info(f" points: {points}")

    ### AERO initialisation -- VSM
    body_aero, vsm_solver = initialize_vsm(
        geometry_csv_path=Path(PROJECT_DIR)
        / "data"
        / "V3_25"
        / "wing_geometry_from_CAD.csv",
        polar_data_dir=Path(PROJECT_DIR) / "data" / "V3_25" / "2D_polars_from_CFD",
        n_panels=9,
        spanwise_panel_distribution="uniform",
        is_half_wing=True,
        is_with_corrected_polar=True,
    )

    ## STRUC initialisation -- PSS
    psystem, params, config = instantiate_psystem(
        PROJECT_DIR, config["kite_name"], config["grav_constant"], config, config_kite
    )

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
    points = config.kite.points_ini
    vel_app = np.array(config.vel_wind - config.vel_kite)
    start_time = time.time()
    is_convergence = False
    residual_f_list = []
    f_tether_drag = np.zeros(3)
    is_residual_below_tol = False
    is_run_only_1_time_step = True

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
    if config.is_with_gravity:
        force_gravity = np.array(
            [np.array(config.grav_constant) * m_pt for m_pt in mass_points]
        )
    else:
        force_gravity = np.zeros(points.shape)

    ## coupling initialisation
    (
        points_wing_segment_corners_aero_orderded,
        index_transformation_struc_to_aero,
    ) = coupling_struc2aero.extract_wingpanel_corners_aero_orderded(
        points, config.kite.connectivity.plate_point_indices
    )
    logging.info(
        f"points_wing_segment_corners_aero_orderded: {points_wing_segment_corners_aero_orderded}"
    )
    n_ribs = int(len(points_wing_segment_corners_aero_orderded) / 2)
    n_panels = int(n_ribs - 1)
    logging.info(f"n_ribs: {n_ribs}, n_panels: {n_panels}")

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

        ## external force
        begin_time_f_ext = time.time()
        ### STRUC --> AERO
        # Ordering the points from left to right (-y to +y), the desired aero-ordering
        points_wing_segment_corners_aero_orderded = points[
            index_transformation_struc_to_aero
        ]
        ### AERO
        force_aero_wing_VSM, body_aero = run_vsm_package(
            body_aero=body_aero,
            solver=vsm_solver,
            le_arr=points_wing_segment_corners_aero_orderded[0::2, :],
            te_arr=points_wing_segment_corners_aero_orderded[1::2, :],
            va_vector=vel_app,
            aero_input_type="reuse_initial_polar_data",
        )
        # reshuffle forces, instead of left-to-right, we make them right-to-left
        force_aero_wing_VSM = force_aero_wing_VSM[::-1]
        ### AERO --> STRUC
        if coupling_method == "NN":
            force_aero_wing = coupling_aero2struc.aero2struc_NN_vsm(
                n_chordwise_aero_nodes,
                body_aero,
                force_aero_wing_VSM,
                points_wing_segment_corners_aero_orderded,
                index_transformation_aero_to_struc_dict,
                points,
            )
        else:
            raise ValueError("Coupling method not recognized; wrong name or typo")

        # TODO: get bridle line forces back in to play
        force_aero_bridle = np.zeros(points.shape)
        force_aero = force_aero_wing + force_aero_bridle

        ## summing up
        force_external = force_aero + force_gravity

        ### STRUC
        ### f_external is flat, and force_external is 2D, ##TODO: could add this to the name?
        f_external = force_external.flatten()
        end_time_f_ext = time.time()
        begin_time_f_int = time.time()
        psystem = run_pss(psystem, params, f_external)
        position.loc[step], _ = psystem.x_v_current
        end_time_f_int = time.time()

        logging.debug(f"position.loc[step].shape: {position.loc[step].shape}")
        logging.debug(f"internal force: {psystem.f_int}")
        logging.debug(f"external force: {f_external}")

        # # TODO: ideally you don't need this here and have aero-also work with the flat format, as this should be faster
        # # saving points in different format
        # points = psystem.x_current_2D
        ## TODO: replacing this function inside the src code to inhere
        points = np.array([particle.x for particle in psystem.particles])
        residual_f = psystem.f_int + f_external
        residual_f_list.append(np.linalg.norm(np.abs(residual_f)))
        print(
            f"i:{i}, t-step:{step:.2f}, residual_f: {np.linalg.norm(residual_f):.3f}N (aero: {end_time_f_ext-begin_time_f_ext:.3f}s, struc: {end_time_f_int-begin_time_f_int:.3f}s)"
        )

        ## calculating delta tape lengths
        delta_depower_tape = (
            psystem.extract_rest_length[index_depower_tape]
            - initial_length_depower_tape
        )

        ### All the convergence checks, are be done in if-elif because only 1 should hold at once
        # if convergence (residual below set tolerance)
        if np.linalg.norm(residual_f) <= aero_structural_tol:
            is_residual_below_tol = True
        if np.linalg.norm(residual_f) <= aero_structural_tol and i > 50:
            is_residual_below_tol = True
            is_convergence = True

        # if residual forces are NaN
        elif np.isnan(np.linalg.norm(residual_f)):
            is_convergence = False
            print("Classic PS diverged - residual force is NaN")
            break
        # if residual forces are not changing anymore
        elif (
            i > 200
            and np.abs(np.mean(residual_f_list[i - 25]) - residual_f_list[i]) < 1
            and np.abs(np.mean(residual_f_list[i - 10]) - residual_f_list[i]) < 1
            and np.abs(np.mean(residual_f_list[i - 5]) - residual_f_list[i]) < 1
            and np.abs(np.mean(residual_f_list[i - 2]) - residual_f_list[i]) < 1
        ):
            is_convergence = False
            print("Classic PS non-converging - residual no longer changes")
            break
        # if to many iterations are needed
        elif i > aero_structural_max_iter:
            is_convergence = False
            print(
                f"Classic PS non-converging - more than max ({aero_structural_max_iter}) iterations needed"
            )
            break
        # special case for running the simulation for only one timestep
        elif is_run_only_1_time_step:
            break

        ## if convergenced and not yet at desired depower-tape length
        elif (
            is_residual_below_tol and delta_depower_tape <= depower_tape_final_extension
        ):
            psystem.update_rest_length(index_depower_tape, depower_tape_extension_step)
            delta_depower_tape = (
                psystem.extract_rest_length[index_depower_tape]
                - initial_length_depower_tape
            )
            print(
                f"||--- delta l_d: {delta_depower_tape:.3f}m | new l_d: {psystem.extract_rest_length[index_depower_tape]:.3f}m"
            )
            is_residual_below_tol = False

        ## if convergenced and all changes are made
        elif is_residual_below_tol:
            is_convergence = True

        if is_convergence:
            break
    ######################################################################
    ## END OF SIMULATION FOR LOOP
    ######################################################################

    # defining post_processing_output
    aero_structural_total_time = time.time() - start_time
    wing_rest_lengths = psystem.extract_rest_length[0:len_wing_rest_length]
    bridle_rest_lengths = psystem.extract_rest_length[len_wing_rest_length:]
    position_without_na = position.dropna(
        how="all",
        subset=position.columns[position.columns.str.contains("[xyz]")],
        axis=0,
    )
    num_of_rows = position_without_na.shape[0]

    sim_output = {
        "points": points,
        "position": position,
        "aero_structural_total_time": aero_structural_total_time,
        "num_of_iterations": i,
        "is_convergence": is_convergence,
        "wing_rest_lengths": wing_rest_lengths,
        "bridle_rest_lengths": bridle_rest_lengths,
        # print data
        "is_convergence": is_convergence,
        "num_of_iterations": i,
        "vel_app": vel_app,
        ## aero_structural_total_time
        "residual_f_including_fixed_nodes": psystem.f,
        "residual_f": residual_f,
        "f_internal": psystem.f_int,
        "f_external": f_external,
        "force_aero": force_aero,
        "force_aero_wing": force_aero_wing,
        "force_aero_bridle": force_aero_bridle,
        "f_tether_drag": f_tether_drag,
        "force_gravity": force_gravity,
        ## wing_rest_lengths
        ## bridle_rest_lengths
        # plot data
        "wingpanels": np.zeros(points.shape),  # wingpanels,
        "controlpoints": np.zeros(points.shape),  # controlpoints,
        "rings": np.zeros(points.shape),  # rings,
        "coord_L": np.zeros(points.shape),  # coord_L,
        "F_rel": np.zeros(points.shape),  # F_rel,
        ## wing_rest_lengths
        ## bridle_rest_lengths
        # animation data
        "position_without_na": position_without_na,
        "num_of_rows": num_of_rows,
        ## wing_rest_lengths
        ## bridle_rest_lengths
    }

    from VSM.interactive import interactive_plot

    points_wing_segment_corners_aero_orderded = points[
        index_transformation_struc_to_aero
    ]
    force_aero_wing_VSM, body_aero = run_vsm_package(
        body_aero=body_aero,
        solver=vsm_solver,
        le_arr=points_wing_segment_corners_aero_orderded[0::2, :],
        te_arr=points_wing_segment_corners_aero_orderded[1::2, :],
        va_vector=vel_app,
        aero_input_type="reuse_initial_polar_data",
    )

    interactive_plot(
        body_aero,
        vel=np.linalg.norm(vel_app),
        angle_of_attack=10,
        side_slip=0,
        yaw_rate=0,
        is_with_aerodynamic_details=True,
        title="TUDELFT_V3_KITE",
    )
    return sim_output


# %%
