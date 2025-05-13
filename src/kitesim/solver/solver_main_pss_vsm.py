# %% Defing the spring-damper system function
# defining function for the spring-damper system

import time
import numpy as np
import pandas as pd
import logging
import os
import sys
from pathlib import Path
from kitesim.VortexStepMethod import VSM_old
from kitesim.aerodynamic import bridle_line_system_aero
from kitesim.coupling import coupling_struc2aero, coupling_aero2struc
from kitesim.solver import solver_utils
from kitesim.aerodynamic import tether_aero, bridle_line_system_aero
from kitesim.particleSystem.ParticleSystem import ParticleSystem
from kitesim.VortexStepMethod import VSM_old
from PSS.particleSystem.SpringDamper import SpringDamperType
import PSS.particleSystem as PSS
import VSM.WingGeometry as VSM_WingGeometry
import VSM.WingAerodynamics as VSM_WingAerodynamics
import VSM.Solver as VSM_Solver
import VSM.plotting as VSM_plotting


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


def prepare_pss_input(input_PSM, grav_constant):
    logging.debug(
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
    pss_initial_conditions = input_PSM["initial_conditions"]
    logging.debug(f"pss_initial_conditions: {pss_initial_conditions}")

    # refining the params_dict
    pss_param_dict = input_PSM["params_dict"]
    pss_param_dict.update({"g": -grav_constant[2]})
    logging.debug(f"pss_param_dict: {pss_param_dict.keys()}")

    # restructuring connectivity matrix
    connectivity_matrix = input_PSM["connectivity_matrix"]

    pss_connectivity_matrix = []
    for idx, _ in enumerate(connectivity_matrix):
        if pss_param_dict["is_compression"][idx] and pss_param_dict["is_tension"][idx]:
            linktype = "default"
        elif (
            pss_param_dict["is_compression"][idx]
            and not pss_param_dict["is_tension"][idx]
        ):
            linktype = "nontensile"
        elif pss_param_dict["is_pulley"][idx]:
            linktype = "pulley"
        elif (
            pss_param_dict["is_tension"][idx]
            and not pss_param_dict["is_compression"][idx]
        ):
            linktype = "noncompressive"

        logging.debug(f"idx: {idx}")
        logging.debug(f"connectivity_matrix[idx]: {connectivity_matrix[idx]}")
        logging.debug(f"pss_param_dict['k'][idx]: {pss_param_dict['k'][idx]}")
        logging.debug(f"pss_param_dict['c']: {pss_param_dict['c']}")
        logging.debug(f"linktype: {linktype}")

        pss_connectivity_matrix.append(
            [
                int(connectivity_matrix[idx][0]),
                int(connectivity_matrix[idx][1]),
                float(pss_param_dict["k"][idx]),
                float(pss_param_dict["c"]),
                string_to_springdampertype(linktype),
            ]
        )

    logging.debug(f"pss_connectivity_matrix: {pss_connectivity_matrix}")

    return pss_connectivity_matrix, pss_initial_conditions, pss_param_dict


def run_pss(psystem, params, f_external):
    f_ext = f_external
    t_vector_internal = np.linspace(
        params["dt"], params["t_steps"] * params["dt"], params["t_steps"]
    )
    E_kin = []
    f_int = []
    E_kin_tol = 1e-3  # 1e-29

    logging.debug(f"Running PS simulation, f_int: {psystem.f_int}")

    # And run the simulation
    for step_internal in t_vector_internal:
        psystem.kin_damp_sim(f_ext)

        E_kin.append(np.linalg.norm(psystem.x_v_current[1] ** 2))
        f_int.append(np.linalg.norm(psystem.f_int))

        converged = False
        if step_internal > 10:
            if np.max(E_kin[-10:-1]) <= E_kin_tol:
                converged = True
        if converged and step_internal > 1:
            print("Kinetic damping PS converged", step_internal)
            break
    return psystem


def grab_V3_25_csv_wing_data():
    # Find the root directory of the repository
    root_dir = os.path.abspath(os.path.dirname(__file__))
    while not os.path.isfile(os.path.join(root_dir, ".gitignore")):
        root_dir = os.path.abspath(os.path.join(root_dir, ".."))
        if root_dir == "/":
            raise FileNotFoundError(
                "Could not find the root directory of the repository."
            )

    ### rib_list_from_CAD_LE_TE_and_surfplan_d_tube_camber_10ribs
    csv_file_path = (
        Path(root_dir)
        / "data"
        / "V3_25"
        / "rib_list_from_CAD_LE_TE_and_surfplan_d_tube_camber_10ribs.csv"
    )
    (
        LE_x_array,
        LE_y_array,
        LE_z_array,
        TE_x_array,
        TE_y_array,
        TE_z_array,
        d_tube_array,
        camber_array,
    ) = np.loadtxt(csv_file_path, delimiter=",", skiprows=1, unpack=True)
    return (
        LE_x_array,
        LE_y_array,
        LE_z_array,
        TE_x_array,
        TE_y_array,
        TE_z_array,
        d_tube_array,
        camber_array,
    )


def run_aerostructural_solver(sim_input):
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

    # logging
    logging.info(
        f"Running aero-structural simulation for {sim_name}, shape: {np.shape(sim_name)}"
    )
    logging.debug(
        f" type of input_PSM: {type(input_PSM)}, shape: {np.shape(input_PSM)}"
    )
    logging.debug(
        f" type of input_VSM: {type(input_VSM)}, shape: {np.shape(input_VSM)}"
    )
    logging.info(
        f" type of input_bridle_aero: {type(input_bridle_aero)}, shape: {np.shape(input_bridle_aero)}"
    )
    logging.info(
        f" type of input_tether_aero: {type(input_tether_aero)}, shape: {np.shape(input_tether_aero)}"
    )
    logging.info(f" type of points: {type(points)}, shape: {np.shape(points)}")
    logging.info(f" type of vel_app: {type(vel_app)}, shape: {np.shape(vel_app)}")
    logging.info(f" type of config: {type(config)}, shape: {np.shape(config)}")
    logging.info(f" type of sim_input: {type(sim_input)}, shape: {np.shape(sim_input)}")
    logging.info(f" points: {points}")

    ## Instantiating the VSM
    n_panels = 18
    spanwise_panel_distribution = "unchanged"
    wing = VSM_WingGeometry.Wing(n_panels, spanwise_panel_distribution)
    vsm_solver = VSM_Solver.Solver(
        aerodynamic_model_type="VSM", is_with_artificial_damping=True
    )
    (
        LE_x_array,
        LE_y_array,
        LE_z_array,
        TE_x_array,
        TE_y_array,
        TE_z_array,
        d_tube_array,
        camber_array,
    ) = grab_V3_25_csv_wing_data()

    ## Instating the PSS
    connectivity_matrx, initial_conditions, params = prepare_pss_input(
        input_PSM, config.grav_constant
    )
    psystem = PSS.ParticleSystem(
        connectivity_matrx,
        initial_conditions,
        params,
    )

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

        ### Struc --> aero
        # Ordering the points from left to right (-y to +y), the desired aero-ordering
        points_wing_segment_corners_aero_orderded = points[
            index_transformation_struc_to_aero
        ]
        # Wing Aerodynamic
        (
            force_aero_wing_VSM_old,
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

        def run_vsm(wing_nodes, wing, vel_app, vsm_solver):
            # flipping the points to the correct order
            wing_nodes_flipped = np.copy(wing_nodes)[::-1]

            # populate the wing object
            for i in range(0, n_ribs):
                index = i * 2
                LE = wing_nodes_flipped[index + 1]
                TE = wing_nodes_flipped[index]
                d_tube = d_tube_array[i]
                camber = camber_array[i]
                print(
                    f"i, {i}, LE: {LE}, TE: {TE}, d_tube: {d_tube:.2f}, camber: {camber:.2f}"
                )
                wing.add_section(LE, TE, ["lei_airfoil_breukels", [d_tube, camber]])

            # instantiate the wing_aero object
            wing_aero = VSM_WingAerodynamics.WingAerodynamics([wing])
            # set va
            wing_aero.va = vel_app
            logging.info(f"wing_aero.va: {wing_aero.va}")

            # plotting the geometry
            VSM_plotting.plot_geometry(
                wing_aero,
                title="wing_aero_CAD_10ribs",
                data_type=".pdf",
                save_path="",
                is_save=False,
                is_show=True,
                view_elevation=15,
                view_azimuth=-120,
            )
            # Solve the wing aerodynamics
            results_aero = vsm_solver.solve(wing_aero)
            force_aero_wing_VSM_flipped = results_aero["F_distribution"]

            # flip it back
            force_aero_wing_VSM = force_aero_wing_VSM_flipped[::-1]
            return force_aero_wing_VSM, wing_aero

        force_aero_wing_VSM, wing_aero = run_vsm(
            points_wing_segment_corners_aero_orderded, wing, vel_app, vsm_solver
        )
        logging.info(
            f"force_aero_wing_VSM: {np.shape(force_aero_wing_VSM)} \n {force_aero_wing_VSM}"
        )
        logging.info(
            f"force_aero_wing_VSM_old: {np.shape(force_aero_wing_VSM_old)} \n {force_aero_wing_VSM_old}"
        )

        # logging.info(
        #     f"force_delta: {np.sum(force_aero_wing_VSM-force_aero_wing_VSM_old)}"
        # )
        # logging.info(f"force_delta: {force_aero_wing_VSM-force_aero_wing_VSM_old}")

        # Aero --> struc
        if coupling_method == "NN":
            # force_aero_wing = coupling_aero2struc.aero2struc_NN(
            #     n_chordwise_aero_nodes,
            #     wingpanels,
            #     force_aero_wing_VSM,
            #     points_wing_segment_corners_aero_orderded,
            #     index_transformation_aero_to_struc_dict,
            #     points,
            # )
            force_aero_wing = coupling_aero2struc.aero2struc_NN_vsm(
                n_chordwise_aero_nodes,
                wing_aero,
                force_aero_wing_VSM,
                points_wing_segment_corners_aero_orderded,
                index_transformation_aero_to_struc_dict,
                points,
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

        ## PSS
        begin_time_f_int = time.time()
        psystem = run_pss(psystem, params, f_external)
        position.loc[step], _ = psystem.x_v_current
        # # TODO: remove
        # # an attempt at relaxation
        # if i > 0:
        #     position.loc[step] = (
        #         0.2 * position.loc[step] + 0.8 * position.loc[t_vector[i - 1]]
        #     )
        end_time_f_int = time.time()

        logging.debug(f"position.loc[step].shape: {position.loc[step].shape}")
        logging.debug(f"internal force: {psystem.f_int}")
        logging.debug(f"external force: {f_external}")

        # # TODO: ideally you don't need this here and have aero-also work with the flat format, as this should be faster
        # # saving points in different format
        # points = psystem.x_current_2D
        ## TODO: replacing this function inside the src code to inhere
        points = np.array([particle.x for particle in psystem.particles])

        # # TODO: this only works if there is one fixed node and it is the first node
        # # remove the fixed nodes residual forces
        # # f_internal = psystem.f_int #this is now without the fixed nodes
        # f_internal = psystem.f_int
        # # TODO: do you really need the copy here? doesn't this blow up the memory?
        # residual_f_including_fixed_nodes = np.copy(f_external + f_internal)
        # residual_f = np.delete(
        #     residual_f_including_fixed_nodes,
        #     slice(int(index_fixed_node * 3), int(index_fixed_node * 3 + 3)),
        # )
        # residual_f_list.append(np.linalg.norm(np.abs(residual_f)))
        residual_f = psystem.f_int + f_external
        residual_f_list.append(np.linalg.norm(np.abs(residual_f)))

        print(
            f"i:{i}, t-step:{step:.2f}, residual_f: {np.linalg.norm(residual_f):.3f}N (aero: {end_time_f_ext-begin_time_f_ext:.3f}s, struc: {end_time_f_int-begin_time_f_int:.3f}s)"
        )
        if is_circular_case and is_print_intermediate_results:
            print(
                f"r_0: {r_0:.2f}m (r_k: {(r_0+center_of_gravity[1]):.2f}m, cg_side: {center_of_gravity[1]:.2f}m)"
            )

            # Should be negative (oriented in negative-side/y direction)
            force_aero_sum = np.sum(force_aero, axis=0)
            f_aero_side = force_aero_sum[1]
            f_tether_side = -force_aero_sum[2] * np.tan(np.arcsin(r_0 / length_tether))
            # Should be positive (oriented in positive-side/y direction)
            f_centrifugal_side = np.sum(mass_points) * (-vel_app[0] ** 2 / r_0)
            print(
                f"Fa_side: {f_aero_side:.3f}N, Ft_side: {f_tether_side:.3f}N, Fc_side: {f_centrifugal_side:.3f}N"
            )

            # # print(
            # #     f"SUM (Fc+Fa+Ft): {f_aero_side+f_tether_side+f_centrifugal_side:.3f}N (should be close to zero)"
            # # )
            # Ma_side = np.sum(np.cross(points, force_aero), axis=0)[0]
            # Mc_side = np.sum(np.cross(points, force_centrifugal), axis=0)[0]
            # print(
            #     f"Mx: {Ma_side+Mc_side:.1f}Nm (Ma_side: {Ma_side:.1f}Nm, Mc_side: {Mc_side:.1f}Nm)"
            # )
            moment_aero_sum = np.sum(np.cross(points, force_aero), axis=0)
            moment_centrifugal_sum = np.sum(np.cross(points, force_centrifugal), axis=0)
            moment_external_sum = np.sum(np.cross(points, force_external), axis=0)
            print(
                f"Mx(roll)  : {moment_external_sum[0]:.1f}Nm (Ma_x: {moment_aero_sum[0]:.1f}Nm, Mc_x: {moment_centrifugal_sum[0]:.1f}Nm)"
            )
            # print(f"My(pitch) : {moment_external_sum[1]:.1f}Nm")
            # print(f"Mz(yaw)   : {moment_external_sum[2]:.1f}Nm")
            print(f" ")

        ## calculating delta tape lengths
        delta_depower_tape = (
            psystem.extract_rest_length[index_depower_tape]
            - initial_length_depower_tape
        )
        delta_steering_tape_right = (
            psystem.extract_rest_length[index_steering_tape_right]
            - initial_length_steering_tape_right
        )

        ##TODO:REMOVE
        ## if convergenced and not yet at desired steering-tape length
        if (
            i > 10
            and is_circular_case
            and np.abs(delta_steering_tape_right)
            <= np.abs(steering_tape_final_extension)
        ):

            # update lenghts, positive extension means left tape shorter, right longer
            ##TODO: changed this, check textual understanding
            psystem.update_rest_length(
                index_steering_tape_left, steering_tape_extension_step
            )
            psystem.update_rest_length(
                index_steering_tape_right, -steering_tape_extension_step
            )
            delta_steering_tape_right = (
                psystem.extract_rest_length[index_steering_tape_right]
                - initial_length_steering_tape_right
            )
            print(
                f"||--- delta l_s: {delta_steering_tape_right:.3f}m | new l_s right: {psystem.extract_rest_length[index_steering_tape_right]:.3f} (& left: {psystem.extract_rest_length[index_steering_tape_left]:.3f}) m"
            )
            is_residual_below_tol = False

        ## Increasing damping by 0.005 after x iterations
        # if i > 200: #and i < 300:
        #     damping_ratio += -0.0001 #-0.00125
        #     psystem.update_damping(damping_ratio)
        # elif i > 2400 and i < 2500:
        #     damping_ratio += -0.005
        #     psystem.update_damping(damping_ratio)

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

        ### Update order - velocity - depower tape - steering tape
        ## if vel not yet at real value, update vel
        if is_with_velocity_initialisation and i < n_vel_initialisation_steps:
            vel_app = vel_app_linspace[i]
            print(f"||--- update vel_app: {vel_app}")

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

        ## if convergenced and not yet at desired steering-tape length
        elif (
            is_residual_below_tol
            and is_circular_case
            and delta_steering_tape_right <= steering_tape_final_extension
        ):
            # update lenghts, positive extension means left tape shorter, right longer
            ##TODO: changed this, check textual understanding
            psystem.update_rest_length(
                index_steering_tape_left, -steering_tape_extension_step
            )
            psystem.update_rest_length(
                index_steering_tape_right, steering_tape_extension_step
            )
            delta_steering_tape_right = (
                psystem.extract_rest_length[index_steering_tape_right]
                - initial_length_steering_tape_right
            )
            print(
                f"||--- delta l_s: {delta_steering_tape_right:.3f}m | new l_s right: {psystem.extract_rest_length[index_steering_tape_right]:.3f} (& left: {psystem.extract_rest_length[index_steering_tape_left]:.3f}) m"
            )
            is_residual_below_tol = False

        ## if convergenced and all changes are made
        elif is_residual_below_tol:
            is_convergence = True
            # if not is_circular_case:
            #     is_convergence = True
            #     break
            # elif is_circular_case and i > 100:
            #     is_convergence = True
            #     break

        if is_convergence:
            break
    ######################################################################
    ## END OF SIMULATION FOR LOOP
    ######################################################################

    print(
        f"delta_steering_tape_right = {(psystem.extract_rest_length[index_steering_tape_right]- initial_length_steering_tape_right)}"
    )

    # ------------------------------------------------------
    # Running Vk_x solver again
    if is_convergence and is_with_vk_optimization:
        print(f"Converged at: i{i}, final update of vel_kite_x")
        # Running Vk_x solver again
        tol_vk_optimization = 1e-1
        vel_app[0] = -solver_utils.optimalisation_of_vk_for_fx_0(
            -vel_app[0],
            vel_app,
            points,
            config,
            input_VSM,
            input_bridle_aero,
            tol_vk_optimization,
        )

    # TODO: this should be one list, not different ones. You are now sending duplicate info

    # defining post_processing_output
    aero_structural_total_time = time.time() - start_time
    wing_rest_lengths = psystem.extract_rest_length[0:len_wing_rest_length]
    bridle_rest_lengths = psystem.extract_rest_length[len_wing_rest_length:]
    position_without_na = position.dropna(
        how="all",
        subset=position.columns[position.columns.str.contains("[xyz]\d+")],
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
        "wingpanels": wingpanels,
        "controlpoints": controlpoints,
        "rings": rings,
        "coord_L": coord_L,
        "F_rel": F_rel,
        ## wing_rest_lengths
        ## bridle_rest_lengths
        # animation data
        "position_without_na": position_without_na,
        "num_of_rows": num_of_rows,
        ## wing_rest_lengths
        ## bridle_rest_lengths
    }
    return sim_output


# %%
