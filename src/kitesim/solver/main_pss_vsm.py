import time
from tqdm import tqdm
import numpy as np
import logging
from pathlib import Path
import copy
from kitesim.solver import (
    aerodynamic,
    initialisation,
    struc2aero,
    aero2struc,
    structural,
    tracking,
    plotting,
)


def run_aerostructural_solver(
    config_dict: dict, config_kite_dict: dict, PROJECT_DIR: Path, results_dir: Path
):
    """
    Runs the aero-structural solver for the given input parameters.

    Args:
        config_dict (dict): Main configuration dictionary.
        config_kite_dict (dict): Kite-specific configuration dictionary.
        PROJECT_DIR (Path): Path to the project directory.
        results_dir (Path): Path to the results directory.

    Returns:
        tracking_data (dict): Dictionary containing time histories of positions, forces, etc.
        meta (dict): Dictionary with meta information about the simulation (timing, convergence, etc).
    """

    ## INIT
    (
        struc_nodes,
        wing_ci,
        wing_cj,
        bridle_ci,
        bridle_cj,
        struc_le_idx_list,
        struc_te_idx_list,
        pulley_point_indices,
        tubular_frame_line_idx_list,
        te_line_idx_list,
        n_struc_ribs,
        wing_connectivity,
        bridle_connectivity,
        kite_connectivity,
        m_array,
        bridle_rest_lengths_initial,
        wing_rest_lengths_initial,
        rest_lengths,
        pulley_data,
    ) = initialisation.main(config_kite_dict)

    ## AERO
    n_aero_panels_per_struc_section = config_dict["aero"][
        "n_aero_panels_per_struc_section"
    ]
    n_struc_sections = n_struc_ribs - 1
    n_aero_panels = n_struc_sections * n_aero_panels_per_struc_section
    body_aero, vsm_solver = aerodynamic.initialize_vsm(
        config_kite=config_kite_dict,
        n_panels=n_aero_panels,
        spanwise_panel_distribution="uniform",
    )
    vel_app = np.array(config_dict["vel_wind"]) - np.array(config_dict["vel_kite"])
    body_aero.va = (vel_app, 0)
    wing = body_aero.wings[0]
    new_sections = wing.refine_aerodynamic_mesh()
    initial_polar_data = []
    for new_section in new_sections:
        initial_polar_data.append(new_section.aero_input)

    ## AERO2STRUC
    aero2struc_mapping = aero2struc.initialize_mapping(
        body_aero.panels,
        struc_nodes,
        struc_le_idx_list,
        struc_te_idx_list,
    )

    ## STRUC
    psystem, params, pss_kite_connectivity = structural.instantiate_psystem(
        config_dict,
        config_kite_dict,
        struc_nodes,
        wing_connectivity,
        kite_connectivity,
        rest_lengths,
        m_array,
        tubular_frame_line_idx_list,
        te_line_idx_list,
        pulley_point_indices,
        pulley_data,
    )

    ## ACTUATION
    len_wing_rest_length = len(wing_rest_lengths_initial)
    index_depower_tape = (
        len_wing_rest_length + config_kite_dict["bridle"]["depower_tape_index"]
    )
    initial_length_depower_tape = params["l0"][index_depower_tape]
    depower_tape_extension_step = config_dict["depower_tape_extension_step"]
    depower_tape_final_extension = config_dict["depower_tape_final_extension"]
    n_depower_tape_steps = int(
        depower_tape_final_extension / depower_tape_extension_step
    )
    logging.info(
        f"Initial depower tape length: {psystem.extract_rest_length[index_depower_tape]:.3f}m"
    )
    logging.info(
        f"Desired depower tape length: {initial_length_depower_tape + depower_tape_final_extension:.3f}m"
    )

    ## PRELOOP
    if config_dict["is_with_gravity"]:
        f_ext_gravity = np.array(
            [np.array(config_dict["grav_constant"]) * m_pt for m_pt in m_array]
        )
    else:
        f_ext_gravity = np.zeros(struc_nodes.shape)
    initial_particles = copy.deepcopy(psystem.particles)
    t_vector = np.linspace(
        params["dt"], params["t_steps"] * params["dt"], params["t_steps"]
    )
    tracking_data = tracking.setup_tracking_arrays(len(struc_nodes), t_vector)
    vel_app = np.array(config_dict["vel_wind"]) - np.array(config_dict["vel_kite"])
    is_convergence = False
    f_residual_list = []
    f_tether_drag = np.zeros(3)
    is_residual_below_tol = False
    struc_nodes_prev = None  # Initialize previous points for tracking
    start_time = time.time()

    ######################################################################
    # SIMULATION LOOP
    ######################################################################
    ## propagating the simulation for each timestep and saving results
    with tqdm(total=len(t_vector), desc="Simulating", leave=True) as pbar:
        for i, step in enumerate(t_vector):
            if i > 0:
                struc_nodes_prev = struc_nodes.copy()

            ## external force
            begin_time_f_ext = time.time()

            ### STRUC --> AERO
            le_arr, te_arr = struc2aero.main(
                struc_nodes,
                struc_le_idx_list,
                struc_te_idx_list,
                n_aero_panels_per_struc_section,
            )

            ### AERO
            f_aero_wing_vsm_format, body_aero, results_aero = (
                aerodynamic.run_vsm_package(
                    body_aero=body_aero,
                    solver=vsm_solver,
                    le_arr=le_arr,
                    te_arr=te_arr,
                    va_vector=vel_app,
                    aero_input_type="reuse_initial_polar_data",
                    initial_polar_data=initial_polar_data,
                    is_with_plot=config_dict["is_with_aero_plot_per_iteration"],
                )
            )
            logging.debug(
                f"Aero symmetry check, f_aero_y: { np.sum([force[1] for force in f_aero_wing_vsm_format])}"
            )

            ### AERO --> STRUC
            f_aero_wing = aero2struc.main(
                config_dict["coupling_method"],
                f_aero_wing_vsm_format,
                struc_nodes,
                np.array(results_aero["panel_cp_locations"]),
                aero2struc_mapping,
                p=2,
                eps=1e-6,
                is_with_coupling_plot=config_dict["is_with_coupling_plot"],
            )

            # TODO: get bridle line forces back in to play
            f_aero_bridle = np.zeros(struc_nodes.shape)
            f_aero = f_aero_wing + f_aero_bridle

            ## summing up
            f_ext = f_aero + f_ext_gravity
            f_ext = np.round(f_ext, 5)

            # Checking symmetry in the forces
            print(
                f"np.sum(f_aero): {np.sum(f_aero[:, 0])}, {np.sum(f_aero[:, 1])}, {np.sum(f_aero[:, 2])}"
            )
            print(
                f"np.sum(f_ext): {np.sum(f_ext[:, 0])}, {np.sum(f_ext[:, 1])}, {np.sum(f_ext[:, 2])}"
            )
            if config_dict["is_with_plot_per_iteration"]:
                plotting.main(
                    struc_nodes, pss_kite_connectivity, f_ext=f_ext, title=f"i: {i}"
                )

            ### STRUC
            f_ext_flat = f_ext.flatten()
            end_time_f_ext = time.time()
            begin_time_f_int = time.time()
            psystem = structural.run_pss(psystem, params, f_ext_flat)
            # position.loc[step], _ = psystem.x_v_current
            end_time_f_int = time.time()

            # logging.debug(f"position.loc[step].shape: {position.loc[step].shape}")
            logging.debug(f"internal force: {psystem.f_int}")
            logging.debug(f"external force: {f_ext}")

            # # TODO: ideally you don't need this here and have aero-also work with the flat format, as this should be faster
            # # saving points in different format
            # points = psystem.x_current_2D
            ## TODO: replacing this function inside the src code to inhere
            # Updating the points
            struc_nodes = np.array([particle.x for particle in psystem.particles])
            f_residual = psystem.f_int + f_ext_flat
            f_residual_list.append(np.linalg.norm(np.abs(f_residual)))

            # Update unified tracking dataframe (replaces position update)
            tracking.update_tracking_arrays(
                tracking_data,
                i,
                psystem,
                struc_nodes,
                struc_nodes_prev,
                f_ext_flat,
                f_residual,
            )
            ## calculating delta tape lengths
            delta_depower_tape = (
                psystem.extract_rest_length[index_depower_tape]
                - initial_length_depower_tape
            )

            pbar.set_postfix(
                {
                    "res": f"{np.linalg.norm(f_residual):.3f}N",
                    # "aero": f"{end_time_f_ext-begin_time_f_ext:.2f}s",
                    # "struc": f"{end_time_f_int-begin_time_f_int:.2f}s",
                }
            )
            pbar.update(1)

            ### All the convergence checks, are be done in if-elif because only 1 should hold at once
            # if convergence (residual below set tolerance)
            if (
                i > n_depower_tape_steps
                and np.linalg.norm(f_residual) <= params["aerostructural_tol"]
            ):
                is_residual_below_tol = True
                is_convergence = True
            # if np.linalg.norm(f_residual) <= params["aerostructural_tol"]l and i > 50:
            #     is_residual_below_tol = True
            #     is_convergence = True

            # if residual forces are NaN
            elif np.isnan(np.linalg.norm(f_residual)):
                is_convergence = False
                logging.info("Classic PS diverged - residual force is NaN")
                break
            # if residual forces are not changing anymore
            elif (
                i > 200
                and np.abs(np.mean(f_residual_list[i - 25]) - f_residual_list[i]) < 1
                and np.abs(np.mean(f_residual_list[i - 10]) - f_residual_list[i]) < 1
                and np.abs(np.mean(f_residual_list[i - 5]) - f_residual_list[i]) < 1
                and np.abs(np.mean(f_residual_list[i - 2]) - f_residual_list[i]) < 1
            ):
                is_convergence = False
                logging.info("Classic PS non-converging - residual no longer changes")
                break
            # if to many iterations are needed
            elif i > config_dict["aero_structural"]["max_iter"]:
                is_convergence = False
                logging.info(
                    f"Classic PS non-converging - more than max ({config_dict["aero_structural"]["max_iter"]}) iterations needed"
                )
                break
            # special case for running the simulation for only one timestep
            elif config_dict["is_run_only_1_time_step"]:
                break

            ## if convergenced and not yet at desired depower-tape length
            elif (
                is_residual_below_tol
                and delta_depower_tape <= depower_tape_final_extension
            ):
                psystem.update_rest_length(
                    index_depower_tape, depower_tape_extension_step
                )
                delta_depower_tape = (
                    psystem.extract_rest_length[index_depower_tape]
                    - initial_length_depower_tape
                )
                logging.info(
                    f"||--- delta l_d: {delta_depower_tape:.3f}m | new l_d: {psystem.extract_rest_length[index_depower_tape]:.3f}m | Steps required: {n_depower_tape_steps}"
                )
                is_residual_below_tol = False

            ## if convergenced and all changes are made
            elif is_residual_below_tol:
                is_convergence = True

            if is_convergence:
                n_iter = i + 1
                break
    ######################################################################
    ## END OF SIMULATION FOR LOOP
    ######################################################################
    if config_dict["is_with_final_plot"]:
        plotting.main(
            np.array([particle.x for particle in psystem.particles]),
            pss_kite_connectivity,
            struc_nodes_initial=np.array([p.x for p in initial_particles]),
            title="PSM: Initial (black) vs Final (red)",
        )
    meta = {
        "total_time_s": time.time() - start_time,
        "n_iter": n_iter,
        "converged": is_convergence,
        "wing_rest_lengths": psystem.extract_rest_length[
            :len_wing_rest_length
        ].tolist(),
        "bridle_rest_lengths": psystem.extract_rest_length[
            len_wing_rest_length:
        ].tolist(),
    }
    return tracking_data, meta
