import time
from tqdm import tqdm
import numpy as np
import logging
from pathlib import Path
import copy
from kitesim import (
    aerodynamic_vsm,
    struc2aero,
    aero2struc,
    structural_pss,
    tracking,
    plotting,
    structural_kite_fem,
    aerodynamic_bridle_line_drag,
)


# Remove hardcoded values, when changing away from V3
def forcing_symmetry(struc_nodes):
    """
    Forcing symmetry in the y-direction for the kite structure nodes.
    This is a temporary solution to ensure symmetry in the simulation.
    """
    symmetry_pairs_dict = {
        1: 19,
        2: 20,
        3: 17,
        4: 18,
        5: 15,
        6: 16,
        7: 13,
        8: 14,
        9: 11,
        10: 12,
        # bridles
        21: 24,
        22: 23,
        25: 26,
        27: 30,
        28: 29,
        31: 32,
        33: 35,
        36: 37,
    }

    for key, value in symmetry_pairs_dict.items():
        struc_nodes[value] = np.array(
            [struc_nodes[key][0], -struc_nodes[key][1], struc_nodes[key][2]]
        )
    struc_nodes[34][1] = 0
    return struc_nodes


def update_power_tape_actuation(
    config,
    psystem,
    kite_fem_structure,
    power_tape_index,
    power_tape_extension_step,
    initial_length_power_tape,
    power_tape_final_extension,
    is_residual_below_tol,
    n_power_tape_steps,
    rest_lengths=None,
):
    """
    Calculate current power tape extension and update if needed for actuation.

    Args:
        config: Configuration dictionary
        psystem: Particle system (for PSS solver)
        kite_fem_structure: FEM structure (for kite_fem solver)
        power_tape_index: Index of power tape in connectivity array
        power_tape_extension_step: Increment for power tape extension
        initial_length_power_tape: Initial length of power tape
        power_tape_final_extension: Final desired power tape extension
        is_residual_below_tol: Flag indicating if residual is below tolerance
        n_power_tape_steps: Number of power tape extension steps
        rest_lengths: Current rest lengths array (for kite_fem solver)

    Returns:
        tuple: (delta_power_tape, is_actuation_finalized)
            - delta_power_tape: Current change in power tape length
            - is_actuation_finalized: True if actuation is complete, False otherwise
    """
    is_actuation_finalized = True

    ## Calculate delta tape lengths based on structural solver
    if config["structural_solver"] == "pss":
        delta_power_tape = (
            psystem.extract_rest_length[power_tape_index] - initial_length_power_tape
        )

        # Check if we need to extend the power tape
        if is_residual_below_tol and delta_power_tape <= power_tape_final_extension:
            psystem.update_rest_length(power_tape_index, power_tape_extension_step)
            delta_power_tape = (
                psystem.extract_rest_length[power_tape_index]
                - initial_length_power_tape
            )
            logging.info(
                f"||--- delta l_d: {delta_power_tape:.3f}m | new l_d: {psystem.extract_rest_length[power_tape_index]:.3f}m | Steps required: {n_power_tape_steps}"
            )
            is_actuation_finalized = False

    elif config["structural_solver"] == "kite_fem":
        # TODO: get power tape actuation working again for kite_fem
        delta_power_tape = rest_lengths[power_tape_index] + power_tape_extension_step
        rest_lengths = kite_fem_structure.modify_get_spring_rest_length(
            spring_ids=[power_tape_index], new_l0s=[delta_power_tape]
        )
        delta_power_tape = 0
        logging.warning("NOT IMPLEMENTED: delta_power_tape for kite_fem")

    return delta_power_tape, is_actuation_finalized


# TODO: this should also use structural is not converging
def check_convergence(
    i,
    f_residual,
    f_residual_list,
    f_aero_wing_vsm_format,
    config,
):
    """
    Check convergence conditions for the aero-structural solver.

    Args:
        i: Current iteration number
        f_residual: Current residual force vector
        f_residual_list: List of residual force norms from all iterations
        f_aero_wing_vsm_format: Aerodynamic forces in VSM format
        config: Configuration dictionary

    Returns:
        tuple: (is_convergence, should_break)
            - is_convergence: True if converged, False otherwise
            - should_break: True if loop should break, False to continue
    """
    is_convergence = False
    should_break = False

    ### All the convergence checks, are be done in if-elif because only 1 should hold at once
    # if convergence (residual below set tolerance)
    if np.linalg.norm(f_residual) <= config["aero_structural_solver"]["tol"]:
        is_convergence = True

    # if residual forces are NaN
    elif np.isnan(np.linalg.norm(f_residual)):
        is_convergence = False
        logging.info("Classic PS diverged - residual force is NaN")
        should_break = True

    # if residual forces are not changing anymore
    elif (
        i > 15
        and np.abs(np.mean(f_residual_list[i - 15]) - f_residual_list[i - 10]) < 1
    ):
        is_convergence = False
        logging.info("Classic PS non-converging - residual no longer changes")
        should_break = True

    # if too many iterations are needed
    elif i > config["aero_structural_solver"]["max_iter"]:
        is_convergence = False
        logging.info(
            f"Classic PS non-converging - more than max ({config['aero_structural_solver']['max_iter']}) iterations needed"
        )
        should_break = True

    # special case for running the simulation for only one timestep
    elif config["is_run_only_1_time_step"]:
        should_break = True

    # when aero does not converge
    elif np.sum([force[1] for force in f_aero_wing_vsm_format]) == np.nan:
        is_convergence = False
        logging.info("Classic PS non-converging - aero forces are NaN")
        should_break = True

    return is_convergence, should_break


def main(
    m_arr=None,
    struc_nodes=None,
    struc_nodes_initial=None,
    config=None,
    ### ACTUATION
    initial_length_power_tape=None,
    n_power_tape_steps=None,
    power_tape_final_extension=None,
    power_tape_extension_step=None,
    ### CONNECTIVITY
    kite_connectivity_arr=None,
    bridle_connectivity_arr=None,
    pulley_line_indices=None,
    pulley_line_to_other_node_pair_dict=None,
    ### STRUC --> AERO
    struc_node_le_indices=None,
    struc_node_te_indices=None,
    ### AERO
    body_aero=None,
    vsm_solver=None,
    vel_app=None,
    initial_polar_data=None,
    bridle_diameter_arr=None,
    ### AERO --> STRUC
    aero2struc_mapping=None,
    power_tape_index=None,
    ### STRUC
    psystem=None,
    kite_fem_structure=None,
):
    """
    Runs the aero-structural solver for the given input parameters.

    Args:
        config (dict): Main configuration dictionary.
        PROJECT_DIR (Path): Path to the project directory.
        results_dir (Path): Path to the results directory.

    Returns:
        tracking_data (dict): Dictionary containing time histories of positions, forces, etc.
        meta (dict): Dictionary with meta information about the simulation (timing, convergence, etc).
    """

    ## PRELOOP
    if config["is_with_gravity"]:
        f_ext_gravity = np.array(
            [np.array(config["grav_constant"]) * m_pt for m_pt in m_arr]
        )
    else:
        f_ext_gravity = np.zeros(struc_nodes.shape)

    if config["structural_solver"] == "kite_fem":
        rest_lengths = kite_fem_structure.modify_get_spring_rest_length()

    t_vector = np.linspace(
        1,
        config["aero_structural_solver"]["max_iter"],
        config["aero_structural_solver"]["max_iter"],
    )
    tracking_data = tracking.setup_tracking_arrays(len(struc_nodes), t_vector)
    is_convergence = False
    f_residual_list = []
    f_tether_drag = np.zeros(3)
    is_residual_below_tol = False
    struc_nodes_prev = None  # Initialize previous points for tracking
    start_time = time.time()
    plotting.set_plot_style()

    ## track initial state
    # Update unified tracking dataframe (replaces position update)
    tracking.update_tracking_arrays(
        tracking_data,
        0,
        struc_nodes,
        np.zeros(np.shape(struc_nodes.flatten())),
        np.zeros(np.shape(struc_nodes.flatten())),
    )

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
                struc_node_le_indices,
                struc_node_te_indices,
                config["aerodynamic"]["n_aero_panels_per_struc_section"],
            )

            ### AERO
            f_aero_wing_vsm_format, body_aero, results_aero = (
                aerodynamic_vsm.run_vsm_package(
                    body_aero=body_aero,
                    solver=vsm_solver,
                    le_arr=le_arr,
                    te_arr=te_arr,
                    va_vector=vel_app,
                    aero_input_type="reuse_initial_polar_data",
                    initial_polar_data=initial_polar_data,
                    is_with_plot=config["is_with_aero_plot_per_iteration"],
                )
            )
            logging.debug(
                f"Aero symmetry check, f_aero_y: {np.sum([force[1] for force in f_aero_wing_vsm_format])}"
            )
            ### AERO --> STRUC
            f_aero_wing = aero2struc.main(
                config["aero2struc"]["coupling_method"],
                f_aero_wing_vsm_format,
                struc_nodes,
                np.array(results_aero["panel_cp_locations"]),
                aero2struc_mapping,
                config["is_with_coupling_plot_per_iteration"],
                config["aero2struc"],
            )

            ### BRIDLE AERO
            f_aero_bridle = aerodynamic_bridle_line_drag.main(
                struc_nodes,
                bridle_connectivity_arr,
                bridle_diameter_arr,
                vel_app,
                config["rho"],
                config["aerodynamic_bridle"]["cd_cable"],
                config["aerodynamic_bridle"]["cf_cable"],
            )
            f_aero = f_aero_wing + f_aero_bridle

            ## EXTERNAL FORCE
            f_ext = f_aero + f_ext_gravity
            f_ext = np.round(f_ext, 5)
            f_ext_flat = f_ext.flatten()

            if config["is_with_struc_plot_per_iteration"]:
                if config["structural_solver"] == "pss":
                    rest_lengths = psystem.extract_rest_length
                elif config["structural_solver"] == "kite_fem":
                    rest_lengths = kite_fem_structure.modify_get_spring_rest_length()
                    kite_fem_structure.plot_convergence()

                plotting.main(
                    struc_nodes,
                    kite_connectivity_arr,
                    rest_lengths,
                    f_ext=f_ext,
                    title=f"i: {i}",
                    body_aero=body_aero,
                    is_with_node_indices=False,
                    pulley_line_indices=pulley_line_indices,
                    pulley_line_to_other_node_pair_dict=pulley_line_to_other_node_pair_dict,
                )
                
            ### RESIDUAL FEM
            if config["structural_solver"] == "kite_fem":
                if i > 0:
                    f_residual = f_int + f_ext_flat
                if i ==0:
                    f_residual = f_ext_flat
                f_residual_list.append(np.linalg.norm(np.abs(f_residual)))
                logging.debug(
                    f"residual force in y-direction: {np.sum([f_residual[1::3]]):.3f}N"
                )
            
            ### STRUC
            end_time_f_ext = time.time()
            begin_time_f_int = time.time()
            if config["structural_solver"] == "pss":
                psystem, is_structural_converged, struc_nodes, f_int = (
                    structural_pss.run_pss(
                        psystem,
                        f_ext_flat,
                        config["structural_pss"],
                    )
                )
            elif config["structural_solver"] == "kite_fem":
                kite_fem_structure, is_structural_converged, struc_nodes, f_int = (
                    structural_kite_fem.run_kite_fem(
                        kite_fem_structure, f_ext_flat, config["structural_kite_fem"]
                    )
                )
            end_time_f_int = time.time()

            ### FORCING SYMMETRY
            if config["is_with_forcing_symmetry"]:
                logging.info("Forcing symmetry in y-direction")
                struc_nodes = forcing_symmetry(struc_nodes)

            ### RESIDUAL PSM
            if config["structural_solver"] == "pss":
                f_residual = f_int + f_ext_flat
                f_residual_list.append(np.linalg.norm(np.abs(f_residual)))
                logging.debug(
                    f"residual force in y-direction: {np.sum([f_residual[1::3]]):.3f}N"
                )

            ### TRACKING
            # Update unified tracking dataframe (replaces position update)
            tracking.update_tracking_arrays(
                tracking_data,
                i,
                struc_nodes,
                f_ext_flat,
                f_residual,
            )

            ### PROGRESS BAR
            pbar.set_postfix(
                {
                    "res": f"{np.linalg.norm(f_residual):.3f}N",
                    "aero": f"{end_time_f_ext-begin_time_f_ext:.2f}s",
                    "struc": f"{end_time_f_int-begin_time_f_int:.2f}s",
                }
            )
            pbar.update(1)

            ### CHECK CONVERGENCE
            is_convergence, should_break = check_convergence(
                i=i,
                f_residual=f_residual,
                f_residual_list=f_residual_list,
                f_aero_wing_vsm_format=f_aero_wing_vsm_format,
                config=config,
            )

            ### ACTUATION (only when converged)
            if is_convergence:
                # Update residual flag for actuation function
                is_residual_below_tol = is_convergence

                delta_power_tape, is_actuation_finalized = update_power_tape_actuation(
                    config=config,
                    psystem=psystem,
                    kite_fem_structure=kite_fem_structure,
                    power_tape_index=power_tape_index,
                    power_tape_extension_step=power_tape_extension_step,
                    initial_length_power_tape=initial_length_power_tape,
                    power_tape_final_extension=power_tape_final_extension,
                    is_residual_below_tol=is_residual_below_tol,
                    n_power_tape_steps=n_power_tape_steps,
                    rest_lengths=(
                        rest_lengths
                        if config["structural_solver"] == "kite_fem"
                        else None
                    ),
                )

                # If actuation not finalized, continue to next iteration
                if not is_actuation_finalized:
                    # ACTUATION PHASE: Continue until power tape reaches final extension
                    continue

            # Check if we should exit the loop
            if should_break or (is_convergence and is_actuation_finalized):
                break
    ######################################################################
    ## END OF SIMULATION FOR LOOP
    ######################################################################
    if config["structural_solver"] == "pss":
        rest_lengths = psystem.extract_rest_length
    elif config["structural_solver"] == "kite_fem":
        rest_lengths = kite_fem_structure.modify_get_spring_rest_length()

    if config["is_with_final_plot"]:
        plotting.main(
            struc_nodes,
            kite_connectivity_arr,
            f_ext=f_ext,
            rest_lengths=rest_lengths,
            struc_nodes_initial=struc_nodes_initial,
            title="Initial vs final",
            pulley_line_indices=pulley_line_indices,
            pulley_line_to_other_node_pair_dict=pulley_line_to_other_node_pair_dict,
        )
    meta = {
        "total_time_s": time.time() - start_time,
        "n_iter": i + 1,
        "converged": is_convergence,
        "rest_lengths": rest_lengths,  # ensure numeric array
        # Convert kite_connectivity to a numeric array for HDF5 compatibility
        "kite_connectivity": np.array(
            [[int(row[0]), int(row[1])] for row in np.array(kite_connectivity_arr)],
            dtype=np.int32,
        ),
    }

    return tracking_data, meta
