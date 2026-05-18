import time
from tqdm import tqdm
import numpy as np
import logging
from pathlib import Path
import copy
from kitesim import (
    aero2struc_level_1,
    aerodynamic_vsm,
    struc2aero,
    structural_kite_fem_level_1,
    structural_pss,
    tracking,
    plotting,
    aerodynamic_bridle_line_drag,
)
from kitesim.utils import calculate_cg, calculate_inertia, load_yaml, rotate_geometry


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


def _compute_power_tape_increment(
    delta_power_tape,
    power_tape_final_extension,
    power_tape_extension_step,
    tol=1e-9,
):
    """
    Compute the signed rest-length increment needed to move toward the target extension.

    Returns:
        tuple: (increment, should_update)
    """
    remaining = power_tape_final_extension - delta_power_tape
    if np.abs(remaining) <= tol:
        return 0.0, False
    if np.abs(power_tape_extension_step) <= tol:
        return 0.0, False

    # Always move toward target and clamp to avoid overshoot.
    increment = np.sign(remaining) * min(
        np.abs(power_tape_extension_step), np.abs(remaining)
    )
    return increment, True


def _find_kite_fem_spring_id_from_connectivity(
    kite_fem_structure,
    kite_connectivity_arr,
    connectivity_idx,
):
    """
    Map ASKITE connectivity index to the matching kite_fem spring element index.
    """
    ci, cj = [int(v) for v in kite_connectivity_arr[connectivity_idx]]
    target_key = (min(ci, cj), max(ci, cj))

    for spring_id, spring_element in enumerate(kite_fem_structure.spring_elements):
        n1 = int(spring_element.spring.n1)
        n2 = int(spring_element.spring.n2)
        if (min(n1, n2), max(n1, n2)) == target_key:
            return spring_id

    raise ValueError(
        f"Could not map power_tape connectivity index {connectivity_idx} "
        f"with nodes ({ci}, {cj}) to a kite_fem spring element."
    )


def update_steering_tape_actuation(
    config,
    psystem,
    kite_fem_structure,
    kite_connectivity_arr,
    steering_tape_indices,
    steering_tape_extension_step,
    initial_length_steering_left,
    initial_length_steering_right,
    steering_tape_final_extension,
):
    """
    Apply steering actuation by changing the relative rest lengths of the two steering tapes.

    The repository convention is to shorten the left steering tape and lengthen the right
    steering tape by the same amount for a positive final extension.
    """
    if np.abs(steering_tape_extension_step) <= 1e-9:
        return False
    if np.abs(steering_tape_final_extension) <= 1e-9:
        return False

    if steering_tape_indices is None or len(steering_tape_indices) < 2:
        raise ValueError("steering_tape_indices must contain at least two entries.")

    steering_tape_left_index = int(steering_tape_indices[0])
    steering_tape_right_index = int(steering_tape_indices[1])

    desired_length_steering_left = (
        initial_length_steering_left - steering_tape_final_extension
    )
    desired_length_steering_right = (
        initial_length_steering_right + steering_tape_final_extension
    )

    if config["structural_solver"] == "pss":
        psystem.update_rest_length(
            steering_tape_left_index,
            desired_length_steering_left
            - float(psystem.extract_rest_length[steering_tape_left_index]),
        )
        psystem.update_rest_length(
            steering_tape_right_index,
            desired_length_steering_right
            - float(psystem.extract_rest_length[steering_tape_right_index]),
        )
    elif config["structural_solver"] == "kite_fem":
        if kite_connectivity_arr is None:
            raise ValueError(
                "kite_connectivity_arr is required for kite_fem steering actuation."
            )

        steering_tape_left_spring_id = _find_kite_fem_spring_id_from_connectivity(
            kite_fem_structure=kite_fem_structure,
            kite_connectivity_arr=kite_connectivity_arr,
            connectivity_idx=steering_tape_left_index,
        )
        steering_tape_right_spring_id = _find_kite_fem_spring_id_from_connectivity(
            kite_fem_structure=kite_fem_structure,
            kite_connectivity_arr=kite_connectivity_arr,
            connectivity_idx=steering_tape_right_index,
        )
        kite_fem_structure.modify_get_spring_rest_length(
            spring_ids=[steering_tape_left_spring_id, steering_tape_right_spring_id],
            new_l0s=[desired_length_steering_left, desired_length_steering_right],
        )
    else:
        raise ValueError("Invalid structural solver specified, either pss or kite_fem")

    logging.info(
        "Steering tapes updated: left %.3fm -> %.3fm, right %.3fm -> %.3fm",
        initial_length_steering_left,
        desired_length_steering_left,
        initial_length_steering_right,
        desired_length_steering_right,
    )
    return True


def update_power_tape_actuation(
    config,
    psystem,
    kite_fem_structure,
    kite_connectivity_arr,
    power_tape_index,
    power_tape_extension_step,
    initial_length_power_tape,
    power_tape_final_extension,
    should_apply_update,
    n_power_tape_steps,
    rest_lengths=None,
):
    """
    Calculate current power tape extension and update if needed for actuation.

    Args:
        config: Configuration dictionary
        psystem: Particle system (for PSS solver)
        kite_fem_structure: FEM structure (for kite_fem solver)
        kite_connectivity_arr: ASKITE connectivity array
        power_tape_index: Index of power tape in connectivity array
        power_tape_extension_step: Increment for power tape extension
        initial_length_power_tape: Initial length of power tape
        power_tape_final_extension: Final desired power tape extension
        should_apply_update: Flag indicating if a depower step should be applied now
        n_power_tape_steps: Number of power tape extension steps
        rest_lengths: Current rest lengths array (for kite_fem solver)

    Returns:
        tuple: (delta_power_tape, is_actuation_finalized, did_update)
            - delta_power_tape: Current change in power tape length
            - is_actuation_finalized: True if actuation is complete, False otherwise
            - did_update: True if a depower step was applied this iteration
    """
    is_actuation_finalized = True
    did_update = False

    ## Calculate delta tape lengths based on structural solver
    if config["structural_solver"] == "pss":
        current_length = float(psystem.extract_rest_length[power_tape_index])
        delta_power_tape = current_length - initial_length_power_tape
        _, should_update = _compute_power_tape_increment(
            delta_power_tape=delta_power_tape,
            power_tape_final_extension=power_tape_final_extension,
            power_tape_extension_step=power_tape_extension_step,
        )
        is_actuation_finalized = not should_update

        if should_apply_update and should_update:
            increment, _ = _compute_power_tape_increment(
                delta_power_tape=delta_power_tape,
                power_tape_final_extension=power_tape_final_extension,
                power_tape_extension_step=power_tape_extension_step,
            )
            psystem.update_rest_length(power_tape_index, increment)
            did_update = True
            current_length = float(psystem.extract_rest_length[power_tape_index])
            delta_power_tape = current_length - initial_length_power_tape
            _, should_update_after = _compute_power_tape_increment(
                delta_power_tape=delta_power_tape,
                power_tape_final_extension=power_tape_final_extension,
                power_tape_extension_step=power_tape_extension_step,
            )
            is_actuation_finalized = not should_update_after
            logging.info(
                f"||--- delta l_d: {delta_power_tape:.3f}m | new l_d: {current_length:.3f}m | Steps required: {n_power_tape_steps}"
            )

    elif config["structural_solver"] == "kite_fem":
        if kite_connectivity_arr is None:
            raise ValueError(
                "kite_connectivity_arr is required for kite_fem power tape actuation."
            )

        spring_id = _find_kite_fem_spring_id_from_connectivity(
            kite_fem_structure=kite_fem_structure,
            kite_connectivity_arr=kite_connectivity_arr,
            connectivity_idx=power_tape_index,
        )
        current_length = float(kite_fem_structure.spring_elements[spring_id].l0)
        delta_power_tape = current_length - initial_length_power_tape
        _, should_update = _compute_power_tape_increment(
            delta_power_tape=delta_power_tape,
            power_tape_final_extension=power_tape_final_extension,
            power_tape_extension_step=power_tape_extension_step,
        )
        is_actuation_finalized = not should_update

        if should_apply_update and should_update:
            increment, _ = _compute_power_tape_increment(
                delta_power_tape=delta_power_tape,
                power_tape_final_extension=power_tape_final_extension,
                power_tape_extension_step=power_tape_extension_step,
            )
            new_length = current_length + increment
            kite_fem_structure.modify_get_spring_rest_length(
                spring_ids=[spring_id],
                new_l0s=[new_length],
            )
            did_update = True
            delta_power_tape = new_length - initial_length_power_tape
            _, should_update_after = _compute_power_tape_increment(
                delta_power_tape=delta_power_tape,
                power_tape_final_extension=power_tape_final_extension,
                power_tape_extension_step=power_tape_extension_step,
            )
            is_actuation_finalized = not should_update_after
            logging.info(
                f"||--- delta l_d: {delta_power_tape:.3f}m | new l_d: {new_length:.3f}m | Steps required: {n_power_tape_steps}"
            )

    return delta_power_tape, is_actuation_finalized, did_update


def update_steering_tape_actuation_progressive(
    config,
    psystem,
    kite_fem_structure,
    kite_connectivity_arr,
    steering_tape_indices,
    steering_tape_extension_step,
    initial_length_steering_left,
    initial_length_steering_right,
    steering_tape_final_extension,
    should_apply_update,
):
    """Progressively apply steering actuation, mirroring depower internal stepping."""
    if steering_tape_indices is None or len(steering_tape_indices) < 2:
        return 0.0, True, False

    target = float(steering_tape_final_extension)
    if np.abs(target) <= 1e-9:
        return 0.0, True, False

    left_idx = int(steering_tape_indices[0])
    right_idx = int(steering_tape_indices[1])

    if config["structural_solver"] == "pss":
        current_left = float(psystem.extract_rest_length[left_idx])
        current_right = float(psystem.extract_rest_length[right_idx])
    elif config["structural_solver"] == "kite_fem":
        if kite_connectivity_arr is None:
            raise ValueError(
                "kite_connectivity_arr is required for kite_fem steering actuation."
            )
        left_spring_id = _find_kite_fem_spring_id_from_connectivity(
            kite_fem_structure=kite_fem_structure,
            kite_connectivity_arr=kite_connectivity_arr,
            connectivity_idx=left_idx,
        )
        right_spring_id = _find_kite_fem_spring_id_from_connectivity(
            kite_fem_structure=kite_fem_structure,
            kite_connectivity_arr=kite_connectivity_arr,
            connectivity_idx=right_idx,
        )
        current_left = float(kite_fem_structure.spring_elements[left_spring_id].l0)
        current_right = float(kite_fem_structure.spring_elements[right_spring_id].l0)
    else:
        raise ValueError("Invalid structural solver specified, either pss or kite_fem")

    current_delta_left = float(initial_length_steering_left) - current_left
    current_delta_right = current_right - float(initial_length_steering_right)
    current_delta = 0.5 * (current_delta_left + current_delta_right)

    if np.abs(target - current_delta) <= 1e-9:
        return current_delta, True, False

    if not should_apply_update:
        return current_delta, False, False

    step = float(steering_tape_extension_step)
    if np.abs(step) <= 1e-9:
        increment = target - current_delta
    else:
        remaining = target - current_delta
        increment = np.sign(remaining) * min(np.abs(step), np.abs(remaining))

    next_delta = current_delta + increment
    update_steering_tape_actuation(
        config=config,
        psystem=psystem,
        kite_fem_structure=kite_fem_structure,
        kite_connectivity_arr=kite_connectivity_arr,
        steering_tape_indices=steering_tape_indices,
        steering_tape_extension_step=step,
        initial_length_steering_left=initial_length_steering_left,
        initial_length_steering_right=initial_length_steering_right,
        steering_tape_final_extension=next_delta,
    )

    is_finalized = np.abs(target - next_delta) <= 1e-9
    return next_delta, is_finalized, True


def compute_adaptive_dt(
    f_residual_list,
    dt_initial,
    dt_max,
    residual_tol,
):
    """
    Compute adaptive dt that increases as residual approaches convergence.

    Strategy: dt scales with the absolute residual level relative to the
    convergence tolerance. Residual near tolerance -> dt approaches dt_max.

    Args:
        f_residual_list: List of residual force norms over iterations
        dt_initial: Initial dt from config
        dt_max: Maximum dt (cap to prevent instability)
        residual_tol: Absolute residual convergence tolerance

    Returns:
        dt: Adaptive timestep value
    """
    if len(f_residual_list) < 1:
        return dt_initial

    current_residual = f_residual_list[-1]

    # Defensive handling for invalid settings/values.
    if (
        residual_tol is None
        or residual_tol <= 0
        or not np.isfinite(residual_tol)
        or not np.isfinite(current_residual)
        or current_residual < 0
    ):
        return dt_initial

    # Scale from dt_initial (large residuals) to dt_max (near tolerance).
    # ratio = 1 at/under tolerance, tends to 0 for large residuals.
    ratio = np.clip(residual_tol / max(current_residual, 1e-12), 0.0, 1.0)
    dt_adaptive = dt_initial + (dt_max - dt_initial) * ratio

    # Ensure dt remains bounded regardless of dt ordering.
    dt_low = min(dt_initial, dt_max)
    dt_high = max(dt_initial, dt_max)
    dt = float(np.clip(dt_adaptive, dt_low, dt_high))

    return dt


def distribute_total_force_by_particle_mass(total_force, m_arr):
    """Distribute a total 3D force over nodes proportional to particle masses."""
    total_force = np.asarray(total_force, dtype=float).reshape(3)
    masses = np.asarray(m_arr, dtype=float).reshape(-1)
    mass_sum = float(np.sum(masses))

    if mass_sum <= 1e-12:
        raise ValueError("Total particle mass must be positive to distribute force.")

    mass_fraction = masses / mass_sum
    return mass_fraction[:, None] * total_force[None, :]


# TODO: this should also use structural is not converging
def check_convergence(
    i,
    f_residual,
    f_residual_list,
    f_aero_wing_vsm_format,
    config,
    stagnation_check_start=0,
):
    """
    Check convergence conditions for the aero-structural solver.

    Args:
        i: Current iteration number
        f_residual: Current residual force vector
        f_residual_list: List of residual force norms from all iterations
        f_aero_wing_vsm_format: Aerodynamic forces in VSM format
        config: Configuration dictionary
        stagnation_check_start: Iteration index from which to check stagnation
            (reset when switching regularization phase)

    Returns:
        tuple: (is_convergence, should_break, is_stagnated)
            - is_convergence: True if converged, False otherwise
            - should_break: True if loop should break, False to continue
            - is_stagnated: True if residual has stagnated (no longer changing)
    """
    is_convergence = False
    should_break = False
    is_stagnated = False

    n_stag = config["aero_structural_solver"]["n_max_constant_residual_force"]
    # Number of iterations since the stagnation check window started
    iters_since_start = i - stagnation_check_start

    ### All the convergence checks, are be done in if-elif because only 1 should hold at once
    # if convergence (residual below set tolerance)
    if np.linalg.norm(f_residual) <= config["aero_structural_solver"]["tol"]:
        is_convergence = True

    # if residual forces are NaN
    elif np.isnan(np.linalg.norm(f_residual)):
        is_convergence = False
        logging.info("Classic PS diverged - residual force is NaN")
        should_break = True

    # if residual forces are not changing anymore across the full recent window
    elif iters_since_start >= n_stag and n_stag > 0:
        window_vals = np.asarray(f_residual_list[i - n_stag : i + 1], dtype=float)
        if window_vals.size > 0 and np.isfinite(window_vals).all():
            residual_span = float(np.max(window_vals) - np.min(window_vals))
            if residual_span < config["aero_structural_solver"]["stagnation_tol"]:
                is_convergence = False
                is_stagnated = True

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

    return is_convergence, should_break, is_stagnated


def main(
    m_arr=None,
    struc_nodes=None,
    struc_nodes_initial=None,
    system_model=None,
    config=None,
    ### ACTUATION
    initial_length_power_tape=None,
    n_power_tape_steps=None,
    power_tape_final_extension=None,
    power_tape_extension_step=None,
    initial_length_steering_left=None,
    initial_length_steering_right=None,
    steering_tape_indices=None,
    steering_tape_final_extension=None,
    steering_tape_extension_step=None,
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
    print(f'--> Running structural_solver: {config["structural_solver"]}')

    ## PRELOOP
    f_ext_gravity = np.zeros(struc_nodes.shape)

    if config["structural_solver"] == "kite_fem":
        rest_lengths = kite_fem_structure.modify_get_spring_rest_length()

    max_iter = config["aero_structural_solver"]["max_iter"]
    # Keep index 0 for the pre-loop initial state and reserve max_iter loop slots.
    t_vector = np.linspace(0, max_iter, max_iter + 1)
    tracking_data = tracking.setup_tracking_arrays(len(struc_nodes), t_vector)
    is_convergence = False
    f_residual_list = []
    f_tether_drag = np.zeros(3)
    is_actuation_finalized = True
    is_steering_finalized = True
    struc_nodes_prev = None  # Initialize previous points for tracking
    start_time = time.time()
    plotting.set_plot_style()

    # Two-phase regularization: phase 1 = with pseudo_dt, phase 2 = without
    reg_phase = 1  # 1 = regularized, 2 = unregularized (polish)
    stagnation_check_start = 0  # iteration at which current phase started

    # Adaptive dt for PSS solver
    dt_initial = config["structural_pss"]["dt"]
    dt_max = config["structural_pss"].get(
        "dt_max", dt_initial * 10.0
    )  # Default to 10x initial dt

    # Quasi-steady stagnation stop: if rounded opt_x stops changing for N iterations.
    qs_stag_decimals = int(
        config["aero_structural_solver"].get("qs_state_stagnation_decimals", 3)
    )
    qs_stag_n_iter = int(
        config["aero_structural_solver"].get("qs_state_stagnation_n_iter", 0)
    )
    steering_actuation_interval_iters = int(
        config["aero_structural_solver"].get("steering_actuation_interval_iters", 5)
    )
    steering_actuation_interval_iters = max(1, steering_actuation_interval_iters)
    power_tape_actuation_interval_iters = int(
        config["aero_structural_solver"].get(
            "power_tape_actuation_interval_iters", steering_actuation_interval_iters
        )
    )
    power_tape_actuation_interval_iters = max(1, power_tape_actuation_interval_iters)
    depower_settle_iterations_after_update = int(
        config["aero_structural_solver"].get(
            "depower_settle_iterations_after_update", 2
        )
    )
    depower_settle_iterations_after_update = max(
        0, depower_settle_iterations_after_update
    )
    qs_opt_prev_rounded = None
    qs_stag_counter = 0
    qs_state_should_break = False
    depower_settle_counter = 0
    logging.info(
        "Steering actuation interval: every %s iterations",
        steering_actuation_interval_iters,
    )
    logging.info(
        "Depower actuation interval: every %s iterations",
        power_tape_actuation_interval_iters,
    )
    logging.info(
        "Depower settle iterations after update: %s",
        depower_settle_iterations_after_update,
    )

    bridle_node_pairs = None

    if config["is_with_aero_bridle"]:
        bridle_node_pairs = (
            aerodynamic_bridle_line_drag.build_bridle_node_pairs_from_line_system(
                struc_nodes,
                getattr(body_aero, "_bridle_line_system", None),
            )
        )

    # Aitken relaxation state
    omega_relaxation = config["aero_structural_solver"].get("relaxation_factor", 0.3)
    r_prev_flat = None

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
    # Initialization of external forces pre-simulation loop
    ######################################################################

    ### STRUC --> AERO
    le_arr, te_arr = struc2aero.main(
        struc_nodes,
        struc_node_le_indices,
        struc_node_te_indices,
        config["aerodynamic"]["n_aero_panels_per_struc_section"],
    )

    cg = calculate_cg(struc_nodes=struc_nodes, m_arr=m_arr)
    ### AERO
    f_aero_wing_vsm_format, body_aero, results_aero = aerodynamic_vsm.run_vsm_package(
        body_aero=body_aero,
        solver=vsm_solver,
        system_model=system_model,
        center_of_gravity=cg,
        le_arr=le_arr,
        te_arr=te_arr,
        # va_vector=vel_app,
        aero_input_type="reuse_initial_polar_data",
        initial_polar_data=initial_polar_data,
        include_gravity=config["is_with_gravity"],
        is_with_plot=config["is_with_aero_plot_per_iteration"],
    )
    logging.debug(
        f"Aero symmetry check, f_aero_y: {np.sum([force[1] for force in f_aero_wing_vsm_format])}"
    )
    roll, pitch, yaw = results_aero["opt_x"][1:4]
    struc_nodes = rotate_geometry(struc_nodes, angle_deg=[roll, pitch, yaw])
    ### AERO --> STRUC
    f_aero_wing = aero2struc_level_1.main(
        config["aero2struc"]["coupling_method"],
        f_aero_wing_vsm_format,
        struc_nodes,
        np.array(results_aero["panel_cp_locations"]),
        aero2struc_mapping,
        config["is_with_coupling_plot_per_iteration"],
        config["aero2struc"],
    )

    # Check moment preservation of aero→struc mapping (pre-loop)
    aero2struc_level_1.check_moment_preservation(
        f_aero_panel=f_aero_wing_vsm_format,
        panel_cps=np.array(results_aero["panel_cp_locations"]),
        f_aero_mapped=f_aero_wing,
        struc_nodes=struc_nodes,
    )

    ### BRIDLE AERO
    # f_aero_bridle = aerodynamic_bridle_line_drag.main(
    #     struc_nodes,
    #     bridle_connectivity_arr,
    #     bridle_diameter_arr,
    #     vel_app,
    #     config["rho"],
    #     config["aerodynamic_bridle"]["cd_cable"],
    #     config["aerodynamic_bridle"]["cf_cable"],
    # )
    f_aero_bridle = np.zeros((len(struc_nodes), 3))
    f_inertial = distribute_total_force_by_particle_mass(
        results_aero.get("inertial_force", np.zeros(3)),
        m_arr,
    )
    f_ext_gravity = distribute_total_force_by_particle_mass(
        results_aero.get("gravity_force", np.zeros(3)),
        m_arr,
    )
    f_aero = f_aero_wing + f_aero_bridle
    ## EXTERNAL FORCE
    f_ext = f_aero + f_inertial + f_ext_gravity
    f_ext = np.round(f_ext, 5)
    f_ext_flat = f_ext.flatten()

    ######################################################################
    # SIMULATION LOOP
    ######################################################################
    ## propagating the simulation for each timestep and saving results
    with tqdm(total=max_iter, desc="Simulating", leave=True) as pbar:
        for i in range(max_iter):
            if i > 0:
                struc_nodes_prev = struc_nodes.copy()

            ########################################################
            ############## INTERNAL FORCE CALCULATION ##############
            ########################################################
            begin_time_f_int = time.time()
            if config["structural_solver"] == "pss":
                # Apply adaptive dt based on convergence progress
                if len(f_residual_list) > 0:
                    adaptive_dt = compute_adaptive_dt(
                        f_residual_list,
                        dt_initial,
                        dt_max,
                        config["aero_structural_solver"]["tol"],
                    )
                    config["structural_pss"]["dt"] = adaptive_dt
                    logging.debug(
                        f"Adaptive dt updated: {adaptive_dt:.6f} (residual: {f_residual_list[-1]:.3f}N)"
                    )
                    print(f"Adaptive dt: {adaptive_dt:.6f} s at iteration {i}")
                psystem, is_structural_converged, struc_nodes, f_int = (
                    structural_pss.run_pss(
                        psystem,
                        f_ext_flat,
                        config["structural_pss"],
                    )
                )
            elif config["structural_solver"] == "kite_fem":
                kite_fem_structure, is_structural_converged, struc_nodes, f_int = (
                    structural_kite_fem_level_1.run_kite_fem(
                        kite_fem_structure, f_ext_flat, config["structural_kite_fem"]
                    )
                )
            end_time_f_int = time.time()

            ### Aitken relaxation of structural nodes
            if struc_nodes_prev is not None:
                r_k = struc_nodes - struc_nodes_prev
                r_k_flat = r_k.flatten()

                if (
                    config["aero_structural_solver"].get(
                        "is_with_aitken_relaxation", True
                    )
                    and r_prev_flat is not None
                ):
                    delta_r = r_k_flat - r_prev_flat
                    denom = np.dot(delta_r, delta_r)
                    if denom > 1e-30:
                        omega_relaxation = -omega_relaxation * (
                            np.dot(r_prev_flat, delta_r) / denom
                        )
                        omega_relaxation = np.clip(omega_relaxation, 0.05, 1.0)

                struc_nodes = struc_nodes_prev + omega_relaxation * r_k
                r_prev_flat = r_k_flat.copy()
                logging.debug(f"Aitken relaxation omega: {omega_relaxation:.4f}")

                # Sync relaxed positions back to structural solver state
                if config["structural_solver"] == "pss":
                    for idx, particle in enumerate(psystem.particles):
                        particle.update_pos(struc_nodes[idx])
                        particle.update_vel(np.zeros(3))
                elif config["structural_solver"] == "kite_fem":
                    # Update kite_fem so the next solve() starts from the
                    # Aitken-relaxed geometry instead of the original construction
                    # geometry.  coords_rotations_init is the reference that
                    # solve() adds displacements to, so moving it here makes the
                    # Newton-Raphson start near the current state.
                    flat_xyz = struc_nodes.flatten()
                    kite_fem_structure.coords_current = flat_xyz.copy()
                    # Build the 6-DOF vector [x,y,z, 0,0,0] per node
                    n_nodes = len(struc_nodes)
                    coords_rot = np.zeros(n_nodes * 6)
                    for ni in range(n_nodes):
                        coords_rot[6 * ni : 6 * ni + 3] = struc_nodes[ni]
                    kite_fem_structure.coords_rotations_init = coords_rot.copy()
                    kite_fem_structure.coords_rotations_current = coords_rot.copy()

            ### PLOT per iteration
            if config["is_with_struc_plot_per_iteration"]:
                if config["structural_solver"] == "pss":
                    rest_lengths = psystem.extract_rest_length
                elif config["structural_solver"] == "kite_fem":
                    rest_lengths = structural_kite_fem_level_1.get_rest_lengths(
                        kite_fem_structure, kite_connectivity_arr
                    )
                    kite_fem_structure.plot_convergence()

                plotting.main(
                    struc_nodes,
                    kite_connectivity_arr,
                    rest_lengths,
                    f_ext=f_ext,
                    f_bridle=f_aero_bridle if config["is_with_aero_bridle"] else None,
                    f_inertial=f_inertial,
                    title=f"i: {i}",
                    body_aero=body_aero,
                    is_with_node_indices=False,
                    pulley_line_indices=pulley_line_indices,
                    pulley_line_to_other_node_pair_dict=pulley_line_to_other_node_pair_dict,
                )

            ########################################################
            ############## INTERNAL FORCE CALCULATION ##############
            ########################################################
            begin_time_f_ext = time.time()

            ### STRUC --> AERO
            le_arr, te_arr = struc2aero.main(
                struc_nodes,
                struc_node_le_indices,
                struc_node_te_indices,
                config["aerodynamic"]["n_aero_panels_per_struc_section"],
            )

            cg = calculate_cg(struc_nodes=struc_nodes, m_arr=m_arr)
            ### AERO
            f_aero_wing_vsm_format, body_aero, results_aero = (
                aerodynamic_vsm.run_vsm_package(
                    body_aero=body_aero,
                    solver=vsm_solver,
                    system_model=system_model,
                    center_of_gravity=cg,
                    le_arr=le_arr,
                    te_arr=te_arr,
                    current_guess=[
                        results_aero["opt_x"][0],
                        0,
                        0,
                        0,
                        results_aero["opt_x"][4],
                    ],
                    # va_vector=vel_app,
                    aero_input_type="reuse_initial_polar_data",
                    initial_polar_data=initial_polar_data,
                    include_gravity=config["is_with_gravity"],
                    is_with_plot=config["is_with_aero_plot_per_iteration"],
                )
            )
            logging.debug(
                f"Aero symmetry check, f_aero_y: {np.sum([force[1] for force in f_aero_wing_vsm_format])}"
            )
            print("Quasi-steady state solver info:")
            # print(f"  Converged: {results_aero['is_converged']}")
            print(f"  Kite_speed: {results_aero['opt_x'][0]:.2f} m/s")
            print(f"  Roll: {results_aero['opt_x'][1]:.2f} deg")
            print(f"  Pitch: {results_aero['opt_x'][2]:.2f} deg")
            print(f"  Yaw: {results_aero['opt_x'][3]:.2f} deg")
            print(f"  Course rate: {results_aero['opt_x'][4]:.2f} rad/s")

            # Stop if quasi-steady state vector has effectively frozen.
            if qs_stag_n_iter > 0:
                qs_opt_current = np.asarray(results_aero.get("opt_x", []), dtype=float)
                if qs_opt_current.size > 0:
                    qs_opt_current_rounded = np.round(qs_opt_current, qs_stag_decimals)
                    if (
                        qs_opt_prev_rounded is not None
                        and qs_opt_prev_rounded.shape == qs_opt_current_rounded.shape
                        and np.array_equal(qs_opt_prev_rounded, qs_opt_current_rounded)
                    ):
                        qs_stag_counter += 1
                    else:
                        qs_stag_counter = 0

                    qs_opt_prev_rounded = qs_opt_current_rounded.copy()

                    if qs_stag_counter >= qs_stag_n_iter:
                        qs_state_should_break = True
                        logging.info(
                            "Stopping: quasi-steady opt_x unchanged up to %s decimals for %s consecutive iterations. "
                            "opt_x_rounded=%s",
                            qs_stag_decimals,
                            qs_stag_n_iter,
                            qs_opt_current_rounded,
                        )
            roll, pitch, yaw = results_aero["opt_x"][1:4]
            struc_nodes = rotate_geometry(struc_nodes, angle_deg=[roll, pitch, yaw])
            ### AERO --> STRUC
            f_aero_wing = aero2struc_level_1.main(
                config["aero2struc"]["coupling_method"],
                f_aero_wing_vsm_format,
                struc_nodes,
                np.array(results_aero["panel_cp_locations"]),
                aero2struc_mapping,
                config["is_with_coupling_plot_per_iteration"],
                config["aero2struc"],
            )

            # Check moment preservation (only first coupling iteration to limit log spam)
            if i == 1:
                aero2struc_level_1.check_moment_preservation(
                    f_aero_panel=f_aero_wing_vsm_format,
                    panel_cps=np.array(results_aero["panel_cp_locations"]),
                    f_aero_mapped=f_aero_wing,
                    struc_nodes=struc_nodes,
                )

            ### BRIDLE AERO
            if config["is_with_aero_bridle"]:
                f_aero_bridle = aerodynamic_bridle_line_drag.main(
                    struc_nodes,
                    bridle_connectivity_arr,
                    bridle_diameter_arr,
                    vel_app,
                    config["rho"],
                    config["aerodynamic_bridle"]["cd_cable"],
                    config["aerodynamic_bridle"]["cf_cable"],
                    body_aero=body_aero,
                    bridle_node_pairs=bridle_node_pairs,
                )
                bridle_force_total = np.sum(f_aero_bridle, axis=0)
                bridle_force_total_norm = np.linalg.norm(bridle_force_total)
                bridle_force_nodal_max = np.max(np.linalg.norm(f_aero_bridle, axis=1))
                print(
                    "Bridle aero force: "
                    f"|sum|={bridle_force_total_norm:.3f}N "
                    f"sum=[{bridle_force_total[0]:.3f}, {bridle_force_total[1]:.3f}, {bridle_force_total[2]:.3f}]N "
                    f"max_node={bridle_force_nodal_max:.3f}N"
                )
                logging.info(
                    "Bridle aero force iter %s: |sum|=%.3fN, sum=(%.3f, %.3f, %.3f)N, max_node=%.3fN",
                    i,
                    bridle_force_total_norm,
                    bridle_force_total[0],
                    bridle_force_total[1],
                    bridle_force_total[2],
                    bridle_force_nodal_max,
                )
            else:
                f_aero_bridle = np.zeros((len(struc_nodes), 3))

            inertial_force_total = np.asarray(
                results_aero.get("inertial_force", np.zeros(3)), dtype=float
            )
            f_inertial = distribute_total_force_by_particle_mass(
                inertial_force_total,
                m_arr,
            )
            gravity_force_total = np.asarray(
                results_aero.get("gravity_force", np.zeros(3)), dtype=float
            )
            f_ext_gravity = distribute_total_force_by_particle_mass(
                gravity_force_total,
                m_arr,
            )
            inertial_force_total_norm = np.linalg.norm(inertial_force_total)
            print(
                "Inertial force (QSM): "
                f"|sum|={inertial_force_total_norm:.3f}N "
                f"sum=[{inertial_force_total[0]:.3f}, {inertial_force_total[1]:.3f}, {inertial_force_total[2]:.3f}]N"
            )
            gravity_force_total_norm = np.linalg.norm(gravity_force_total)
            print(
                "Gravity force (QSM): "
                f"|sum|={gravity_force_total_norm:.3f}N "
                f"sum=[{gravity_force_total[0]:.3f}, {gravity_force_total[1]:.3f}, {gravity_force_total[2]:.3f}]N"
            )
            logging.info(
                "Inertial force iter %s: |sum|=%.3fN, sum=(%.3f, %.3f, %.3f)N",
                i,
                inertial_force_total_norm,
                inertial_force_total[0],
                inertial_force_total[1],
                inertial_force_total[2],
            )
            logging.info(
                "Gravity force iter %s: |sum|=%.3fN, sum=(%.3f, %.3f, %.3f)N",
                i,
                gravity_force_total_norm,
                gravity_force_total[0],
                gravity_force_total[1],
                gravity_force_total[2],
            )
            f_aero = f_aero_wing + f_aero_bridle

            ## EXTERNAL FORCE
            f_ext = f_aero + f_inertial + f_ext_gravity
            f_ext = np.round(f_ext, 5)
            f_ext_flat = f_ext.flatten()
            end_time_f_ext = time.time()

            ### FORCING SYMMETRY
            if config["is_with_forcing_symmetry"]:
                logging.info("Forcing symmetry in y-direction")
                struc_nodes = forcing_symmetry(struc_nodes)

            ### RESIDUAL
            f_residual = f_int + f_ext_flat

            # Zero out residual at fixed (constrained) nodes — their imbalance
            # is carried by the constraint reaction force, not by f_int.
            # Without this, the residual includes e.g. the weight of node 0
            # (~92 N) which can never converge to zero.
            if config["structural_solver"] == "pss":
                for fix_idx in config["structural_pss"]["fixed_point_indices"]:
                    f_residual[3 * fix_idx : 3 * fix_idx + 3] = 0.0

            f_residual_list.append(np.linalg.norm(np.abs(f_residual)))
            if config["structural_solver"] == "pss":
                logging.debug(
                    f"residual force in y-direction: {np.sum([f_residual[1::3]]):.3f}N"
                )

            ### TRACKING
            # Update unified tracking dataframe (replaces position update)
            # Use i+1 so that positions[0] retains the true initial geometry
            # stored in the pre-loop call.
            tracking.update_tracking_arrays(
                tracking_data,
                i + 1,
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
            is_convergence, should_break, is_stagnated = check_convergence(
                i=i,
                f_residual=f_residual,
                f_residual_list=f_residual_list,
                f_aero_wing_vsm_format=f_aero_wing_vsm_format,
                config=config,
                stagnation_check_start=stagnation_check_start,
            )

            should_apply_steering_now = (
                steering_tape_extension_step != 0
                and steering_tape_final_extension != 0
                and ((i + 1) % steering_actuation_interval_iters == 0)
            )

            delta_steering, is_steering_finalized, did_update_steering = (
                update_steering_tape_actuation_progressive(
                    config=config,
                    psystem=psystem,
                    kite_fem_structure=kite_fem_structure,
                    kite_connectivity_arr=kite_connectivity_arr,
                    steering_tape_indices=steering_tape_indices,
                    steering_tape_extension_step=steering_tape_extension_step,
                    initial_length_steering_left=initial_length_steering_left,
                    initial_length_steering_right=initial_length_steering_right,
                    steering_tape_final_extension=steering_tape_final_extension,
                    should_apply_update=should_apply_steering_now,
                )
            )

            # Two-phase regularization: on stagnation in phase 1, disable
            # pseudo_dt and continue to let the solver polish to true equilibrium.
            if is_stagnated:
                if reg_phase == 1 and config["structural_solver"] == "kite_fem":
                    reg_phase = 2
                    stagnation_check_start = i  # reset stagnation window
                    config["structural_kite_fem"]["pseudo_dt"] = None
                    logging.info(
                        f"Phase 1 stagnated at iter {i} (res={np.linalg.norm(f_residual):.1f}N). "
                        f"Switching to phase 2: pseudo_dt=None (no regularization)."
                    )
                else:
                    logging.info(
                        "Classic PS non-converging - residual no longer changes"
                    )
                    should_break = True

            if qs_state_should_break:
                # Do not allow quasi-steady stagnation stopping to interrupt
                # progressive actuation before targets are fully applied.
                if (
                    is_actuation_finalized
                    and is_steering_finalized
                    and depower_settle_counter == 0
                ):
                    break
                qs_state_should_break = False

            ### ACTUATION (depower & steering checked every iteration; applied at enforced cadence)
            should_apply_depower_now = (
                power_tape_extension_step != 0
                and power_tape_final_extension != 0
                and ((i + 1) % power_tape_actuation_interval_iters == 0)
            )
            (
                delta_power_tape,
                is_actuation_finalized,
                did_update_depower,
            ) = update_power_tape_actuation(
                config=config,
                psystem=psystem,
                kite_fem_structure=kite_fem_structure,
                kite_connectivity_arr=kite_connectivity_arr,
                power_tape_index=power_tape_index,
                power_tape_extension_step=power_tape_extension_step,
                initial_length_power_tape=initial_length_power_tape,
                power_tape_final_extension=power_tape_final_extension,
                should_apply_update=should_apply_depower_now,
                n_power_tape_steps=n_power_tape_steps,
                rest_lengths=(
                    rest_lengths if config["structural_solver"] == "kite_fem" else None
                ),
            )

            if did_update_depower:
                depower_settle_counter = depower_settle_iterations_after_update
            elif depower_settle_counter > 0:
                depower_settle_counter -= 1

            # If either actuation is not finalized, continue actuation phase.
            if (
                (not is_actuation_finalized)
                or (not is_steering_finalized)
                or (depower_settle_counter > 0)
            ):
                continue

            # Check if we should exit the loop
            if should_break or (
                is_convergence and is_actuation_finalized and is_steering_finalized
            ):
                break
    ######################################################################
    ## END OF SIMULATION FOR LOOP
    ######################################################################

    # print out the geometric angle of attack of the mid panel
    panels = body_aero.panels

    # Select middle panel
    mid_idx = len(panels) // 2
    panel = panels[mid_idx]

    # Midpoints of leading and trailing edges
    le_mid = 0.5 * (panel.LE_point_1 + panel.LE_point_2)
    te_mid = 0.5 * (panel.TE_point_1 + panel.TE_point_2)

    # Chord direction vector
    vec_chord = te_mid - le_mid
    vec_chord /= np.linalg.norm(vec_chord)

    # Apparent wind direction (normalize)
    vec_wind = body_aero.va
    # Project onto plane of interest (optional: usually x-z plane)
    # Remove spanwise component if needed
    vec_chord_2d = np.array([vec_chord[0], vec_chord[2]])
    vec_wind_2d = np.array([vec_wind[0], vec_wind[2]])

    vec_chord_2d /= np.linalg.norm(vec_chord_2d)
    vec_wind_2d /= np.linalg.norm(vec_wind_2d)

    # Angle between vectors (signed)
    dot = np.clip(np.dot(vec_chord_2d, vec_wind_2d), -1.0, 1.0)
    cross = np.cross(vec_chord_2d, vec_wind_2d)

    angle = np.arctan2(cross, dot)
    alpha_at_ac_mid = np.ravel(results_aero["alpha_at_ac"])[mid_idx]

    print(f"alpha = {np.degrees(angle):.2f}° (va vs mid-span chord)")
    print(
        f'alpha = {np.rad2deg(alpha_at_ac_mid):.2f}° (incl. induced velocity, from results_aero["alpha_at_ac"])'
    )
    # print(
    #     f'results_aero["alpha_uncorrected"]: {float(np.rad2deg(results_aero["alpha_uncorrected"][mid_idx])):.2f}°'
    # )
    # print(
    #     f'results_aero["alpha_geometric"]: wrt horizontal {results_aero["alpha_geometric"][mid_idx]:.2f}°'
    # )

    if config["structural_solver"] == "pss":
        rest_lengths = psystem.extract_rest_length
    elif config["structural_solver"] == "kite_fem":
        rest_lengths = structural_kite_fem_level_1.get_rest_lengths(
            kite_fem_structure, kite_connectivity_arr
        )

    if config["is_with_final_plot"]:
        plotting.main(
            struc_nodes,
            kite_connectivity_arr,
            f_ext=f_ext,
            f_bridle=f_aero_bridle if config["is_with_aero_bridle"] else None,
            f_inertial=f_inertial,
            rest_lengths=rest_lengths,
            struc_nodes_initial=struc_nodes_initial,
            title="Initial vs final",
            pulley_line_indices=pulley_line_indices,
            pulley_line_to_other_node_pair_dict=pulley_line_to_other_node_pair_dict,
            vel_app=vel_app,
        )

    # Calculate cl, cd, tether force, and va for output
    opt_x = np.asarray(results_aero.get("opt_x", np.full(5, np.nan)), dtype=float)
    kite_speed = opt_x[0] if opt_x.size > 0 else np.nan

    # Apparent wind speed from quasi-steady kinematics (opt_x[0]).
    # Fallback to body_aero.va for backward compatibility.
    if opt_x.size > 0 and np.isfinite(opt_x[0]):
        va = float(opt_x[0])
    else:
        try:
            va_vec = np.asarray(body_aero.va, dtype=float)
            va = float(np.linalg.norm(va_vec))
        except:
            va = np.nan

    # Extract Cl and Cd from VSM solution
    # Try from results_aero dict first (cl_distribution, cd_distribution)
    cl = np.nan
    cd = np.nan

    try:
        cl_dist = np.asarray(results_aero.get("cl", []), dtype=float).ravel()
        cd_dist = np.asarray(results_aero.get("cd", []), dtype=float).ravel()
        if cl_dist.size > 0 and np.isfinite(cl_dist).any():
            cl = float(np.nanmean(cl_dist))  # Wing-averaged Cl
        if cd_dist.size > 0 and np.isfinite(cd_dist).any():
            cd = float(np.nanmean(cd_dist))  # Wing-averaged Cd
    except Exception as e:
        logging.warning(
            f"Could not extract cl/cd from cl_distribution/cd_distribution: {e}"
        )

    # Fallback: try to extract from body_aero panels after solve
    if np.isnan(cl) or np.isnan(cd):
        try:
            panels = body_aero.panels
            if panels and len(panels) > 0:
                # Try to access cl/cd attributes from panels
                cl_vals = [
                    p.cl for p in panels if hasattr(p, "cl") and np.isfinite(p.cl)
                ]
                cd_vals = [
                    p.cd for p in panels if hasattr(p, "cd") and np.isfinite(p.cd)
                ]
                if cl_vals:
                    cl = float(np.mean(cl_vals))
                if cd_vals:
                    cd = float(np.mean(cd_vals))
                logging.info(f"Extracted Cl={cl:.4f}, Cd={cd:.4f} from panel objects")
        except Exception as e:
            logging.warning(f"Could not extract cl/cd from panel objects: {e}")

    # Also debug: log what's in results_aero
    logging.debug(f"results_aero keys: {list(results_aero.keys())}")

    # Tether force: sum of all forces in Z direction at equilibrium
    try:
        f_aero_total = np.sum(f_aero_wing, axis=0)  # f_aero_wing is shape (n_nodes, 3)
        tether_force = float(
            f_aero_total[2] + gravity_force_total[2] + inertial_force_total[2]
        )
    except:
        tether_force = np.nan

    meta = {
        "total_time_s": time.time() - start_time,
        "n_iter": i + 2,  # +2: 1 for pre-loop initial state + (i+1) loop entries
        "converged": is_convergence,
        "qs_success": bool(results_aero.get("success", False)),
        "opt_x": opt_x,
        "aero_roll_deg": float(results_aero.get("aero_roll_deg", np.nan)),
        "aoa_deg": float(results_aero.get("aoa_deg", np.nan)),
        "side_slip_deg": float(results_aero.get("side_slip_deg", np.nan)),
        "va": va,
        "cl": float(cl),
        "cd": float(cd),
        "tether_force": float(tether_force),
        "rest_lengths": rest_lengths,  # ensure numeric array
        # Convert kite_connectivity to a numeric array for HDF5 compatibility
        "kite_connectivity": np.array(
            [[int(row[0]), int(row[1])] for row in np.array(kite_connectivity_arr)],
            dtype=np.int32,
        ),
    }

    return tracking_data, meta
