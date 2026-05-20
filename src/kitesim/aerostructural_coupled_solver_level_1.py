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
    analysis_metrics,
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


def update_power_tape_actuation(
    config,
    psystem,
    kite_fem_structure,
    kite_connectivity_arr,
    power_tape_index,
    power_tape_extension_step,
    initial_length_power_tape,
    power_tape_final_extension,
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
        n_power_tape_steps: Number of power tape extension steps
        rest_lengths: Current rest lengths array (for kite_fem solver)

    Returns:
        tuple: (delta_power_tape, is_actuation_finalized)
            - delta_power_tape: Current change in power tape length
            - is_actuation_finalized: True if actuation is complete, False otherwise
    """
    is_actuation_finalized = True

    if config["structural_solver"] == "pss":
        current_length = float(psystem.extract_rest_length[power_tape_index])
        delta_power_tape = current_length - initial_length_power_tape

        increment, should_update = _compute_power_tape_increment(
            delta_power_tape=delta_power_tape,
            power_tape_final_extension=power_tape_final_extension,
            power_tape_extension_step=power_tape_extension_step,
        )
        if should_update:
            psystem.update_rest_length(power_tape_index, increment)
            current_length = float(psystem.extract_rest_length[power_tape_index])
            delta_power_tape = current_length - initial_length_power_tape
            logging.info(
                f"||--- delta l_d: {delta_power_tape:.3f}m | new l_d: {current_length:.3f}m | Steps required: {n_power_tape_steps}"
            )
            is_actuation_finalized = False

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

        increment, should_update = _compute_power_tape_increment(
            delta_power_tape=delta_power_tape,
            power_tape_final_extension=power_tape_final_extension,
            power_tape_extension_step=power_tape_extension_step,
        )
        if should_update:
            new_length = current_length + increment
            kite_fem_structure.modify_get_spring_rest_length(
                spring_ids=[spring_id],
                new_l0s=[new_length],
            )
            delta_power_tape = new_length - initial_length_power_tape
            logging.info(
                f"||--- delta l_d: {delta_power_tape:.3f}m | new l_d: {new_length:.3f}m | Steps required: {n_power_tape_steps}"
            )
            is_actuation_finalized = False

    return delta_power_tape, is_actuation_finalized


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

    # if residual forces are not changing anymore (compare start of window vs current)
    elif (
        iters_since_start >= n_stag
        and np.abs(f_residual_list[i - n_stag] - f_residual_list[i])
        < config["aero_structural_solver"]["stagnation_tol"]
    ):
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


def _elongation_control_settings(config):
    solver_config = config.get("aero_structural_solver", {})
    return {
        "enabled": bool(
            solver_config.get("is_with_elongation_stiffness_control", True)
        ),
        "limit": float(solver_config.get("max_positive_elongation_ratio", 0.01)),
        "stiffness_factor": float(
            solver_config.get("elongation_stiffness_increase_factor", 1.1)
        ),
        "max_stiffness_factor": float(
            solver_config.get("elongation_stiffness_max_factor", 100.0)
        ),
        "print_report": bool(solver_config.get("is_printing_elongations", True)),
    }


def _get_structural_stiffnesses(config, psystem, kite_fem_structure):
    if config["structural_solver"] == "pss":
        return np.asarray([link.k for link in psystem.springdampers], dtype=float)
    if config["structural_solver"] == "kite_fem":
        return np.asarray(
            [spring_element.k for spring_element in kite_fem_structure.spring_elements],
            dtype=float,
        )
    raise ValueError(f"Unsupported structural solver: {config['structural_solver']}")


def _set_structural_stiffness(config, psystem, kite_fem_structure, idx, stiffness):
    if config["structural_solver"] == "pss":
        psystem.springdampers[idx].k = float(stiffness)
        return
    if config["structural_solver"] == "kite_fem":
        spring_element = kite_fem_structure.spring_elements[idx]
        spring_element.k = float(stiffness)
        spring_element.spring.kxe = float(stiffness)
        return
    raise ValueError(f"Unsupported structural solver: {config['structural_solver']}")


def _safe_elongation_ratio(current_length, rest_length):
    if rest_length is None or rest_length <= 0 or not np.isfinite(rest_length):
        return np.nan
    return (current_length - rest_length) / rest_length


def _pss_effective_element_length(psystem, idx, link):
    current_length = float(link.l)
    rest_length = float(link.l0)
    linktype = str(link.linktype).split(".")[-1].lower()

    if "pulley" in linktype:
        pulley_pairs = getattr(psystem, "_ParticleSystem__pulley_other_line_pair", {})
        other_pair = pulley_pairs.get(str(idx))
        if other_pair is not None:
            idx_p3, idx_p4, rest_len_other = other_pair[:3]
            p3 = psystem.particles[int(idx_p3)].x
            p4 = psystem.particles[int(idx_p4)].x
            current_other = np.linalg.norm(p3 - p4)
            current_length += float(current_other - rest_len_other)

    return current_length, rest_length, linktype


def _element_elongation_rows(
    config,
    psystem,
    kite_fem_structure,
    kite_connectivity_arr,
    initial_stiffnesses,
):
    rows = []
    current_stiffnesses = _get_structural_stiffnesses(
        config, psystem, kite_fem_structure
    )

    if config["structural_solver"] == "pss":
        connectivity = list(kite_connectivity_arr)
        for idx, link in enumerate(psystem.springdampers):
            current_length, rest_length, linktype = _pss_effective_element_length(
                psystem, idx, link
            )
            if idx < len(connectivity):
                node_i, node_j = int(connectivity[idx][0]), int(connectivity[idx][1])
            else:
                node_i, node_j = -1, -1
            rows.append(
                {
                    "idx": idx,
                    "node_i": node_i,
                    "node_j": node_j,
                    "linktype": linktype,
                    "current_length": current_length,
                    "rest_length": rest_length,
                    "elongation_ratio": _safe_elongation_ratio(
                        current_length, rest_length
                    ),
                    "stiffness": float(current_stiffnesses[idx]),
                    "stiffness_factor": (
                        float(current_stiffnesses[idx] / initial_stiffnesses[idx])
                        if initial_stiffnesses[idx] > 0
                        else np.nan
                    ),
                }
            )
        return rows

    if config["structural_solver"] == "kite_fem":
        coords = kite_fem_structure.coords_current
        for idx, spring_element in enumerate(kite_fem_structure.spring_elements):
            _, current_length = spring_element.unit_vector(coords)
            linktype = str(spring_element.springtype).lower()
            rest_length = float(spring_element.l0)
            if linktype == "pulley":
                other = kite_fem_structure.spring_elements[
                    spring_element.i_other_pulley
                ]
                _, other_length = other.unit_vector(coords)
                current_length += other_length

            rows.append(
                {
                    "idx": idx,
                    "node_i": int(spring_element.spring.n1),
                    "node_j": int(spring_element.spring.n2),
                    "linktype": linktype,
                    "current_length": float(current_length),
                    "rest_length": rest_length,
                    "elongation_ratio": _safe_elongation_ratio(
                        float(current_length), rest_length
                    ),
                    "stiffness": float(current_stiffnesses[idx]),
                    "stiffness_factor": (
                        float(current_stiffnesses[idx] / initial_stiffnesses[idx])
                        if initial_stiffnesses[idx] > 0
                        else np.nan
                    ),
                }
            )
        return rows

    raise ValueError(f"Unsupported structural solver: {config['structural_solver']}")


def _print_elongation_report(rows, limit):
    max_positive = max(
        [
            row["elongation_ratio"]
            for row in rows
            if np.isfinite(row["elongation_ratio"])
        ]
        + [-np.inf]
    )
    print("\nElement elongations after residual convergence")
    print(f"Required: positive elongation < {100.0 * limit:.2f}%")
    print(
        f"{'idx':>4} {'nodes':>9} {'type':>16} {'current_l':>11} {'rest_l':>11} "
        f"{'elong':>10} {'k':>12} {'k/k0':>8} {'status':>8}"
    )
    print(
        f"{'-' * 4} {'-' * 9} {'-' * 16} {'-' * 11} {'-' * 11} "
        f"{'-' * 10} {'-' * 12} {'-' * 8} {'-' * 8}"
    )
    for row in rows:
        elong = row["elongation_ratio"]
        status = "FAIL" if np.isfinite(elong) and elong >= limit else "ok"
        elong_str = f"{100.0 * elong:+.3f}%" if np.isfinite(elong) else "nan"
        print(
            f"{row['idx']:>4d} "
            f"{row['node_i']:>3d}-{row['node_j']:<3d} "
            f"{row['linktype'][:16]:>16} "
            f"{row['current_length']:>10.4f}m "
            f"{row['rest_length']:>10.4f}m "
            f"{elong_str:>10} "
            f"{row['stiffness']:>12.3g} "
            f"{row['stiffness_factor']:>8.2f} "
            f"{status:>8}"
        )
    max_positive = max(max_positive, 0.0)
    print(f"Max positive elongation: {100.0 * max_positive:.3f}%\n")


def check_elongation_and_update_stiffness(
    config,
    psystem,
    kite_fem_structure,
    kite_connectivity_arr,
    initial_stiffnesses,
):
    settings = _elongation_control_settings(config)
    rows = _element_elongation_rows(
        config,
        psystem,
        kite_fem_structure,
        kite_connectivity_arr,
        initial_stiffnesses,
    )
    if settings["print_report"]:
        _print_elongation_report(rows, settings["limit"])

    violation_indices = [
        row["idx"]
        for row in rows
        if (
            not np.isfinite(row["elongation_ratio"])
            or row["elongation_ratio"] >= settings["limit"]
        )
    ]
    if not violation_indices:
        return True, False, rows

    current_stiffnesses = _get_structural_stiffnesses(
        config, psystem, kite_fem_structure
    )
    max_allowed = initial_stiffnesses * settings["max_stiffness_factor"]
    if np.all(current_stiffnesses >= max_allowed):
        logging.info(
            "Elongation control non-converged: all stiffnesses reached "
            f"{settings['max_stiffness_factor']:.1f}x their initial values."
        )
        return False, True, rows

    updated = []
    update_indices = [
        idx
        for idx, current_k in enumerate(current_stiffnesses)
        if np.isfinite(current_k)
        and initial_stiffnesses[idx] > 0
        and current_k < max_allowed[idx]
    ]
    for idx in update_indices:
        current_k = current_stiffnesses[idx]
        new_k = min(current_k * settings["stiffness_factor"], max_allowed[idx])
        _set_structural_stiffness(config, psystem, kite_fem_structure, idx, new_k)
        updated.append((idx, current_k, new_k))

    if updated:
        logging.info(
            "Elongation control stiffened all eligible elements: "
            + ", ".join(
                f"{idx} {old_k:.3g}->{new_k:.3g}" for idx, old_k, new_k in updated
            )
        )
        print(
            "Elongation control: increased stiffness by "
            f"{settings['stiffness_factor']:.3g} for all eligible elements "
            f"({len(updated)} updated); "
            "continuing coupled solve."
        )
        return False, False, rows

    logging.info(
        "Elongation control non-converged: overstretched elements remain after "
        "all eligible stiffnesses reached their configured cap. "
        f"Overstretched indices: {violation_indices}"
    )
    return False, True, rows


def increase_stiffness_after_stagnation(
    config,
    psystem,
    kite_fem_structure,
    initial_stiffnesses,
):
    """
    Increase structural stiffness globally after residual stagnation.

    Returns:
        bool: True if at least one element stiffness was increased, False otherwise.
    """

    def format_percent(value):
        if value is None:
            return "n/a"
        return f"{100.0 * value:.3f}%"

    def format_factor(value):
        if value is None:
            return "n/a"
        return f"{value:.2f}x"

    def stagnation_control_diagnostics():
        kite_connectivity_arr = config.get("kite_connectivity_arr", None)

        try:
            rows = _element_elongation_rows(
                config,
                psystem,
                kite_fem_structure,
                kite_connectivity_arr,
                initial_stiffnesses,
            )
        except Exception:
            logging.exception("Could not compute stagnation-control diagnostics.")
            return None

        if not rows:
            return None

        elongations = [
            row["elongation_ratio"]
            for row in rows
            if np.isfinite(row.get("elongation_ratio", np.nan))
        ]

        stiffness_factors = [
            row["stiffness_factor"]
            for row in rows
            if np.isfinite(row.get("stiffness_factor", np.nan))
        ]

        return {
            "max_elongation": max(elongations) if elongations else None,
            "max_stiffness_factor": (
                max(stiffness_factors) if stiffness_factors else None
            ),
        }

    settings = _elongation_control_settings(config)

    stiffness_factor = settings["stiffness_factor"]
    max_stiffness_factor = settings["max_stiffness_factor"]

    current_stiffnesses = _get_structural_stiffnesses(
        config, psystem, kite_fem_structure
    )

    max_allowed = initial_stiffnesses * max_stiffness_factor

    updated = []
    capped = 0
    skipped = 0

    for idx, current_k in enumerate(current_stiffnesses):
        initial_k = initial_stiffnesses[idx]
        max_k = max_allowed[idx]

        if not np.isfinite(current_k) or initial_k <= 0:
            skipped += 1
            continue

        if current_k >= max_k:
            capped += 1
            continue

        new_k = min(current_k * stiffness_factor, max_k)

        _set_structural_stiffness(
            config,
            psystem,
            kite_fem_structure,
            idx,
            new_k,
        )

        updated.append((idx, current_k, new_k, new_k / initial_k))

    if not updated:
        logging.info(
            "Stagnation control: no stiffnesses increased " "(%d capped, %d skipped).",
            capped,
            skipped,
        )
        return False

    min_factor = min(item[3] for item in updated)
    max_factor = max(item[3] for item in updated)

    logging.info(
        "Stagnation control: increased stiffness of %d elements by %.3g "
        "(factor range relative to initial: %.3g–%.3g; %d capped, %d skipped).",
        len(updated),
        stiffness_factor,
        min_factor,
        max_factor,
        capped,
        skipped,
    )

    logging.debug(
        "Updated element stiffnesses: %s",
        ", ".join(
            f"{idx}: {old_k:.3g}->{new_k:.3g} ({factor:.2f}x initial)"
            for idx, old_k, new_k, factor in updated
        ),
    )

    diagnostics = stagnation_control_diagnostics()

    if diagnostics is not None:
        logging.info(
            "Stagnation control diagnostics: max elongation = %s, "
            "max stiffness factor = %s.",
            format_percent(diagnostics["max_elongation"]),
            format_factor(diagnostics["max_stiffness_factor"]),
        )

    return True


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
    print(f'--> Running structural_solver: {config["structural_solver"]}')

    ## PRELOOP
    if config["is_with_gravity"]:
        f_ext_gravity = np.array(
            [np.array(config["grav_constant"]) * m_pt for m_pt in m_arr]
        )
    else:
        f_ext_gravity = np.zeros(struc_nodes.shape)

    if config["structural_solver"] == "kite_fem":
        rest_lengths = kite_fem_structure.modify_get_spring_rest_length()

    max_iter = config["aero_structural_solver"]["max_iter"]
    # Keep index 0 for the pre-loop initial state and reserve max_iter loop slots.
    t_vector = np.linspace(0, max_iter, max_iter + 1)
    tracking_data = tracking.setup_tracking_arrays(len(struc_nodes), t_vector)
    s_ref_m2, s_ref_source = analysis_metrics.resolve_reference_area(
        config,
        struc_nodes,
        struc_node_le_indices,
        struc_node_te_indices,
    )
    is_convergence = False
    f_residual_list = []
    f_tether_drag = np.zeros(3)
    struc_nodes_prev = None  # Initialize previous points for tracking
    start_time = time.time()
    plotting.set_plot_style()
    elongation_settings = _elongation_control_settings(config)
    initial_stiffnesses = _get_structural_stiffnesses(
        config, psystem, kite_fem_structure
    )
    elongation_rows = []
    elongation_control_failed = False
    is_actuation_finalized = False  # becomes True once tape reaches final extension

    # Two-phase regularization: phase 1 = with pseudo_dt, phase 2 = without
    reg_phase = 1  # 1 = regularized, 2 = unregularized (polish)
    stagnation_check_start = 0  # iteration at which current phase started

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

    ### AERO
    f_aero_wing_vsm_format, body_aero, results_aero = aerodynamic_vsm.run_vsm_direct(
        body_aero=body_aero,
        solver=vsm_solver,
        le_arr=le_arr,
        te_arr=te_arr,
        va_vector=vel_app,
        aero_input_type="reuse_initial_polar_data",
        initial_polar_data=initial_polar_data,
        is_with_plot=config["is_with_aero_plot_per_iteration"],
    )
    logging.debug(
        f"Aero symmetry check, f_aero_y: {np.sum([force[1] for force in f_aero_wing_vsm_format])}"
    )
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

    ######################################################################
    # SIMULATION LOOP
    ######################################################################
    ## propagating the simulation for each timestep and saving results
    with tqdm(total=max_iter, desc="Simulating", leave=True) as pbar:
        for i in range(max_iter):
            struc_nodes_before_update = struc_nodes.copy()
            if i > 0:
                struc_nodes_prev = struc_nodes_before_update

            ########################################################
            ############## INTERNAL FORCE CALCULATION ##############
            ########################################################
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

            ### AERO
            f_aero_wing_vsm_format, body_aero, results_aero = (
                aerodynamic_vsm.run_vsm_direct(
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
                )
            else:
                f_aero_bridle = np.zeros((len(struc_nodes), 3))
            f_aero = f_aero_wing + f_aero_bridle

            ## EXTERNAL FORCE
            f_ext = f_aero + f_ext_gravity
            f_ext = np.round(f_ext, 5)
            f_ext_flat = f_ext.flatten()
            end_time_f_ext = time.time()
            aero_metrics = analysis_metrics.compute_global_coefficients(
                f_wing_total=np.sum(f_aero_wing_vsm_format, axis=0),
                f_bridle_total=np.sum(f_aero_bridle, axis=0),
                vel_app=vel_app,
                rho=config["rho"],
                s_ref_m2=s_ref_m2,
            )
            if i == 0:
                aero_force_update_norm = 0.0
            else:
                aero_force_update_norm = np.linalg.norm(
                    aero_metrics["aero_force_total"]
                    - tracking_data["aero_force_total"][i]
                )

            geometry_update_norm = np.linalg.norm(
                struc_nodes - struc_nodes_before_update
            )
            geometry_update_rel_norm = geometry_update_norm / max(
                np.linalg.norm(struc_nodes), 1e-12
            )
            geometry_metrics = analysis_metrics.compute_geometry_metrics(
                struc_nodes,
                struc_node_le_indices,
                struc_node_te_indices,
            )

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
                f_int,
                f_residual_flat=f_residual,
                **aero_metrics,
                **geometry_metrics,
                geometry_update_norm=geometry_update_norm,
                geometry_update_rel_norm=geometry_update_rel_norm,
                aero_force_update_norm=aero_force_update_norm,
                omega_aitken=omega_relaxation,
                regularization_phase=reg_phase,
                is_structural_converged=float(is_structural_converged),
                is_aero_converged=float(results_aero.get("success", True)),
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

            # Two-phase regularization: on stagnation in phase 1, disable
            # pseudo_dt and continue to let the solver polish to true equilibrium.
            if is_stagnated:
                if reg_phase == 1 and config["structural_solver"] == "kite_fem":
                    reg_phase = 2
                    config["structural_kite_fem"]["pseudo_dt"] = None
                    logging.info(
                        f"Phase 1 stagnated at iter {i} (res={np.linalg.norm(f_residual):.1f}N). "
                        f"Switching to phase 2: pseudo_dt=None (no regularization)."
                    )

                # On residual stagnation, increase stiffness globally using the same
                # configured stiffness factor/cap as elongation control.
                stiffness_updated = increase_stiffness_after_stagnation(
                    config=config,
                    psystem=psystem,
                    kite_fem_structure=kite_fem_structure,
                    initial_stiffnesses=initial_stiffnesses,
                )

                if stiffness_updated:
                    stagnation_check_start = i  # reset stagnation window
                    r_prev_flat = None
                    continue

                logging.info(
                    "Classic PS non-converging - residual no longer changes and all stiffnesses are at configured cap"
                )
                elongation_control_failed = True
                should_break = True

            ### ELONGATION CHECK (only when converged)
            if is_convergence:
                if elongation_settings["enabled"]:
                    (
                        is_elongation_ok,
                        elongation_should_break,
                        elongation_rows,
                    ) = check_elongation_and_update_stiffness(
                        config=config,
                        psystem=psystem,
                        kite_fem_structure=kite_fem_structure,
                        kite_connectivity_arr=kite_connectivity_arr,
                        initial_stiffnesses=initial_stiffnesses,
                    )

                    if not is_elongation_ok:
                        is_convergence = False
                        r_prev_flat = None
                        stagnation_check_start = i
                        if elongation_should_break:
                            elongation_control_failed = True
                            should_break = True
                            break
                        else:
                            continue

            ### ACTUATION — applied every iteration until tape reaches final extension
            if not is_actuation_finalized:
                _, is_actuation_finalized = update_power_tape_actuation(
                    config=config,
                    psystem=psystem,
                    kite_fem_structure=kite_fem_structure,
                    kite_connectivity_arr=kite_connectivity_arr,
                    power_tape_index=power_tape_index,
                    power_tape_extension_step=power_tape_extension_step,
                    initial_length_power_tape=initial_length_power_tape,
                    power_tape_final_extension=power_tape_final_extension,
                    n_power_tape_steps=n_power_tape_steps,
                    rest_lengths=(
                        rest_lengths
                        if config["structural_solver"] == "kite_fem"
                        else None
                    ),
                )

            # Check if we should exit the loop
            if should_break or is_convergence:
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
    vec_wind = vel_app / np.linalg.norm(vel_app)

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

    print(f"alpha = {np.degrees(angle):.2f}° (va vs mid-span chord)")
    alpha_at_ac_values = np.ravel(results_aero.get("alpha_at_ac", []))
    if alpha_at_ac_values.size > mid_idx:
        print(
            f'alpha = {float(np.rad2deg(alpha_at_ac_values[mid_idx])):.2f}° (incl. induced velocity, from results_aero["alpha_at_ac"])'
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
            rest_lengths=rest_lengths,
            struc_nodes_initial=struc_nodes_initial,
            title="Initial vs final",
            pulley_line_indices=pulley_line_indices,
            pulley_line_to_other_node_pair_dict=pulley_line_to_other_node_pair_dict,
            vel_app=vel_app,
        )
    final_idx = min(i + 1, len(tracking_data["positions"]) - 1)
    final_geom_metrics = analysis_metrics.compute_geometry_metrics(
        struc_nodes,
        struc_node_le_indices,
        struc_node_te_indices,
    )
    final_depower_tape_length_m = (
        float(rest_lengths[power_tape_index])
        if power_tape_index is not None and len(rest_lengths) > power_tape_index
        else np.nan
    )
    final_power_tape_extension_m = (
        final_depower_tape_length_m - float(initial_length_power_tape)
        if np.isfinite(final_depower_tape_length_m)
        and initial_length_power_tape is not None
        else np.nan
    )
    final_alpha_at_ac = np.asarray(results_aero.get("alpha_at_ac", []), dtype=float)
    final_midspan_alpha_deg = (
        float(np.rad2deg(final_alpha_at_ac.reshape(-1)[mid_idx]))
        if final_alpha_at_ac.size > mid_idx
        else np.nan
    )
    try:
        final_elongation_rows = _element_elongation_rows(
            config,
            psystem,
            kite_fem_structure,
            kite_connectivity_arr,
            initial_stiffnesses,
        )
    except Exception:
        final_elongation_rows = elongation_rows
    final_elongation_ratios = np.asarray(
        [row["elongation_ratio"] for row in final_elongation_rows], dtype=float
    )
    finite_final_elongations = final_elongation_ratios[
        np.isfinite(final_elongation_ratios)
    ]
    final_max_positive_elongation_ratio = (
        float(max(np.max(finite_final_elongations), 0.0))
        if finite_final_elongations.size
        else np.nan
    )
    final_stiffnesses = _get_structural_stiffnesses(config, psystem, kite_fem_structure)
    final_stiffness_factors = np.divide(
        final_stiffnesses,
        initial_stiffnesses,
        out=np.full_like(final_stiffnesses, np.nan, dtype=float),
        where=initial_stiffnesses > 0,
    )

    meta = {
        "total_time_s": time.time() - start_time,
        "n_iter": i + 2,  # +2: 1 for pre-loop initial state + (i+1) loop entries
        "converged": is_convergence,
        "final_residual_norm": float(tracking_data["residual_norm"][final_idx]),
        "final_geometry_update_norm": float(
            tracking_data["geometry_update_norm"][final_idx]
        ),
        "final_C_L_wing": float(tracking_data["C_L_wing"][final_idx]),
        "final_C_D_wing": float(tracking_data["C_D_wing"][final_idx]),
        "final_glide_ratio_wing": float(tracking_data["glide_ratio_wing"][final_idx]),
        "final_C_L_total_aero": float(tracking_data["C_L_total_aero"][final_idx]),
        "final_C_D_total_aero": float(tracking_data["C_D_total_aero"][final_idx]),
        "final_glide_ratio_total_aero": float(
            tracking_data["glide_ratio_total_aero"][final_idx]
        ),
        "final_V_a": float(np.linalg.norm(vel_app)),
        "S_ref_m2": float(s_ref_m2),
        "S_ref_source": str(s_ref_source),
        "rho": float(config["rho"]),
        "coefficient_force_source": "wing_only_panel_forces_and_optional_bridle_total",
        "wind_axis_convention": "drag_positive_along_apparent_wind_lift_positive_z_component",
        "final_depower_tape_length_m": final_depower_tape_length_m,
        "final_u_dp": (
            (final_depower_tape_length_m - 0.2) / 5.0
            if np.isfinite(final_depower_tape_length_m)
            else np.nan
        ),
        "final_power_tape_extension_m": final_power_tape_extension_m,
        "final_pitch_or_trim_deg": float(np.degrees(angle)),
        "final_center_or_midspan_alpha_deg": final_midspan_alpha_deg,
        "final_projected_span_m": float(final_geom_metrics["projected_span_m"]),
        "final_projected_area_m2": float(final_geom_metrics["projected_area_m2"]),
        "final_mean_twist_deg": float(final_geom_metrics["mean_twist_deg"]),
        "final_max_abs_twist_deg": float(final_geom_metrics["max_abs_twist_deg"]),
        "elongation_control_enabled": bool(elongation_settings["enabled"]),
        "elongation_control_failed": bool(elongation_control_failed),
        "max_positive_elongation_limit_ratio": float(elongation_settings["limit"]),
        "final_max_positive_elongation_ratio": final_max_positive_elongation_ratio,
        "final_elongation_all_below_limit": bool(
            np.isfinite(final_max_positive_elongation_ratio)
            and final_max_positive_elongation_ratio < elongation_settings["limit"]
        ),
        "structural_stiffnesses_initial": initial_stiffnesses,
        "structural_stiffnesses_final": final_stiffnesses,
        "structural_stiffness_factors_final": final_stiffness_factors,
        "rest_lengths": rest_lengths,  # ensure numeric array
        # Convert kite_connectivity to a numeric array for HDF5 compatibility
        "kite_connectivity": np.array(
            [[int(row[0]), int(row[1])] for row in np.array(kite_connectivity_arr)],
            dtype=np.int32,
        ),
        "struc_node_le_indices": np.asarray(struc_node_le_indices, dtype=np.int32),
        "struc_node_te_indices": np.asarray(struc_node_te_indices, dtype=np.int32),
    }

    return tracking_data, meta
