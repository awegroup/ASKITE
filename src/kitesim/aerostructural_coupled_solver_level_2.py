import time
from tqdm import tqdm
import numpy as np
import logging
from pathlib import Path
import copy
from kitesim import (
    aerodynamic_vsm,
    struc2aero,
    aero2struc_level_2,
    structural_pss,
    tracking,
    plotting,
    structural_kite_fem_level_2,
    aerodynamic_bridle_line_drag,
)


def build_symmetry_mapping(struc_nodes, tol=1e-5):
    """Build a mapping of symmetrical node pairs from the *initial* geometry.

    For every node with y > 0 (positive-span side), find its mirror partner
    at (x, −y, z) on the negative-span side.  Also identify nodes that sit
    on the symmetry plane (|y| ≤ tol).

    Call this **once** during initialisation and pass the result to
    :func:`forcing_symmetry` at every iteration.

    Parameters
    ----------
    struc_nodes : np.ndarray, shape (n, 3)
        Initial (undeformed) node coordinates.
    tol : float
        Absolute tolerance for matching mirror coordinates.

    Returns
    -------
    symmetry_mapping : dict
        ``{"pairs": np.ndarray shape (m, 2),
           "center_indices": list[int]}``
        *pairs[:, 0]* = positive-y node index (source),
        *pairs[:, 1]* = negative-y node index (mirror).
    """
    pos_indices = [i for i, pt in enumerate(struc_nodes) if pt[1] > tol]
    neg_indices = [i for i, pt in enumerate(struc_nodes) if pt[1] < -tol]
    center_indices = [i for i, pt in enumerate(struc_nodes) if abs(pt[1]) <= tol]

    pairs = []
    for pi in pos_indices:
        mirrored = np.array(
            [struc_nodes[pi][0], -struc_nodes[pi][1], struc_nodes[pi][2]]
        )
        for ni in neg_indices:
            if np.allclose(struc_nodes[ni], mirrored, atol=tol):
                pairs.append((pi, ni))
                break

    pairs = np.array(pairs) if pairs else np.empty((0, 2), dtype=int)

    logging.info(
        f"Symmetry mapping: {len(pairs)} pairs, {len(center_indices)} center nodes"
    )
    return {"pairs": pairs, "center_indices": center_indices}


def forcing_symmetry(struc_nodes, symmetry_mapping):
    """Force y-symmetry on the structural nodes using a pre-built mapping.

    For each (source, mirror) pair the mirror node is set to
    ``[x_source, -y_source, z_source]``.  Centre-plane nodes are forced to
    ``y = 0``.

    Parameters
    ----------
    struc_nodes : np.ndarray, shape (n, 3)
    symmetry_mapping : dict
        Output of :func:`build_symmetry_mapping`.

    Returns
    -------
    struc_nodes : np.ndarray
        The (modified in-place) node array.
    """
    for src, mir in symmetry_mapping["pairs"]:
        struc_nodes[mir] = np.array(
            [struc_nodes[src][0], -struc_nodes[src][1], struc_nodes[src][2]]
        )
    for ci in symmetry_mapping["center_indices"]:
        struc_nodes[ci][1] = 0.0
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
        kite_connectivity_arr: ASKITE connectivity array
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
        current_length = float(psystem.extract_rest_length[power_tape_index])
        delta_power_tape = current_length - initial_length_power_tape

        if is_residual_below_tol:
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

        if is_residual_below_tol:
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

    n_stag = config["aero_structural_solver"].get("n_max_constant_residual_force", 15)
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
    elif iters_since_start > n_stag and np.abs(
        f_residual_list[i - n_stag] - f_residual_list[i]
    ) < config["aero_structural_solver"].get("stagnation_tol", 1.0):
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
    canopy_sections=None,
    strut_sections=None,
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
    is_convergence = False
    f_residual_list = []
    f_tether_drag = np.zeros(3)
    is_residual_below_tol = False
    struc_nodes_prev = None  # Initialize previous points for tracking
    start_time = time.time()
    plotting.set_plot_style()

    # Two-phase regularization: phase 1 = with pseudo_dt, phase 2 = without
    reg_phase = 1  # 1 = regularized, 2 = unregularized (polish)
    stagnation_check_start = 0  # iteration at which current phase started

    # Aitken relaxation state
    omega_relaxation = config["aero_structural_solver"].get("relaxation_factor", 0.3)
    r_prev_flat = None

    # Build symmetry mapping once from initial (undeformed) geometry
    if config["is_with_forcing_symmetry"]:
        symmetry_mapping = build_symmetry_mapping(struc_nodes_initial)

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
    f_aero_wing_vsm_format, body_aero, results_aero = aerodynamic_vsm.run_vsm_package(
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

    # # TODO: debuggin here
    # # print out the cp locations as % of le to te for debugging purposes
    # alpha_arr = np.array(results_aero["alpha_at_ac"]).ravel()
    # alpha_arr_geom = np.array(results_aero["alpha_geometric"]).ravel()
    # cl_arr = np.array(results_aero["cl_distribution"]).ravel()
    # for i, (panel, alpha, alpha_geom, cl) in enumerate(
    #     zip(body_aero.panels, alpha_arr, alpha_arr_geom, cl_arr)
    # ):
    #     cp = np.array(results_aero["panel_cp_locations"][i])
    #     le_mid = 0.5 * (panel.LE_point_1 + panel.LE_point_2)
    #     te_mid = 0.5 * (panel.TE_point_1 + panel.TE_point_2)
    #     # chordwise fraction along panel chord axis
    #     cp_rel = np.dot(cp - le_mid, panel.y_airf) / panel.chord
    #     print(
    #         f"i:{i}, CP: {cp_rel:.3f}, alpha_corr: {alpha:.2f}deg, alpha_geom: {alpha_geom:.2f}deg, cl: {cl:.3f}, le: {le_mid}, te: {te_mid}"
    #     )
    #     F = np.array(results_aero["F_distribution"][i])
    #     M = np.array(results_aero["M_distribution"][i])
    #     ac = panel.aerodynamic_center
    #     y_airf = panel.y_airf
    #     z_airf = panel.z_airf
    #     c = panel.chord

    #     r = ac  # reference_point is [0,0,0] in config
    #     M_local = M - np.cross(r, F)
    #     m_pitch = np.dot(M_local, z_airf)
    #     F_perp = F - np.dot(F, z_airf) * z_airf
    #     F_perp_mag = np.linalg.norm(F_perp)
    #     lever_raw = m_pitch / max(F_perp_mag, 1e-12)
    #     lever_clamped = np.clip(lever_raw, -0.25 * c, 0.75 * c)
    #     cp_rel = 0.25 + lever_clamped / c

    #     print(
    #         f"i:{i}, cp_rel:{cp_rel:.3f}, lever_raw:{lever_raw/c:.3f}, F_perp:{F_perp_mag:.3e}"
    #     )

    #     cd, cm = panel.compute_cd_cm(alpha)
    #     print(
    #         f"i:{i}, cm:{cm:.4f}, alpha_corr:{alpha:.2f}, alpha_geom:{alpha_geom:.2f}, cp:{cp_rel:.3f}"
    #     )

    ### AERO --> STRUC
    f_aero_wing = aero2struc_level_2.main(
        config["aero2struc"]["coupling_method"],
        f_aero_wing_vsm_format,
        struc_nodes,
        np.array(results_aero["panel_cp_locations"]),
        aero2struc_mapping,
        config["is_with_coupling_plot_per_iteration"],
        config["aero2struc"],
        canopy_sections,
        strut_sections,
        body_aero.panels,
    )

    # Check moment preservation of aero→struc mapping (pre-loop)
    aero2struc_level_2.check_moment_preservation(
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
            if i > 0:
                struc_nodes_prev = struc_nodes.copy()

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
                    structural_kite_fem_level_2.run_kite_fem(
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
                    rest_lengths = structural_kite_fem_level_2.get_rest_lengths(
                        kite_fem_structure, kite_connectivity_arr
                    )
                    # kite_fem_structure.plot_convergence()  # not available in kite_fem

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
            f_aero_wing = aero2struc_level_2.main(
                config["aero2struc"]["coupling_method"],
                f_aero_wing_vsm_format,
                struc_nodes,
                np.array(results_aero["panel_cp_locations"]),
                aero2struc_mapping,
                config["is_with_coupling_plot_per_iteration"],
                config["aero2struc"],
                canopy_sections,
                strut_sections,
                body_aero.panels,
            )

            # Check moment preservation (only first coupling iteration to limit log spam)
            if i == 1:
                aero2struc_level_2.check_moment_preservation(
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

            ### FORCING SYMMETRY
            if config["is_with_forcing_symmetry"]:
                logging.info("Forcing symmetry in y-direction")
                struc_nodes = forcing_symmetry(struc_nodes, symmetry_mapping)

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

            # Two-phase regularization: on stagnation in phase 1, disable
            # pseudo_dt and continue to let the solver polish to true equilibrium.
            if is_stagnated:
                if reg_phase == 1 and config["structural_solver"] == "kite_fem":
                    reg_phase = 2
                    stagnation_check_start = i  # reset stagnation window
                    config["structural_kite_fem"]["pseudo_dt"] = None
                    logging.info(
                        f"Phase 1 stagnated at iter {i} "
                        f"(res={np.linalg.norm(f_residual):.1f}N). "
                        f"Switching to phase 2: pseudo_dt=None (no regularization)."
                    )
                else:
                    logging.info(
                        "Classic PS non-converging - residual no longer changes"
                    )
                    should_break = True

            ### ACTUATION (only when converged)
            if is_convergence:
                # Update residual flag for actuation function
                is_residual_below_tol = is_convergence

                delta_power_tape, is_actuation_finalized = update_power_tape_actuation(
                    config=config,
                    psystem=psystem,
                    kite_fem_structure=kite_fem_structure,
                    kite_connectivity_arr=kite_connectivity_arr,
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
    print(
        f'alpha = {float(np.rad2deg(results_aero["alpha_at_ac"][mid_idx])):.2f}° (incl. induced velocity, from results_aero["alpha_at_ac"])'
    )

    if config["structural_solver"] == "pss":
        rest_lengths = psystem.extract_rest_length
    elif config["structural_solver"] == "kite_fem":
        rest_lengths = structural_kite_fem_level_2.get_rest_lengths(
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
