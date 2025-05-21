# %% Defing the spring-damper system function
# defining function for the spring-damper system

import time
from tqdm import tqdm
import numpy as np
import pandas as pd
import logging
import matplotlib.pyplot as plt
import os
import sys
from pathlib import Path
import copy
import h5py
from kitesim.utils import save_results
from kitesim.solver import (
    aerodynamic,
    initialisation,
    struc2aero,
    aero2struc,
    structural,
    tracking,
)


def plot_psm(particles, pss_kite_connectivity, f_ext=None, title="PSM State"):
    """
    Plot a ParticleSystem snapshot in 3D.

    Args:
        particles: list of particle objects with attributes
                   .x (np.array shape (3,)), .fixed (bool)
        connectivity_matrix: list of [i, j, ...] giving springs
        f_ext: optional external forces, either:
               – flat array shape (n_pts*3,)
               – array shape (n_pts,3)
        title: figure title
    """
    connectivity_matrix = pss_kite_connectivity
    # unpack positions & fixed mask
    xs = np.array([p.x for p in particles])
    fixed = np.array([p.fixed for p in particles])

    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    # scatter free vs fixed
    ax.scatter(*(xs[~fixed].T), color="blue", label="Free nodes", s=20)
    ax.scatter(*(xs[fixed].T), color="black", label="Fixed nodes", s=20)

    # add the index to each scatter point
    for i, pos in enumerate(xs):
        ax.text(pos[0], pos[1], pos[2], str(i), color="black", fontsize=8)

    # draw springs
    for i, j, *rest in connectivity_matrix:
        p1, p2 = xs[i], xs[j]
        ax.plot(
            [p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], color="gray", linewidth=1
        )

    # optionally draw external forces
    if f_ext is not None:
        arr = np.array(f_ext)
        if arr.ndim == 1:
            arr = arr.reshape(-1, 3)
        for pos, frc in zip(xs, arr):
            ax.quiver(*pos, *frc, length=1, normalize=True, color="red")

    # set aspect ratio to equal
    bb = xs.max(axis=0) - xs.min(axis=0)
    ax.set_box_aspect(bb)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title(title)
    ax.legend()
    plt.show()


def plot_psm_final_vs_initial(initial_positions, final_positions, connectivity_matrix):
    """
    Compare initial vs final PSM states in one 3D plot.

    Args:
        initial_positions: list of [x, v, m, fixed] or objects with .x/.fixed
        final_positions:   same format as initial_positions
        connectivity_matrix: list of [i, j, ...] giving springs
    """

    # unpack coords
    def unpack(pos_list):
        out = []
        for item in pos_list:
            if hasattr(item, "x"):
                out.append(item.x)
            else:
                out.append(item[0])
        return np.array(out)

    xi = unpack(initial_positions)
    xf = unpack(final_positions)

    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    ax.scatter(*(xi.T), color="black", marker="o", s=10, label="Initial")
    ax.scatter(*(xf.T), color="red", marker="o", s=10, label="Final")

    # draw lines
    for i, j, *rest in connectivity_matrix:
        p1i, p2i = xi[i], xi[j]
        p1f, p2f = xf[i], xf[j]
        ax.plot(
            [p1i[0], p2i[0]],
            [p1i[1], p2i[1]],
            [p1i[2], p2i[2]],
            color="black",
            linewidth=1,
        )
        ax.plot(
            [p1f[0], p2f[0]],
            [p1f[1], p2f[1]],
            [p1f[2], p2f[2]],
            color="red",
            linewidth=1,
        )

    # equal aspect
    all_pts = np.vstack((xi, xf))
    bb = all_pts.max(axis=0) - all_pts.min(axis=0)
    ax.set_box_aspect(bb)

    ax.set_title("PSM: Initial (black) vs Final (red)")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.legend()
    plt.show()


def run_aerostructural_solver(config_dict, config_kite_dict, PROJECT_DIR, results_dir):
    """Runs the aero-structural solver for the given input parameters"""

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
        force_gravity = np.array(
            [np.array(config_dict["grav_constant"]) * m_pt for m_pt in m_array]
        )
    else:
        force_gravity = np.zeros(struc_nodes.shape)
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
            f_aero_wing_VSM, body_aero, results_aero = aerodynamic.run_vsm_package(
                body_aero=body_aero,
                solver=vsm_solver,
                le_arr=le_arr,
                te_arr=te_arr,
                va_vector=vel_app,
                aero_input_type="reuse_initial_polar_data",
                initial_polar_data=initial_polar_data,
                is_with_plot=config_dict["is_with_aero_plot_per_iteration"],
            )
            logging.debug(
                f"Aero symmetry check, aero_force_y: { np.sum([force[1] for force in f_aero_wing_VSM])}"
            )

            ### AERO --> STRUC
            f_aero_wing = aero2struc.main(
                config_dict["coupling_method"],
                f_aero_wing_VSM,
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
            f_external = f_aero + force_gravity
            ## rounding the forces to 5 decimal points
            f_external = np.round(f_external, 5)

            # Checking symmetry in the forces
            print(
                f"np.sum(f_aero_wing_VSM): {np.sum(f_aero_wing_VSM[:, 0])}, {np.sum(f_aero_wing_VSM[:, 1])}, {np.sum(f_aero_wing_VSM[:, 2])}"
            )

            print(
                f"np.sum(f_external): {np.sum(f_external[:, 0])}, {np.sum(f_external[:, 1])}, {np.sum(f_external[:, 2])}"
            )

            if config_dict["is_with_plot_per_iteration"]:
                # Compute rest-lengths between the points, using the kite_connectivity
                plot_psm(
                    psystem.particles,
                    pss_kite_connectivity,
                    f_ext=f_external,
                    title=f"i: {i}",
                )

            ### STRUC
            ### f_external is flat, and f_external is 2D, ##TODO: could add this to the name?
            f_external = f_external.flatten()
            end_time_f_ext = time.time()
            begin_time_f_int = time.time()
            psystem = structural.run_pss(psystem, params, f_external)
            # position.loc[step], _ = psystem.x_v_current
            end_time_f_int = time.time()

            # logging.debug(f"position.loc[step].shape: {position.loc[step].shape}")
            logging.debug(f"internal force: {psystem.f_int}")
            logging.debug(f"external force: {f_external}")

            # # TODO: ideally you don't need this here and have aero-also work with the flat format, as this should be faster
            # # saving points in different format
            # points = psystem.x_current_2D
            ## TODO: replacing this function inside the src code to inhere
            # Updating the points
            struc_nodes = np.array([particle.x for particle in psystem.particles])
            f_residual = psystem.f_int + f_external
            f_residual_list.append(np.linalg.norm(np.abs(f_residual)))

            # Update unified tracking dataframe (replaces position update)
            tracking.update_tracking_arrays(
                tracking_data,
                i,
                psystem,
                struc_nodes,
                struc_nodes_prev,
                f_external,
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
        plot_psm_final_vs_initial(
            initial_particles, psystem.particles, pss_kite_connectivity
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
