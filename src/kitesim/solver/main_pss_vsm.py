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
from kitesim.coupling import coupling_struc2aero, coupling_aero2struc
from kitesim.solver.vsm_functions import initialize_vsm, run_vsm_package
from kitesim.solver.pss_functions import instantiate_psystem, run_pss
from kitesim.solver.initialisation import initialising_solver


def setup_tracking_arrays(n_pts, t_vector):
    nt = len(t_vector)
    return {
        "positions": np.zeros((nt, n_pts, 3)),
        "f_external": np.zeros((nt, n_pts, 3)),
        "f_residual": np.zeros((nt, n_pts, 3)),
        "residual_norm": np.zeros(nt),
        "max_residual": np.zeros(nt),
        "pos_change": np.zeros(nt),
        "vel_change": np.zeros(nt),
    }


def update_tracking_arrays(
    tracking, step, psystem, points, points_prev, f_external_flat, f_residual_flat
):
    # Unpack 3D storage
    pos3d = tracking["positions"]
    ext3d = tracking["f_external"]
    res3d = tracking["f_residual"]

    n_pts = pos3d.shape[1]

    # 1) Positions
    pos_flat, _ = psystem.x_v_current  # shape (n_pts*3,)
    pos3d[step] = pos_flat.reshape(n_pts, 3)

    # 2) External & residual forces: reshape before storing
    ext3d[step] = f_external_flat.reshape(n_pts, 3)
    res3d[step] = f_residual_flat.reshape(n_pts, 3)

    # 3) Norms
    tracking["residual_norm"][step] = np.linalg.norm(f_residual_flat)
    tracking["max_residual"][step] = np.max(np.abs(f_residual_flat))

    # 4) Deltas
    if points_prev is not None:
        dp = points - points_prev
        tracking["pos_change"][step] = np.linalg.norm(dp)

        if hasattr(psystem, "v_current_2D") and hasattr(psystem, "v_prev_2D"):
            dv = psystem.v_current_2D - psystem.v_prev_2D
            tracking["vel_change"][step] = np.linalg.norm(dv)


def save_results(tracking, meta, filename):
    with h5py.File(filename, "w") as f:
        grp = f.create_group("tracking")
        for name, arr in tracking.items():
            grp.create_dataset(name, data=arr[: meta["n_iter"]], compression="gzip")
        for k, v in meta.items():
            grp.attrs[k] = v


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
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_box_aspect(bb)
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

    ### SOLVER INITIALISATION
    (
        points,
        le_idx,
        te_idx,
        wing_connectivity,
        bridle_connectivity,
        kite_connectivity,
        wing_rest_lengths_initial,
        rest_lengths,
        m_array,
    ) = initialising_solver(config_kite_dict)

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
    wing = body_aero.wings[0]
    new_sections = wing.refine_aerodynamic_mesh()
    initial_polar_data = []
    for new_section in new_sections:
        initial_polar_data.append(new_section.aero_input)
    # print(f"new_sections: {new_sections}")
    # print(f"aero_input_arr: {len(aero_input_arr)}")
    # breakpoint()

    ## STRUC initialisation -- PSS
    psystem, params, pss_kite_connectivity = instantiate_psystem(
        config_dict,
        config_kite_dict,
        points,
        wing_connectivity,
        kite_connectivity,
        rest_lengths,
        m_array,
    )
    initial_particles = copy.deepcopy(psystem.particles)
    if config_dict["is_with_initial_plot"]:
        plot_psm(
            initial_particles, pss_kite_connectivity, f_ext=None, title="PSM State"
        )

    # INITIALISATION OF VARIABLES
    aero_structural_tol = params["aerostructural_tol"]
    is_run_only_1_time_step = config_dict["is_run_only_1_time_step"]
    coupling_method = config_dict["coupling_method"]
    n_chordwise_aero_nodes = config_dict["aero"]["n_chordwise_aero_nodes"]
    aero_structural_max_iter = config_dict["aero_structural"]["max_iter"]
    # actuation
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

    if config_dict["is_with_gravity"]:
        force_gravity = np.array(
            [np.array(config_dict["grav_constant"]) * m_pt for m_pt in m_array]
        )
    else:
        force_gravity = np.zeros(points.shape)

    ## coupling initialisation
    (
        points_wing_segment_corners_aero_orderded,
        index_transformation_struc_to_aero,
    ) = coupling_struc2aero.extract_wingpanel_corners_L_to_R(
        points, np.array(config_kite_dict["plate_point_indices"], dtype=int)
    )
    logging.debug(
        f"points_wing_segment_corners_aero_orderded: {points_wing_segment_corners_aero_orderded}"
    )
    ### creating a dict with key value pairs, to transform from aero to struc
    index_transformation_aero_to_struc_dict = {}
    for i, value in enumerate(index_transformation_struc_to_aero):
        index_transformation_aero_to_struc_dict[value] = i

    # PreLoop Initialisation
    t_vector = np.linspace(
        params["dt"], params["t_steps"] * params["dt"], params["t_steps"]
    )
    tracking = setup_tracking_arrays(len(points), t_vector)
    vel_app = np.array(config_dict["vel_wind"]) - np.array(config_dict["vel_kite"])
    is_convergence = False
    f_residual_list = []
    f_tether_drag = np.zeros(3)
    is_residual_below_tol = False
    points_prev = None  # Initialize previous points for tracking
    start_time = time.time()

    # le_arr = points_wing_segment_corners_aero_orderded[0::2, :]
    # te_arr = points_wing_segment_corners_aero_orderded[1::2, :]
    # print(f"le_arr: len(le_arr): {len(le_arr)}")
    # print(f"te_arr: len(te_arr): {len(te_arr)}")
    # f_aero_wing_VSM, body_aero = run_vsm_package(
    #     body_aero=body_aero,
    #     solver=vsm_solver,
    #     le_arr=points_wing_segment_corners_aero_orderded[0::2, :],
    #     te_arr=points_wing_segment_corners_aero_orderded[1::2, :],
    #     va_vector=vel_app,
    #     aero_input_type="reuse_initial_polar_data",
    #     initial_polar_data=initial_polar_data,
    # )
    # print(f"f_aero_wing_VSM: {f_aero_wing_VSM}")
    # breakpoint()

    ######################################################################
    # SIMULATION LOOP
    ######################################################################
    ## propagating the simulation for each timestep and saving results
    with tqdm(total=len(t_vector), desc="Simulating", leave=True) as pbar:
        for i, step in enumerate(t_vector):
            if i > 0:
                points_prev = points.copy()

            ## external force
            begin_time_f_ext = time.time()
            ### STRUC --> AERO
            # Ordering the points from left to right (-y to +y), the desired aero-ordering
            points_wing_segment_corners_aero_orderded = points[
                index_transformation_struc_to_aero
            ]
            ### AERO
            if i == 0:
                f_aero_wing_VSM, body_aero = run_vsm_package(
                    body_aero=body_aero,
                    solver=vsm_solver,
                    le_arr=points_wing_segment_corners_aero_orderded[0::2, :],
                    te_arr=points_wing_segment_corners_aero_orderded[1::2, :],
                    va_vector=vel_app,
                    aero_input_type="reuse_initial_polar_data",
                    initial_polar_data=initial_polar_data,
                )
                ### AERO --> STRUC
                if coupling_method == "NN":
                    f_aero_wing = coupling_aero2struc.aero2struc_NN_vsm(
                        n_chordwise_aero_nodes,
                        body_aero,
                        f_aero_wing_VSM,
                        points_wing_segment_corners_aero_orderded,
                        points,
                        config_kite_dict["plate_point_indices"],
                        config_dict["is_with_coupling_plot"],
                    )
                else:
                    raise ValueError(
                        "Coupling method not recognized; wrong name or typo"
                    )

                # TODO: get bridle line forces back in to play
                f_aero_bridle = np.zeros(points.shape)
                f_aero = f_aero_wing + f_aero_bridle

                ## summing up
                f_external = f_aero + force_gravity
                ## rounding the forces to 5 decimal points
                f_external = np.round(f_external, 5)

                # Checking symmetry in the forces
                print(
                    f"np.sum(f_external): {np.sum(f_external[:, 0])}, {np.sum(f_external[:, 1])}, {np.sum(f_external[:, 2])}"
                )

            # Compute rest-lengths between the points, using the kite_connectivity
            plot_psm(
                psystem.particles,
                pss_kite_connectivity,
                f_ext=f_external,
                title=f"i: {i}",
            )

            # def compute_delta_rest_lengths(psystem):
            #     rest_lengths = []
            #     delta_rest_lengths = []
            #     for i, element in enumerate(psystem.springdampers):
            #         delta = element.l - element.l0
            #         rest_lengths.append(element.l0)
            #         delta_rest_lengths.append(delta)
            #         if delta > 0:
            #             print(f"i: {i} | delta: {delta} | l0: {element.l0}")
            #     # print(f"delta_rest_lengths: {delta_rest_lengths}")

            # compute_delta_rest_lengths(psystem)

            ### STRUC
            ### f_external is flat, and f_external is 2D, ##TODO: could add this to the name?
            f_external = f_external.flatten()
            end_time_f_ext = time.time()
            begin_time_f_int = time.time()
            psystem = run_pss(psystem, params, f_external)
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
            points = np.array([particle.x for particle in psystem.particles])
            f_residual = psystem.f_int + f_external
            f_residual_list.append(np.linalg.norm(np.abs(f_residual)))

            # Update unified tracking dataframe (replaces position update)
            update_tracking_arrays(
                tracking,
                i,
                psystem,
                points,
                points_prev,
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
                and np.linalg.norm(f_residual) <= aero_structural_tol
            ):
                is_residual_below_tol = True
                is_convergence = True
            # if np.linalg.norm(f_residual) <= aero_structural_tol and i > 50:
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
            elif i > aero_structural_max_iter:
                is_convergence = False
                logging.info(
                    f"Classic PS non-converging - more than max ({aero_structural_max_iter}) iterations needed"
                )
                break
            # special case for running the simulation for only one timestep
            elif is_run_only_1_time_step:
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
    meta = {
        "total_time_s": time.time() - start_time,
        "n_iter": i,
        "converged": is_convergence,
        "wing_rest_lengths": psystem.extract_rest_length[
            :len_wing_rest_length
        ].tolist(),
        "bridle_rest_lengths": psystem.extract_rest_length[
            len_wing_rest_length:
        ].tolist(),
    }
    h5_path = Path(results_dir) / "sim_output.h5"
    save_results(tracking, meta, h5_path)

    if config_dict["is_with_final_plot"]:
        plot_psm_final_vs_initial(
            initial_particles, psystem.particles, pss_kite_connectivity
        )

    return {"hdf5_path": str(h5_path)}
