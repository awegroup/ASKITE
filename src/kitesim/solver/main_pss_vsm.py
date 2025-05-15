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
from kitesim.solver.vsm_functions import initialize_vsm, run_vsm_package
from kitesim.solver.pss_functions import instantiate_psystem, run_pss
from kitesim.solver.initialisation import initialising_solver


# Set up unified tracking dataframe including position, forces, and residuals
def setup_tracking_dataframe(params, t_vector):
    """
    Initialize unified dataframe to track position, forces, and residuals across iterations.

    Args:
        params: Dictionary with simulation parameters
        t_vector: Time vector for simulation steps

    Returns:
        Unified tracking dataframe
    """
    n_points = params["n"]
    columns = {}

    # Position columns (per node)
    for i in range(n_points):
        columns[f"x{i+1}"] = np.zeros(len(t_vector))
        columns[f"y{i+1}"] = np.zeros(len(t_vector))
        columns[f"z{i+1}"] = np.zeros(len(t_vector))

    # Force columns (per node)
    for i in range(n_points):
        # Aerodynamic forces
        columns[f"aero_x{i+1}"] = np.zeros(len(t_vector))
        columns[f"aero_y{i+1}"] = np.zeros(len(t_vector))
        columns[f"aero_z{i+1}"] = np.zeros(len(t_vector))

        # Internal forces
        columns[f"int_x{i+1}"] = np.zeros(len(t_vector))
        columns[f"int_y{i+1}"] = np.zeros(len(t_vector))
        columns[f"int_z{i+1}"] = np.zeros(len(t_vector))

        # External forces (sum of all external forces)
        columns[f"ext_x{i+1}"] = np.zeros(len(t_vector))
        columns[f"ext_y{i+1}"] = np.zeros(len(t_vector))
        columns[f"ext_z{i+1}"] = np.zeros(len(t_vector))

    # Residual and convergence columns
    columns["residual_norm"] = np.zeros(len(t_vector))  # Overall residual norm
    columns["max_residual"] = np.zeros(len(t_vector))  # Maximum residual component
    columns["pos_change_norm"] = np.zeros(
        len(t_vector)
    )  # Position change between iterations
    columns["vel_change_norm"] = np.zeros(
        len(t_vector)
    )  # Velocity change between iterations

    # Create unified dataframe
    tracking_df = pd.DataFrame(index=t_vector, columns=columns)

    return tracking_df


# Update the unified tracking dataframe during simulation
def update_tracking(
    tracking_df,
    step,
    psystem,
    points,
    points_prev,
    force_aero,
    force_external,
    residual_f,
):
    """
    Update tracking dataframe with current iteration data

    Args:
        tracking_df: Unified tracking dataframe
        step: Current time step
        psystem: Particle system object
        points: Current position of points
        points_prev: Position of points from previous iteration
        force_aero: Aerodynamic forces
        force_external: Total external forces
        residual_f: Residual forces
    """
    n_points = len(points)

    # Reshape forces to match points shape if needed
    f_int_2d = psystem.f_int.reshape(-1, 3)

    # Update position data (replace existing implementation)
    pos_data, _ = psystem.x_v_current
    for i in range(n_points):
        tracking_df.loc[step, f"x{i+1}"] = pos_data[i * 3]
        tracking_df.loc[step, f"y{i+1}"] = pos_data[i * 3 + 1]
        tracking_df.loc[step, f"z{i+1}"] = pos_data[i * 3 + 2]

    # Update force data
    for i in range(n_points):
        # Aerodynamic forces
        tracking_df.loc[step, f"aero_x{i+1}"] = force_aero[i, 0]
        tracking_df.loc[step, f"aero_y{i+1}"] = force_aero[i, 1]
        tracking_df.loc[step, f"aero_z{i+1}"] = force_aero[i, 2]

        # Internal forces
        tracking_df.loc[step, f"int_x{i+1}"] = f_int_2d[i, 0]
        tracking_df.loc[step, f"int_y{i+1}"] = f_int_2d[i, 1]
        tracking_df.loc[step, f"int_z{i+1}"] = f_int_2d[i, 2]

        # External forces
        tracking_df.loc[step, f"ext_x{i+1}"] = force_external[i, 0]
        tracking_df.loc[step, f"ext_y{i+1}"] = force_external[i, 1]
        tracking_df.loc[step, f"ext_z{i+1}"] = force_external[i, 2]

    # Update residual data
    residual_norm = np.linalg.norm(residual_f)
    max_residual = np.max(np.abs(residual_f))

    # Calculate position and velocity changes if previous points are available
    if points_prev is not None:
        pos_change = points - points_prev
        pos_change_norm = np.linalg.norm(pos_change)

        # Velocity change (if velocity is tracked)
        if hasattr(psystem, "v_current_2D") and hasattr(psystem, "v_prev_2D"):
            vel_change = psystem.v_current_2D - psystem.v_prev_2D
            vel_change_norm = np.linalg.norm(vel_change)
        else:
            vel_change_norm = 0.0
    else:
        pos_change_norm = 0.0
        vel_change_norm = 0.0

    tracking_df.loc[step, "residual_norm"] = residual_norm
    tracking_df.loc[step, "max_residual"] = max_residual
    tracking_df.loc[step, "pos_change_norm"] = pos_change_norm
    tracking_df.loc[step, "vel_change_norm"] = vel_change_norm


def run_aerostructural_solver(config_dict, config_kite_dict, PROJECT_DIR, results_dir):
    """Runs the aero-structural solver for the given input parameters"""

    ### SOLVER INITIALISATION
    (
        points,
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

    ## STRUC initialisation -- PSS
    psystem, params = instantiate_psystem(
        config_dict,
        config_kite_dict,
        points,
        wing_connectivity,
        kite_connectivity,
        rest_lengths,
        m_array,
    )

    ##TODO: add a function to set the initial position of the kite
    # if config_dict["is_with_initial_plot"]:
    #     plot the geometry

    ## setting up the position-dataframe
    t_vector = np.linspace(
        params["dt"], params["t_steps"] * params["dt"], params["t_steps"]
    )
    tracking_df = setup_tracking_dataframe(params, t_vector)

    # INITIALISATION OF VARIABLES
    ## parameters used in the loop
    points_prev = None  # Initialize previous points for tracking
    aero_structural_tol = params["aerostructural_tol"]
    vel_app = np.array(config_dict["vel_wind"]) - np.array(config_dict["vel_kite"])
    start_time = time.time()
    is_convergence = False
    residual_f_list = []
    f_tether_drag = np.zeros(3)
    is_residual_below_tol = False
    is_run_only_1_time_step = True
    coupling_method = config_dict["coupling_method"]
    n_chordwise_aero_nodes = config_dict["aero"]["n_chordwise_aero_nodes"]
    aero_structural_max_iter = config_dict["aero_structural"]["max_iter"]

    ## actuation
    len_wing_rest_length = len(wing_rest_lengths_initial)
    index_depower_tape = (
        len_wing_rest_length + config_kite_dict["bridle"]["depower_tape_index"]
    )
    initial_length_depower_tape = params["l0"][index_depower_tape]
    depower_tape_extension_step = config_dict["depower_tape_extension_step"]
    depower_tape_final_extension = config_dict["depower_tape_final_extension"]
    index_steering_tape_left = (
        len_wing_rest_length + config_kite_dict["bridle"]["left_steering_tape_index"]
    )
    index_steering_tape_right = (
        len_wing_rest_length + config_kite_dict["bridle"]["right_steering_tape_index"]
    )
    initial_length_steering_tape_right = params["l0"][index_steering_tape_right]
    steering_tape_extension_step = config_dict["steering_tape_extension_step"]
    steering_tape_final_extension = config_dict["steering_tape_final_extension"]
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
    ) = coupling_struc2aero.extract_wingpanel_corners_aero_orderded(
        points, np.array(config_kite_dict["plate_point_indices"])
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
        residual_f = psystem.f_int + f_external
        residual_f_list.append(np.linalg.norm(np.abs(residual_f)))
        print(
            f"i:{i}, t-step:{step:.2f}, residual_f: {np.linalg.norm(residual_f):.3f}N (aero: {end_time_f_ext-begin_time_f_ext:.3f}s, struc: {end_time_f_int-begin_time_f_int:.3f}s)"
        )

        # Update unified tracking dataframe (replaces position update)
        update_tracking(
            tracking_df,
            step,
            psystem,
            points,
            points_prev,
            force_aero,
            force_external,
            residual_f,
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

    # 1. Filter out rows with NaNs in all [xyz] columns
    tracking_without_na = tracking_df.dropna(
        how="all",
        subset=tracking_df.columns[tracking_df.columns.str.contains("[xyz]")],
        axis=0,
    )

    # 2. Gather scalar metadata
    aero_structural_total_time = time.time() - start_time
    num_of_rows = tracking_without_na.shape[0]
    wing_rest_lengths = psystem.extract_rest_length[:len_wing_rest_length]
    bridle_rest_lengths = psystem.extract_rest_length[len_wing_rest_length:]

    # 3. Write to HDF5 (one file holds table + metadata)
    h5_path = results_dir / "sim_output.h5"
    with pd.HDFStore(h5_path, mode="w") as store:
        # store the filtered tracking table
        store["tracking_no_na"] = tracking_without_na

        # attach metadata as attributes on that node
        meta = {
            "total_time_s": aero_structural_total_time,
            "num_iterations": i,
            "converged": is_convergence,
            "n_data_points": num_of_rows,
            "wing_rest_lengths": wing_rest_lengths.tolist(),
            "bridle_rest_lengths": bridle_rest_lengths.tolist(),
        }
        store.get_storer("tracking_no_na").attrs.metadata = meta

    print(f"\nWrote {num_of_rows} valid snapshots and metadata to {h5_path}")

    # 4. Return path (or load back if you need Python objects)
    return {"hdf5_path": str(h5_path)}

    # from VSM.interactive import interactive_plot

    # points_wing_segment_corners_aero_orderded = points[
    #     index_transformation_struc_to_aero
    # ]
    # force_aero_wing_VSM, body_aero = run_vsm_package(
    #     body_aero=body_aero,
    #     solver=vsm_solver,
    #     le_arr=points_wing_segment_corners_aero_orderded[0::2, :],
    #     te_arr=points_wing_segment_corners_aero_orderded[1::2, :],
    #     va_vector=vel_app,
    #     aero_input_type="reuse_initial_polar_data",
    # )

    # interactive_plot(
    #     body_aero,
    #     vel=np.linalg.norm(vel_app),
    #     angle_of_attack=10,
    #     side_slip=0,
    #     yaw_rate=0,
    #     is_with_aerodynamic_details=True,
    #     title="TUDELFT_V3_KITE",
    # )


# Example of post-processing function to analyze convergence
def analyze_convergence(tracking_df):
    """
    Analyze the convergence behavior of the simulation

    Args:
        tracking_df: Unified tracking dataframe with residual data

    Returns:
        Dictionary with convergence metrics
    """
    # Calculate convergence rate in the last iterations
    final_iterations = min(20, len(tracking_df))
    if final_iterations > 5:
        residual_values = tracking_df["residual_norm"].iloc[-final_iterations:].values
        iterations = np.arange(final_iterations)

        # Fit exponential decay to estimate convergence rate
        # log(residual) = log(initial) - rate * iteration
        if np.all(residual_values > 0):
            log_residuals = np.log(residual_values)

            # Simple linear regression
            A = np.vstack([iterations, np.ones(len(iterations))]).T
            rate, log_initial = np.linalg.lstsq(A, log_residuals, rcond=None)[0]

            convergence_rate = -rate  # Negative because we expect decreasing residuals
        else:
            convergence_rate = np.nan
    else:
        convergence_rate = np.nan

    # Calculate oscillation metrics (standard deviation of residual changes)
    if len(tracking_df) > 5:
        residual_changes = np.diff(tracking_df["residual_norm"].values)
        oscillation_metric = (
            np.std(residual_changes) / np.mean(np.abs(residual_changes))
            if np.mean(np.abs(residual_changes)) > 0
            else 0
        )
    else:
        oscillation_metric = np.nan

    # Calculate force balance metrics
    # Average ratio of internal to external forces in final state
    last_idx = tracking_df.index[-1]
    force_columns = tracking_df.filter(regex="^(int|ext)_[xyz]").columns
    int_cols = [col for col in force_columns if col.startswith("int_")]
    ext_cols = [col for col in force_columns if col.startswith("ext_")]

    # Sum of force magnitudes in final state
    int_force_sum = np.sqrt(
        np.sum([tracking_df.loc[last_idx, col] ** 2 for col in int_cols])
    )
    ext_force_sum = np.sqrt(
        np.sum([tracking_df.loc[last_idx, col] ** 2 for col in ext_cols])
    )
    force_balance_ratio = int_force_sum / ext_force_sum if ext_force_sum > 0 else np.nan

    return {
        "convergence_rate": convergence_rate,
        "oscillation_metric": oscillation_metric,
        "min_residual": tracking_df["residual_norm"].min(),
        "max_residual": tracking_df["residual_norm"].max(),
        "final_residual": tracking_df["residual_norm"].iloc[-1],
        "iterations_to_converge": len(tracking_df),
        "force_balance_ratio": force_balance_ratio,
    }


# Example of using the tracking data for visualization
def plot_simulation_data(tracking_df, node_indices=None, save_path=None):
    """
    Create plots to visualize the simulation results

    Args:
        tracking_df: Unified tracking dataframe
        node_indices: List of node indices to plot (if None, plots summary data)
        save_path: Path to save figures (if None, displays figures)
    """
    import matplotlib.pyplot as plt

    # If no specific nodes are requested, use a default selection or summary
    if node_indices is None:
        node_indices = [1]  # Default to first node

    # Create figure for residuals
    plt.figure(figsize=(10, 6))
    plt.semilogy(
        tracking_df.index, tracking_df["residual_norm"], "b-", label="Residual Norm"
    )
    plt.semilogy(
        tracking_df.index, tracking_df["max_residual"], "r--", label="Max Residual"
    )
    plt.semilogy(
        tracking_df.index,
        tracking_df["pos_change_norm"],
        "g-.",
        label="Position Change",
    )
    plt.grid(True, which="both", ls="--")
    plt.xlabel("Time Step")
    plt.ylabel("Residual (log scale)")
    plt.title("Convergence History")
    plt.legend()
    if save_path:
        plt.savefig(f"{save_path}/residuals.png", dpi=300, bbox_inches="tight")

    # Create figure for position trajectory of selected nodes
    plt.figure(figsize=(10, 8))
    for node_idx in node_indices:
        x = tracking_df[f"x{node_idx}"]
        y = tracking_df[f"y{node_idx}"]
        z = tracking_df[f"z{node_idx}"]

        ax = plt.axes(projection="3d")
        ax.plot3D(x, y, z, label=f"Node {node_idx}")
        ax.set_xlabel("X Position")
        ax.set_ylabel("Y Position")
        ax.set_zlabel("Z Position")
        ax.set_title(f"Position Trajectory for Selected Nodes")
        ax.legend()

    if save_path:
        plt.savefig(
            f"{save_path}/position_trajectory.png", dpi=300, bbox_inches="tight"
        )

    # Create figure for force components for selected nodes
    for node_idx in node_indices:
        plt.figure(figsize=(12, 8))

        # Create 3 subplots for x, y, z components
        for i, component in enumerate(["x", "y", "z"]):
            plt.subplot(3, 1, i + 1)
            plt.plot(
                tracking_df.index,
                tracking_df[f"aero_{component}{node_idx}"],
                "b-",
                label="Aero Force",
            )
            plt.plot(
                tracking_df.index,
                tracking_df[f"int_{component}{node_idx}"],
                "r--",
                label="Internal Force",
            )
            plt.plot(
                tracking_df.index,
                tracking_df[f"ext_{component}{node_idx}"],
                "g-.",
                label="External Force",
            )
            plt.grid(True)
            plt.xlabel("Time Step")
            plt.ylabel(f"{component.upper()} Force Component")
            if i == 0:
                plt.title(f"Force Components for Node {node_idx}")
            plt.legend()

        plt.tight_layout()
        if save_path:
            plt.savefig(
                f"{save_path}/forces_node{node_idx}.png", dpi=300, bbox_inches="tight"
            )

    if not save_path:
        plt.show()


# %%


# def run_aerostructural_solver(config_dict, config_kite_dict, PROJECT_DIR):
#     """Runs the aero-structural solver for the given input parameters"""

#     ### SOLVER INITIALISATION
#     (
#         points,
#         wing_connectivity,
#         bridle_connectivity,
#         kite_connectivity,
#         wing_rest_lengths_initial,
#         rest_lengths,
#         m_array,
#     ) = initialising_solver(config_kite_dict)

#     ### AERO initialisation -- VSM
#     body_aero, vsm_solver = initialize_vsm(
#         geometry_csv_path=Path(PROJECT_DIR)
#         / "data"
#         / "V3_25"
#         / "wing_geometry_from_CAD.csv",
#         polar_data_dir=Path(PROJECT_DIR) / "data" / "V3_25" / "2D_polars_from_CFD",
#         n_panels=9,
#         spanwise_panel_distribution="uniform",
#         is_half_wing=True,
#         is_with_corrected_polar=True,
#     )

#     ## STRUC initialisation -- PSS
#     psystem, params = instantiate_psystem(
#         config_dict,
#         config_kite_dict,
#         points,
#         wing_connectivity,
#         kite_connectivity,
#         rest_lengths,
#         m_array,
#     )

#     ##TODO: add a function to set the initial position of the kite
#     # if config_dict["is_with_initial_plot"]:
#     #     plot the geometry

#     ## setting up the position-dataframe
#     t_vector = np.linspace(
#         params["dt"], params["t_steps"] * params["dt"], params["t_steps"]
#     )
#     x = {}
#     for i in range(params["n"]):
#         x[f"x{i + 1}"] = np.zeros(len(t_vector))
#         x[f"y{i + 1}"] = np.zeros(len(t_vector))
#         x[f"z{i + 1}"] = np.zeros(len(t_vector))

#     position = pd.DataFrame(index=t_vector, columns=x)
#     aero_structural_tol = params["aerostructural_tol"]

#     # INITIALISATION OF VARIABLES
#     ## parameters used in the loop
#     vel_app = np.array(config_dict["vel_wind"]) - np.array(config_dict["vel_kite"])
#     start_time = time.time()
#     is_convergence = False
#     residual_f_list = []
#     f_tether_drag = np.zeros(3)
#     is_residual_below_tol = False
#     is_run_only_1_time_step = True
#     coupling_method = config_dict["coupling_method"]
#     n_chordwise_aero_nodes = config_dict["aero"]["n_chordwise_aero_nodes"]
#     aero_structural_max_iter = config_dict["aero_structural"]["max_iter"]

#     ## actuation
#     len_wing_rest_length = len(wing_rest_lengths_initial)
#     index_depower_tape = (
#         len_wing_rest_length + config_kite_dict["bridle"]["depower_tape_index"]
#     )
#     initial_length_depower_tape = params["l0"][index_depower_tape]
#     depower_tape_extension_step = config_dict["depower_tape_extension_step"]
#     depower_tape_final_extension = config_dict["depower_tape_final_extension"]
#     index_steering_tape_left = (
#         len_wing_rest_length + config_kite_dict["bridle"]["left_steering_tape_index"]
#     )
#     index_steering_tape_right = (
#         len_wing_rest_length + config_kite_dict["bridle"]["right_steering_tape_index"]
#     )
#     initial_length_steering_tape_right = params["l0"][index_steering_tape_right]
#     steering_tape_extension_step = config_dict["steering_tape_extension_step"]
#     steering_tape_final_extension = config_dict["steering_tape_final_extension"]
#     ### printing initial lengths
#     print(
#         f"Initial depower tape length: {psystem.extract_rest_length[index_depower_tape]:.3f}m"
#     )
#     print(
#         f"Desired depower tape length: {initial_length_depower_tape + depower_tape_final_extension:.3f}m"
#     )
#     print(
#         f"Initial steering tape length right: {psystem.extract_rest_length[index_steering_tape_right]:.3f}m"
#     )
#     print(
#         f"Desired steering tape length right: {initial_length_steering_tape_right + steering_tape_final_extension:.3f}m"
#     )

#     if config_dict["is_with_gravity"]:
#         force_gravity = np.array(
#             [np.array(config_dict["grav_constant"]) * m_pt for m_pt in m_array]
#         )
#     else:
#         force_gravity = np.zeros(points.shape)

#     ## coupling initialisation
#     (
#         points_wing_segment_corners_aero_orderded,
#         index_transformation_struc_to_aero,
#     ) = coupling_struc2aero.extract_wingpanel_corners_aero_orderded(
#         points, np.array(config_kite_dict["plate_point_indices"])
#     )
#     logging.info(
#         f"points_wing_segment_corners_aero_orderded: {points_wing_segment_corners_aero_orderded}"
#     )
#     n_ribs = int(len(points_wing_segment_corners_aero_orderded) / 2)
#     n_panels = int(n_ribs - 1)
#     logging.info(f"n_ribs: {n_ribs}, n_panels: {n_panels}")

#     ### creating a dict with key value pairs, to transform from aero to struc
#     index_transformation_aero_to_struc_dict = {}
#     for i, value in enumerate(index_transformation_struc_to_aero):
#         index_transformation_aero_to_struc_dict[value] = i

#     print(f" ")
#     print(f"Running aero-structural simulation")
#     print(f"----------------------------------- ")

#     ######################################################################
#     # SIMULATION LOOP
#     ######################################################################
#     ## propagating the simulation for each timestep and saving results
#     for i, step in enumerate(t_vector):

#         ## external force
#         begin_time_f_ext = time.time()
#         ### STRUC --> AERO
#         # Ordering the points from left to right (-y to +y), the desired aero-ordering
#         points_wing_segment_corners_aero_orderded = points[
#             index_transformation_struc_to_aero
#         ]
#         ### AERO
#         force_aero_wing_VSM, body_aero = run_vsm_package(
#             body_aero=body_aero,
#             solver=vsm_solver,
#             le_arr=points_wing_segment_corners_aero_orderded[0::2, :],
#             te_arr=points_wing_segment_corners_aero_orderded[1::2, :],
#             va_vector=vel_app,
#             aero_input_type="reuse_initial_polar_data",
#         )
#         # reshuffle forces, instead of left-to-right, we make them right-to-left
#         force_aero_wing_VSM = force_aero_wing_VSM[::-1]
#         ### AERO --> STRUC
#         if coupling_method == "NN":
#             force_aero_wing = coupling_aero2struc.aero2struc_NN_vsm(
#                 n_chordwise_aero_nodes,
#                 body_aero,
#                 force_aero_wing_VSM,
#                 points_wing_segment_corners_aero_orderded,
#                 index_transformation_aero_to_struc_dict,
#                 points,
#             )
#         else:
#             raise ValueError("Coupling method not recognized; wrong name or typo")

#         # TODO: get bridle line forces back in to play
#         force_aero_bridle = np.zeros(points.shape)
#         force_aero = force_aero_wing + force_aero_bridle

#         ## summing up
#         force_external = force_aero + force_gravity

#         ### STRUC
#         ### f_external is flat, and force_external is 2D, ##TODO: could add this to the name?
#         f_external = force_external.flatten()
#         end_time_f_ext = time.time()
#         begin_time_f_int = time.time()
#         psystem = run_pss(psystem, params, f_external)
#         position.loc[step], _ = psystem.x_v_current
#         end_time_f_int = time.time()

#         logging.debug(f"position.loc[step].shape: {position.loc[step].shape}")
#         logging.debug(f"internal force: {psystem.f_int}")
#         logging.debug(f"external force: {f_external}")

#         # # TODO: ideally you don't need this here and have aero-also work with the flat format, as this should be faster
#         # # saving points in different format
#         # points = psystem.x_current_2D
#         ## TODO: replacing this function inside the src code to inhere
#         points = np.array([particle.x for particle in psystem.particles])
#         residual_f = psystem.f_int + f_external
#         residual_f_list.append(np.linalg.norm(np.abs(residual_f)))
#         print(
#             f"i:{i}, t-step:{step:.2f}, residual_f: {np.linalg.norm(residual_f):.3f}N (aero: {end_time_f_ext-begin_time_f_ext:.3f}s, struc: {end_time_f_int-begin_time_f_int:.3f}s)"
#         )

#         ## calculating delta tape lengths
#         delta_depower_tape = (
#             psystem.extract_rest_length[index_depower_tape]
#             - initial_length_depower_tape
#         )

#         ### All the convergence checks, are be done in if-elif because only 1 should hold at once
#         # if convergence (residual below set tolerance)
#         if np.linalg.norm(residual_f) <= aero_structural_tol:
#             is_residual_below_tol = True
#         if np.linalg.norm(residual_f) <= aero_structural_tol and i > 50:
#             is_residual_below_tol = True
#             is_convergence = True

#         # if residual forces are NaN
#         elif np.isnan(np.linalg.norm(residual_f)):
#             is_convergence = False
#             print("Classic PS diverged - residual force is NaN")
#             break
#         # if residual forces are not changing anymore
#         elif (
#             i > 200
#             and np.abs(np.mean(residual_f_list[i - 25]) - residual_f_list[i]) < 1
#             and np.abs(np.mean(residual_f_list[i - 10]) - residual_f_list[i]) < 1
#             and np.abs(np.mean(residual_f_list[i - 5]) - residual_f_list[i]) < 1
#             and np.abs(np.mean(residual_f_list[i - 2]) - residual_f_list[i]) < 1
#         ):
#             is_convergence = False
#             print("Classic PS non-converging - residual no longer changes")
#             break
#         # if to many iterations are needed
#         elif i > aero_structural_max_iter:
#             is_convergence = False
#             print(
#                 f"Classic PS non-converging - more than max ({aero_structural_max_iter}) iterations needed"
#             )
#             break
#         # special case for running the simulation for only one timestep
#         elif is_run_only_1_time_step:
#             break

#         ## if convergenced and not yet at desired depower-tape length
#         elif (
#             is_residual_below_tol and delta_depower_tape <= depower_tape_final_extension
#         ):
#             psystem.update_rest_length(index_depower_tape, depower_tape_extension_step)
#             delta_depower_tape = (
#                 psystem.extract_rest_length[index_depower_tape]
#                 - initial_length_depower_tape
#             )
#             print(
#                 f"||--- delta l_d: {delta_depower_tape:.3f}m | new l_d: {psystem.extract_rest_length[index_depower_tape]:.3f}m"
#             )
#             is_residual_below_tol = False

#         ## if convergenced and all changes are made
#         elif is_residual_below_tol:
#             is_convergence = True

#         if is_convergence:
#             break
#     ######################################################################
#     ## END OF SIMULATION FOR LOOP
#     ######################################################################

#     # defining post_processing_output
#     aero_structural_total_time = time.time() - start_time
#     wing_rest_lengths = psystem.extract_rest_length[0:len_wing_rest_length]
#     bridle_rest_lengths = psystem.extract_rest_length[len_wing_rest_length:]
#     position_without_na = position.dropna(
#         how="all",
#         subset=position.columns[position.columns.str.contains("[xyz]")],
#         axis=0,
#     )
#     num_of_rows = position_without_na.shape[0]

#     sim_output = {
#         "points": points,
#         "position": position,
#         "aero_structural_total_time": aero_structural_total_time,
#         "num_of_iterations": i,
#         "is_convergence": is_convergence,
#         "wing_rest_lengths": wing_rest_lengths,
#         "bridle_rest_lengths": bridle_rest_lengths,
#         # print data
#         "is_convergence": is_convergence,
#         "num_of_iterations": i,
#         "vel_app": vel_app,
#         ## aero_structural_total_time
#         "residual_f_including_fixed_nodes": psystem.f,
#         "residual_f": residual_f,
#         "f_internal": psystem.f_int,
#         "f_external": f_external,
#         "force_aero": force_aero,
#         "force_aero_wing": force_aero_wing,
#         "force_aero_bridle": force_aero_bridle,
#         "f_tether_drag": f_tether_drag,
#         "force_gravity": force_gravity,
#         ## wing_rest_lengths
#         ## bridle_rest_lengths
#         # plot data
#         "wingpanels": np.zeros(points.shape),  # wingpanels,
#         "controlpoints": np.zeros(points.shape),  # controlpoints,
#         "rings": np.zeros(points.shape),  # rings,
#         "coord_L": np.zeros(points.shape),  # coord_L,
#         "F_rel": np.zeros(points.shape),  # F_rel,
#         ## wing_rest_lengths
#         ## bridle_rest_lengths
#         # animation data
#         "position_without_na": position_without_na,
#         "num_of_rows": num_of_rows,
#         ## wing_rest_lengths
#         ## bridle_rest_lengths
#     }

#     # from VSM.interactive import interactive_plot

#     # points_wing_segment_corners_aero_orderded = points[
#     #     index_transformation_struc_to_aero
#     # ]
#     # force_aero_wing_VSM, body_aero = run_vsm_package(
#     #     body_aero=body_aero,
#     #     solver=vsm_solver,
#     #     le_arr=points_wing_segment_corners_aero_orderded[0::2, :],
#     #     te_arr=points_wing_segment_corners_aero_orderded[1::2, :],
#     #     va_vector=vel_app,
#     #     aero_input_type="reuse_initial_polar_data",
#     # )

#     # interactive_plot(
#     #     body_aero,
#     #     vel=np.linalg.norm(vel_app),
#     #     angle_of_attack=10,
#     #     side_slip=0,
#     #     yaw_rate=0,
#     #     is_with_aerodynamic_details=True,
#     #     title="TUDELFT_V3_KITE",
#     # )
#     return sim_output
