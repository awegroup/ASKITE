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
from kitesim.solver.vsm_functions import (
    initialize_vsm,
    run_vsm_package,
    plot_vsm_geometry,
)
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


def initialize_aero2struc_mapping(
    panels: np.ndarray,
    struc_nodes: np.ndarray,
    struc_le_idx_list: np.ndarray,
    struc_te_idx_list: np.ndarray,
) -> np.ndarray:
    """
    For each panel CP, find the two LE and two TE structural‐node indices
    whose y-coordinates bracket the CP’s y. Returns an (n_panels, 4) array
    of [le_lo, le_hi, te_lo, te_hi].
    """
    # extract and sort LE candidates by their y
    le_coords = struc_nodes[struc_le_idx_list]
    le_y = le_coords[:, 1]
    le_order = np.argsort(le_y)
    le_sorted_idx = np.array(struc_le_idx_list)[le_order]
    le_sorted_y = le_y[le_order]

    # same for TE
    te_coords = struc_nodes[struc_te_idx_list]
    te_y = te_coords[:, 1]
    te_order = np.argsort(te_y)
    te_sorted_idx = np.array(struc_te_idx_list)[te_order]
    te_sorted_y = te_y[te_order]

    n = len(panels)
    mapping = np.zeros((n, 4), dtype=int)

    for i, panel in enumerate(panels):
        y = panel.aerodynamic_center[1]
        # LE insertion point
        hi_le = np.searchsorted(le_sorted_y, y)
        lo_le = np.clip(hi_le - 1, 0, len(le_sorted_y) - 1)
        hi_le = np.clip(hi_le, 0, len(le_sorted_y) - 1)

        # TE insertion
        hi_te = np.searchsorted(te_sorted_y, y)
        lo_te = np.clip(hi_te - 1, 0, len(te_sorted_y) - 1)
        hi_te = np.clip(hi_te, 0, len(te_sorted_y) - 1)

        mapping[i, :] = [
            le_sorted_idx[lo_le],
            le_sorted_idx[hi_le],
            te_sorted_idx[lo_te],
            te_sorted_idx[hi_te],
        ]

    return mapping


def create_struc_nodes(config_kite_dict):
    """Creates the nodes for the structural system"""

    # Extract airfoil
    af = config_kite_dict["airfoils"]
    headers = af["headers"]  # e.g. ["LE_x", "LE_y", ...]
    data = af["data"]  # list of lists
    airfoil_data = config_kite_dict["airfoils"]["data"]
    df_airfoil = pd.DataFrame(data, columns=headers)

    wing_ci = []
    wing_cj = []

    # 1) Pick out exactly the panel‐indices you’ll represent in the structure:
    n = len(df_airfoil)
    selected = [
        i for i in range(n) if i == 0 or df_airfoil.loc[i, "is_strut"] or i == n - 1
    ]

    # 2) Build your structural‐node list and remember for each panel which
    struc_nodes_list = [[0, 0, 0]]
    #    structural index is its LE and which is its TE
    le_map = {}
    te_map = {}
    struc_le_idx_list = []
    struc_te_idx_list = []
    for i in selected:
        row = df_airfoil.loc[i]
        # append LE
        le_idx = len(struc_nodes_list)
        struc_nodes_list.append(row[["LE_x", "LE_y", "LE_z"]].tolist())
        # append TE
        te_idx = len(struc_nodes_list)
        struc_nodes_list.append(row[["TE_x", "TE_y", "TE_z"]].tolist())
        le_map[i] = le_idx
        te_map[i] = te_idx
        # append to the list of indices
        struc_le_idx_list.append(le_idx)
        struc_te_idx_list.append(te_idx)

    tubular_frame_line_idx_list = []
    te_line_idx_list = []
    # 3) Now build wing_ci/wing_cj with explicit if‐branches
    for pos, i in enumerate(selected):
        curr_le = le_map[i]
        curr_te = te_map[i]

        if pos == 0:
            # --- first panel only has a forward neighbor ---
            # connect LE->TE
            wing_ci.append(curr_le)
            wing_cj.append(curr_te)
            tubular_frame_line_idx_list.append(len(wing_ci) - 1)
            # connect to next panel’s LE (LE[i] -- LE[i+1])
            nxt = selected[pos + 1]
            wing_ci.append(curr_le)
            wing_cj.append(le_map[nxt])
            tubular_frame_line_idx_list.append(len(wing_ci) - 1)
            # connect to next panel’s TE (TE[i] -- TE[i+1])
            wing_ci.append(curr_te)
            wing_cj.append(te_map[nxt])
            te_line_idx_list.append(len(wing_ci) - 1)

        elif pos == len(selected) - 1:
            # --- last panel only has a backward neighbor ---
            # connect LE->TE
            wing_ci.append(curr_le)
            wing_cj.append(curr_te)
            tubular_frame_line_idx_list.append(len(wing_ci) - 1)
            # connect to previous panel’s LE (LE[i] -- LE[i-1])
            prev = selected[pos - 1]
            wing_ci.append(curr_le)
            wing_cj.append(le_map[prev])
            tubular_frame_line_idx_list.append(len(wing_ci) - 1)
            # connect to previous panel’s TE (TE[i] -- TE[i-1])
            wing_ci.append(curr_te)
            wing_cj.append(te_map[prev])
            te_line_idx_list.append(len(wing_ci) - 1)

        else:
            # --- interior strut panel: both backward & forward ---
            # connect LE->TE
            wing_ci.append(curr_le)
            wing_cj.append(curr_te)
            tubular_frame_line_idx_list.append(len(wing_ci) - 1)

            prev = selected[pos - 1]
            nxt = selected[pos + 1]
            # LE -> prev_LE and LE -> next_LE
            wing_ci.append(curr_le)
            wing_cj.append(le_map[prev])
            tubular_frame_line_idx_list.append(len(wing_ci) - 1)
            wing_ci.append(curr_le)
            wing_cj.append(le_map[nxt])
            tubular_frame_line_idx_list.append(len(wing_ci) - 1)
            # TE -> prev_TE and TE -> next_TE
            wing_ci.append(curr_te)
            wing_cj.append(te_map[prev])
            te_line_idx_list.append(len(wing_ci) - 1)
            wing_ci.append(curr_te)
            wing_cj.append(te_map[nxt])
            te_line_idx_list.append(len(wing_ci) - 1)

    for i, struc_node in enumerate(struc_nodes_list):
        if i == 0:
            continue
        struc_nodes_list[i][0] = struc_node[0] + 0.7
        struc_nodes_list[i][2] = struc_node[2] + 7.3

    if config_kite_dict["is_with_bridle_line_system_surfplan"]:
        df_bridle_line_system_surfplan = pd.DataFrame(
            config_kite_dict["df_bridle_line_system_surfplan"]["data"],
            columns=config_kite_dict["df_bridle_line_system_surfplan"]["headers"],
        )

        # Initialize connectivity lists
        bridle_ci = []
        bridle_cj = []

        # Helper function to compare 3D points (rounded to tolerance)
        def point_key(p, tol=1e-8):
            return tuple(
                np.round(p, decimals=8)
            )  # prevent floating point uniqueness issues

        # Build a dictionary of existing nodes for fast lookup
        existing_node_map = {
            point_key(node): idx for idx, node in enumerate(struc_nodes_list)
        }

        # Iterate over each bridle line
        for _, row in df_bridle_line_system_surfplan.iterrows():
            # Extract p1 and p2
            p1 = [row["p1_x"], row["p1_y"], row["p1_z"]]
            p2 = [row["p2_x"], row["p2_y"], row["p2_z"]]

            # Deduplicate and register p1
            key1 = point_key(p1)
            if key1 in existing_node_map:
                idx1 = existing_node_map[key1]
            else:
                idx1 = len(struc_nodes_list)
                struc_nodes_list.append(p1)
                existing_node_map[key1] = idx1

            # Deduplicate and register p2
            key2 = point_key(p2)
            if key2 in existing_node_map:
                idx2 = existing_node_map[key2]
            else:
                idx2 = len(struc_nodes_list)
                struc_nodes_list.append(p2)
                existing_node_map[key2] = idx2

            # Store connectivity
            bridle_ci.append(idx1)
            bridle_cj.append(idx2)
    else:
        # Extract brilde_nodes
        bridle_nodes = config_kite_dict["bridle_nodes"]
        headers = bridle_nodes["headers"]  # e.g. ["LE_x", "LE_y", ...]
        data = bridle_nodes["data"]  # list of lists
        bridle_nodes_data = config_kite_dict["bridle_nodes"]["data"]
        df_bridle_nodes = pd.DataFrame(data, columns=headers)

        pulley_point_indices = []
        for _, row in df_bridle_nodes.iterrows():
            struc_nodes_list.append(row[["x", "y", "z"]].tolist())
            if row["type"] == "pulley":
                pulley_point_indices.append(len(struc_nodes_list) - 1)

        bridle_ci = config_kite_dict["bridle_lines"]["ci"]
        bridle_cj = config_kite_dict["bridle_lines"]["cj"]

        # Given that we know

    return (
        np.array(struc_nodes_list),
        wing_ci,
        wing_cj,
        bridle_ci,
        bridle_cj,
        struc_le_idx_list,
        struc_te_idx_list,
        pulley_point_indices,
        tubular_frame_line_idx_list,
        te_line_idx_list,
    )


def distribute_mass(
    points_ini,
    bridle_ci,
    bridle_cj,
    wing_ci,
    wing_cj,
    pulley_point_indices,
    config_kite_dict,
):

    WING_MASS = config_kite_dict["wing_mass"]
    BRIDLE_RHO = config_kite_dict["bridle"]["density"]
    BRIDLE_DIAMETER = config_kite_dict["bridle"]["diameter"]
    KCU_MASS = config_kite_dict["kcu"]["mass"]
    KCU_INDEX = config_kite_dict["kcu"]["index"]
    PULLEY_MASS = config_kite_dict["pulley"]["mass"]

    # Define a spring force matrix of the right size
    node_masses = np.zeros(
        points_ini.shape[0]
    )  # Initialising with zero matrix in same shape as points

    ## Bridle lines
    for idx, (idx_bridle_node_i, idx_bridle_node_j) in enumerate(
        zip(bridle_ci, bridle_cj)
    ):  # loop through each bridle line
        # Calculate the length of the bridle line
        length_bridle = np.linalg.norm(
            points_ini[idx_bridle_node_i] - points_ini[idx_bridle_node_j]
        )
        # Calculate the mass of the bridle line
        mass_bridle = BRIDLE_RHO * np.pi * (BRIDLE_DIAMETER / 2) ** 2 * length_bridle
        # Add the mass of the bridle line to the nodes
        node_masses[idx_bridle_node_i] += mass_bridle / 2
        node_masses[idx_bridle_node_j] += mass_bridle / 2

    # print(f'bridle node_masses: {np.sum(node_masses)}')
    ## Pulleys
    for idx in pulley_point_indices:
        node_masses[idx] += PULLEY_MASS

    # print(f'pulley& bridle node_masses: {np.sum(node_masses)}')
    ## KCU
    node_masses[KCU_INDEX] += KCU_MASS

    # print(f'pulley& bridle & KCU node_masses: {np.sum(node_masses)}')

    ## Wing
    for idx, (idx_wing_node_i, idx_wing_node_i) in enumerate(
        zip(set(wing_ci), set(wing_cj))
    ):  # making them sets, to have it unique
        node_masses[idx] += WING_MASS / len(set(wing_ci))

    # print(f'pulley& bridle & KCU & wing node_masses: {np.sum(node_masses)}')

    return node_masses


def calculate_edge_lengths(ci, cj, pos):
    """returns the edge lengths between the nodes with index ci and cj
    for the given positions pos
    input : ci,cj,pos
    output: springL"""
    springL = np.zeros(ci.shape)
    for idx, (ci, cj) in enumerate(zip(ci, cj)):
        springL[idx] = np.linalg.norm(pos[cj, :] - pos[ci, :])
    return springL


def run_aerostructural_solver(config_dict, config_kite_dict, PROJECT_DIR, results_dir):
    """Runs the aero-structural solver for the given input parameters"""

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
    ) = create_struc_nodes(config_kite_dict)

    wing_connectivity = np.column_stack(
        (
            wing_ci,
            wing_cj,
        )
    )
    bridle_connectivity = np.column_stack(
        (
            bridle_ci,
            bridle_cj,
        )
    )
    kite_connectivity = np.vstack((wing_connectivity, bridle_connectivity))

    ### AERO initialisation -- VSM
    # body_aero, vsm_solver = initialize_vsm(
    #     geometry_csv_path=Path(PROJECT_DIR)
    #     / "data"
    #     / "V3_25"
    #     / "wing_geometry_from_CAD.csv",
    #     polar_data_dir=Path(PROJECT_DIR) / "data" / "V3_25" / "2D_polars_from_CFD",
    #     n_panels=9,
    #     spanwise_panel_distribution="uniform",
    #     is_half_wing=True,
    #     is_with_corrected_polar=True,
    # )
    # body_aero, vsm_solver = initialize_vsm(
    #     geometry_csv_path=Path(PROJECT_DIR)
    #     / "data"
    #     / "V3_25"
    #     / "wing_geometry_CAD.csv",
    #     polar_data_dir=Path(PROJECT_DIR) / "data" / "V3_25" / "2D_polars_from_CFD",
    #     n_panels=9,
    #     spanwise_panel_distribution="uniform",
    #     is_half_wing=False,
    #     is_with_corrected_polar=False,
    # )
    body_aero, vsm_solver = initialize_vsm(
        config_kite=config_kite_dict,
        n_panels=9,
        spanwise_panel_distribution="uniform",
    )
    vel_app = np.array(config_dict["vel_wind"]) - np.array(config_dict["vel_kite"])
    body_aero.va = (vel_app, 0)
    # solve the problem
    results = vsm_solver.solve(body_aero)
    f_aero = np.array(results["F_distribution"])
    f_y = np.sum([force[1] for force in f_aero])
    logging.info(f"Sum of forces in y-direction: {f_y}")

    wing = body_aero.wings[0]
    new_sections = wing.refine_aerodynamic_mesh()
    initial_polar_data = []
    for new_section in new_sections:
        initial_polar_data.append(new_section.aero_input)

    # print(f"initial_polar_data: {initial_polar_data}")
    aero2struc_mapping = initialize_aero2struc_mapping(
        body_aero.panels,
        struc_nodes,
        struc_le_idx_list,
        struc_te_idx_list,
    )

    m_array = distribute_mass(
        struc_nodes,
        bridle_ci,
        bridle_cj,
        wing_ci,
        wing_cj,
        pulley_point_indices,
        config_kite_dict,
    )
    ## computing rest lengths
    bridle_rest_lengths_initial = np.array(
        calculate_edge_lengths(
            np.array(bridle_ci),
            np.array(bridle_cj),
            struc_nodes,
        )
    )
    wing_rest_lengths_initial = np.array(
        calculate_edge_lengths(
            np.array(wing_ci),
            np.array(wing_cj),
            struc_nodes,
        )
    )
    rest_lengths = np.concatenate(
        (wing_rest_lengths_initial, bridle_rest_lengths_initial)
    )

    # ## INITIALISATION
    # (
    #     points,
    #     le_idx,
    #     te_idx,
    #     wing_connectivity,
    #     bridle_connectivity,
    #     kite_connectivity,
    #     wing_rest_lengths_initial,
    #     rest_lengths,
    #     m_array,
    # ) = initialising_solver(config_kite_dict)

    ## STRUC initialisation -- PSS
    psystem, params, pss_kite_connectivity = instantiate_psystem(
        config_dict,
        config_kite_dict,
        struc_nodes,
        wing_connectivity,
        kite_connectivity,
        rest_lengths,
        m_array,
        tubular_frame_line_idx_list,
        te_line_idx_list,
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
        force_gravity = np.zeros(struc_nodes.shape)

    # PreLoop Initialisation
    t_vector = np.linspace(
        params["dt"], params["t_steps"] * params["dt"], params["t_steps"]
    )
    tracking = setup_tracking_arrays(len(struc_nodes), t_vector)
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
            # TODO: Write something more complex to handle more than 1 on 1 mapping
            le_arr = struc_nodes[struc_le_idx_list]
            te_arr = struc_nodes[struc_te_idx_list]

            ### AERO
            f_aero_wing_VSM, body_aero, results_aero = run_vsm_package(
                body_aero=body_aero,
                solver=vsm_solver,
                le_arr=le_arr,
                te_arr=te_arr,
                va_vector=vel_app,
                aero_input_type="reuse_initial_polar_data",
                initial_polar_data=initial_polar_data,
            )
            f_y = np.sum([force[1] for force in f_aero_wing_VSM])
            logging.info(f"Sum of forces in y-direction: {f_y}")

            ### AERO --> STRUC
            if coupling_method == "NN":
                # TODO: remove hardcoded values
                f_aero_wing = coupling_aero2struc.aero2struc_NN_vsm(
                    f_aero_wing_VSM,  # (n_panels,3)
                    struc_nodes,  # (n_struc,3)
                    np.array(results_aero["panel_cp_locations"]),  # (n_panels,3)
                    aero2struc_mapping,  # (n_panels,4)
                    p=2,
                    eps=1e-6,
                    is_with_coupling_plot=config_dict["is_with_coupling_plot"],
                )
            else:
                raise ValueError("Coupling method not recognized; wrong name or typo")

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
            struc_nodes = np.array([particle.x for particle in psystem.particles])
            f_residual = psystem.f_int + f_external
            f_residual_list.append(np.linalg.norm(np.abs(f_residual)))

            # Update unified tracking dataframe (replaces position update)
            update_tracking_arrays(
                tracking,
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
        "n_iter": n_iter,
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
