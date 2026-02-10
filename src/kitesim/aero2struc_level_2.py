import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from kitesim.plotting import plot_aerodynamic_forces_chordwise_distributed


def build_ordered_sections(struc_nodes, canopy_sections, strut_sections):
    """
    Build spanwise-ordered section lists and their coordinates.

    Args:
        struc_nodes (np.ndarray): Structural node positions (n_nodes, 3).
        canopy_sections (list[list[int]]): Chordwise node indices per canopy section.
        strut_sections (list[list[int]]): Chordwise node indices per strut section.

    Returns:
        tuple: (sections, section_coords)
            - sections: list of index lists ordered by LE y-coordinate.
            - section_coords: list of arrays with coordinates for each section.
    """
    sections = canopy_sections + strut_sections
    # sort by LE y (or by section[0] if that encodes span order)
    sections = sorted(sections, key=lambda sec: struc_nodes[sec[0]][1])
    # return indices + coords
    section_coords = [struc_nodes[np.array(sec)] for sec in sections]
    return sections, section_coords


def map_aero_forces_to_struct_nodes(aero_points, aero_forces, struc_nodes, sections):
    """
    Map distributed aerodynamic forces onto structural beam nodes using
    bilinear interpolation across two spanwise-bracketing sections.

    For each aerodynamic force application point:
      1. Find the two sections that bracket it in spanwise (y) direction.
      2. Compute spanwise weight eta in [0, 1] between those sections.
      3. Within each section, find the closest beam segment and project
         onto it → chordwise weight xi in [0, 1].
      4. Distribute force to the resulting 4 nodes (2 per section) via
         combined weights:  w_section * (1-xi) and w_section * xi.

    Symmetry preservation: because sections are sorted by y and bracketing
    is deterministic, mirror-symmetric aero points receive mirror-symmetric
    weights, so symmetric loads produce exactly symmetric nodal forces.

    Force is exactly preserved.  Moment error is proportional to the
    out-of-plane offset between the aero application point and the beam
    segments (same as the single-segment version).

    Args:
        aero_points (np.ndarray): (N_aero, 3) coordinates where forces act.
        aero_forces (np.ndarray): (N_aero, 3) force vectors.
        struc_nodes (np.ndarray): (N_struc, 3) structural node coordinates.
        sections (list[list[int]]): Each element is a list of node indices
            forming a chordwise beam (consecutive pairs are segments).
            Must be sorted by ascending LE y-coordinate (as returned by
            build_ordered_sections).

    Returns:
        np.ndarray: (N_struc, 3) lumped force vectors at structural nodes.
    """
    aero_points = np.asarray(aero_points, dtype=float)
    aero_forces = np.asarray(aero_forces, dtype=float)
    n_struc = len(struc_nodes)
    nodal_forces = np.zeros((n_struc, 3), dtype=float)

    n_sections = len(sections)
    if n_sections == 0:
        return nodal_forces

    # LE y-coordinate of each section (sections are assumed sorted by y)
    section_y = np.array([struc_nodes[sec[0]][1] for sec in sections])

    # Pre-build segments per section for fast lookup
    section_segments = []
    for sec in sections:
        segs = []
        for j in range(len(sec) - 1):
            idx_a, idx_b = sec[j], sec[j + 1]
            p_a = struc_nodes[idx_a]
            seg_vec = struc_nodes[idx_b] - p_a
            seg_len_sq = np.dot(seg_vec, seg_vec)
            segs.append((idx_a, idx_b, p_a, seg_vec, seg_len_sq))
        section_segments.append(segs)

    for k in range(len(aero_points)):
        p = aero_points[k]
        f = aero_forces[k]
        y_p = p[1]

        # --- spanwise bracketing ---
        right = int(np.searchsorted(section_y, y_p))
        left = right - 1

        if right >= n_sections:
            right = n_sections - 1
            left = right - 1 if right > 0 else right
        if left < 0:
            left = 0
            right = 1 if n_sections > 1 else 0

        if left == right:
            # single section available
            pairs = [(left, 1.0)]
        else:
            dy = section_y[right] - section_y[left]
            if dy < 1e-30:
                eta = 0.5
            else:
                eta = (y_p - section_y[left]) / dy
            eta = max(0.0, min(1.0, eta))
            pairs = [(left, 1.0 - eta), (right, eta)]

        # --- distribute to both sections ---
        for s_idx, s_w in pairs:
            if s_w < 1e-15:
                continue
            segs = section_segments[s_idx]
            if len(segs) == 0:
                # section has only 1 node → lump everything there
                nodal_forces[sections[s_idx][0]] += s_w * f
                continue

            # find closest segment in this section
            best_dist_sq = np.inf
            best_idx_a = 0
            best_idx_b = 0
            best_xi = 0.0

            for idx_a, idx_b, p_a, seg_vec, seg_len_sq in segs:
                if seg_len_sq < 1e-30:
                    xi = 0.0
                else:
                    xi = np.dot(p - p_a, seg_vec) / seg_len_sq
                xi_c = max(0.0, min(1.0, xi))
                proj = p_a + xi_c * seg_vec
                d_sq = np.sum((p - proj) ** 2)

                if d_sq < best_dist_sq:
                    best_dist_sq = d_sq
                    best_idx_a = idx_a
                    best_idx_b = idx_b
                    best_xi = xi_c

            nodal_forces[best_idx_a] += s_w * (1.0 - best_xi) * f
            nodal_forces[best_idx_b] += s_w * best_xi * f

    return nodal_forces


def verify_force_moment_conservation(
    aero_points, aero_forces, struc_nodes, nodal_forces, ref_point=None
):
    """
    Print a check of total force and moment conservation.

    Args:
        aero_points (np.ndarray): (N_aero, 3)
        aero_forces (np.ndarray): (N_aero, 3)
        struc_nodes (np.ndarray): (N_struc, 3)
        nodal_forces (np.ndarray): (N_struc, 3)
        ref_point (np.ndarray, optional): (3,) reference for moment. Default: origin.
    """
    if ref_point is None:
        ref_point = np.zeros(3)

    F_aero = np.sum(aero_forces, axis=0)
    F_struc = np.sum(nodal_forces, axis=0)
    M_aero = np.sum(np.cross(aero_points - ref_point, aero_forces), axis=0)
    M_struc = np.sum(np.cross(struc_nodes - ref_point, nodal_forces), axis=0)

    print("=== Force & Moment Conservation Check ===")
    print(f"  Total aero force:    {F_aero}")
    print(f"  Total struct force:  {F_struc}")
    print(f"  Force error:         {np.linalg.norm(F_struc - F_aero):.2e}")
    print(f"  Total aero moment:   {M_aero}")
    print(f"  Total struct moment: {M_struc}")
    print(f"  Moment error:        {np.linalg.norm(M_struc - M_aero):.2e}")
    print("==========================================")


def _load_cp_distribution(cp_path):
    """
    Load Cp distribution data from a file.

    Args:
        cp_path (str or Path): Path to Cp file with columns: x y Cp.

    Returns:
        tuple: (x, y, cp) arrays.
    """
    rows = []
    with open(cp_path, "r") as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if len(parts) < 3:
                continue
            rows.append([float(parts[0]), float(parts[1]), float(parts[2])])
    data = np.array(rows, dtype=float)
    if data.size == 0:
        raise ValueError(f"No Cp data found in {cp_path}")
    return data[:, 0], data[:, 1], data[:, 2]


def _average_duplicate_x(x, cp):
    """
    Average Cp values for duplicate x positions.

    Args:
        x (np.ndarray): Chordwise positions.
        cp (np.ndarray): Cp values at those positions.

    Returns:
        tuple: (uniq_x, cp_mean) with averaged Cp per unique x.
    """
    order = np.argsort(x)
    x_sorted = x[order]
    cp_sorted = cp[order]
    uniq_x, inv = np.unique(x_sorted, return_inverse=True)
    cp_sum = np.zeros_like(uniq_x, dtype=float)
    counts = np.zeros_like(uniq_x, dtype=float)
    np.add.at(cp_sum, inv, cp_sorted)
    np.add.at(counts, inv, 1.0)
    cp_mean = cp_sum / np.maximum(counts, 1.0)
    return uniq_x, cp_mean


def _split_surfaces_by_order(x, y, cp):
    """
    Split Cp data into upper and lower surfaces using file order.

    Assumes the data is ordered around the airfoil (LE -> TE -> LE). The split
    is performed at the first crossing where x >= 1.0 (TE). The first segment
    is treated as upper, the second as lower. Only x-ordering is used.

    Args:
        x (np.ndarray): Chordwise positions.
        y (np.ndarray): Surface-normal positions.
        cp (np.ndarray): Cp values.

    Returns:
        tuple: (x_u, y_u, cp_u), (x_l, y_l, cp_l)
    """
    x = np.asarray(x)
    y = np.asarray(y)
    cp = np.asarray(cp)

    te_candidates = np.where(x >= 1.0)[0]
    if te_candidates.size > 0:
        te_idx = int(te_candidates[0])
    else:
        te_idx = int(np.argmax(x))

    if te_idx == 0 or te_idx >= len(x) - 1:
        mid = len(x) // 2
        return (x[:mid], y[:mid], cp[:mid]), (x[mid:], y[mid:], cp[mid:])

    x_u, y_u, cp_u = x[: te_idx + 1], y[: te_idx + 1], cp[: te_idx + 1]
    x_l, y_l, cp_l = x[te_idx:], y[te_idx:], cp[te_idx:]
    return (x_u, y_u, cp_u), (x_l, y_l, cp_l)


def chordwise_weights_from_cp_file(cp_path, x_targets):
    """
    Compute normalized chordwise weights from a Cp distribution file.

    The weights are based on Delta Cp = Cp_lower - Cp_upper interpolated
    onto x_targets. If Delta Cp is invalid, returns uniform weights.

    Args:
        cp_path (str or Path): Path to Cp file with columns: x y Cp.
        x_targets (np.ndarray): Chordwise positions in [0, 1] to weight.

    Returns:
        np.ndarray: Normalized weights for each x_target.
    """
    x, y, cp = _load_cp_distribution(cp_path)
    (x_u, y_u, cp_u), (x_l, y_l, cp_l) = _split_surfaces_by_order(x, y, cp)

    if len(x_u) == 0 or len(x_l) == 0:
        return np.full_like(x_targets, 1.0 / len(x_targets), dtype=float)

    x_u, cp_u = _average_duplicate_x(x_u, cp_u)
    x_l, cp_l = _average_duplicate_x(x_l, cp_l)

    if len(x_u) < 2 or len(x_l) < 2:
        return np.full_like(x_targets, 1.0 / len(x_targets), dtype=float)

    cp_upper = np.interp(x_targets, x_u, cp_u, left=cp_u[0], right=cp_u[-1])
    cp_lower = np.interp(x_targets, x_l, cp_l, left=cp_l[0], right=cp_l[-1])

    delta_cp = cp_lower - cp_upper
    delta_cp = np.maximum(delta_cp, 0.0)
    total = np.sum(delta_cp)
    if total <= 0.0:
        return np.full_like(x_targets, 1.0 / len(x_targets), dtype=float)
    return delta_cp / total


def plot_delta_cp_and_weights(cp_path, n_bins=10):
    """
    Plot Cp distributions, Delta Cp, and chordwise weights.

    Args:
        cp_path (str or Path): Path to Cp file with columns: x y Cp.
        n_bins (int): Number of chordwise bins/targets for weights.

    Returns:
        None. Displays a 3x1 plot.
    """
    x, y, cp = _load_cp_distribution(cp_path)
    (x_u, y_u, cp_u), (x_l, y_l, cp_l) = _split_surfaces_by_order(x, y, cp)
    if len(x_u) == 0 or len(x_l) == 0:
        raise ValueError("Cp file must contain both upper and lower surface points.")

    x_u, cp_u = _average_duplicate_x(x_u, cp_u)
    x_l, cp_l = _average_duplicate_x(x_l, cp_l)

    x_full = np.union1d(x_u, x_l)
    cp_upper_full = np.interp(x_full, x_u, cp_u, left=cp_u[0], right=cp_u[-1])
    cp_lower_full = np.interp(x_full, x_l, cp_l, left=cp_l[0], right=cp_l[-1])
    delta_cp = cp_lower_full - cp_upper_full

    x_targets = np.linspace(0.0, 1.0, n_bins)
    weights = chordwise_weights_from_cp_file(cp_path, x_targets)

    fig, axes = plt.subplots(3, 1, sharex=True, figsize=(8, 10))
    n_steps = 6

    axes[0].plot(x_u, cp_u, "o-", color="tab:blue", label=r"$C_p$ upper", markersize=2)
    axes[0].plot(
        x_l, cp_l, "o-", color="tab:orange", label=r"$C_p$ lower", markersize=2
    )
    axes[0].plot([0, 1], [0, 0], "--", color="k", lw=1)
    axes[0].set_ylabel(r"$C_p$ (-)")
    axes[0].set_title("Cp Distribution")
    axes[0].invert_yaxis()
    axes[0].grid(True, linestyle="--", linewidth=0.5)
    axes[0].legend()

    axes[1].plot(
        x_full, delta_cp, "o-", color="tab:green", label=r"$\Delta C_p$", markersize=2
    )
    bar_width = 0.8 / max(n_bins - 1, 1)
    axes[1].bar(
        x_targets,
        weights,
        width=bar_width,
        color="tab:gray",
        alpha=0.4,
        label="Chordwise weights",
    )
    axes[1].plot([0, 1], [0, 0], "--", color="k", lw=1)
    axes[1].set_ylabel(r"$\Delta C_p$ / weights (-)")
    axes[1].set_title("Delta Cp and Weights")
    axes[1].grid(True, linestyle="--", linewidth=0.5)
    axes[1].legend()

    axes[2].plot(
        x_targets,
        weights,
        "o-",
        color="tab:gray",
        label="Chordwise weights",
        markersize=2,
    )
    axes[2].set_xlabel(r"$x/c$ (-)")
    axes[2].set_ylabel("Weight (-)")
    axes[2].set_title("Weights")
    axes[2].grid(True, linestyle="--", linewidth=0.5)
    axes[2].legend()

    for ax in axes:
        ax.set_xlim(0, 1)
        ax.set_xticks(np.linspace(0, 1, n_steps))

    plt.tight_layout()
    plt.show()


def main(
    coupling_method: str,
    f_aero_wing_vsm_format: np.ndarray,
    struc_nodes: np.ndarray,
    panel_cp_locations: np.ndarray,
    aero2struc_mapping: np.ndarray,
    is_with_coupling_plot: bool,
    config_aer2struc: dict,
    canopy_sections,
    strut_sections,
    panels,
    section_ids=None,
    cp_distribution_path=None,
    is_with_delta_cp_and_weights_plot=False,
):
    """
    Main interface for mapping aerodynamic panel forces to structural nodes.

    Args:
        coupling_method (str): Coupling method name (e.g., "NN").
        f_aero_wing_vsm_format (np.ndarray): Aerodynamic forces per panel (n_panels,3).
        struc_nodes (np.ndarray): Structural node positions (n_struc,3).
        panel_cp_locations (np.ndarray): Panel control points (n_panels,3).
        aero2struc_mapping (np.ndarray): Mapping from panels to 4 node indices (n_panels,4).
        is_with_coupling_plot (bool): If True, plot the mapping.
        p (float): Power for inverse-distance weighting.
        eps (float): Small value to avoid division by zero.
        cp_distribution_path (str or Path, optional): Path to Cp distribution file
            used for chordwise weighting. Defaults to cp_AOA_8.dat in the data folder.

    Returns:
        np.ndarray: Forces on structural nodes (n_struc,3).
    """
    # Displacing the single spanwise aero force over 10 chordwise nodes
    n_chordwise_nodes = 10
    vsm_wing_nodes_distributed_chordwise = []
    vsm_wing_forces_distributed_chordwise = []

    if cp_distribution_path is None:
        cp_distribution_path = (
            Path(__file__).resolve().parents[2]
            / "data"
            / "TUDELFT_V3_KITE"
            / "cp_distributions"
            / "cp_AOA_8.dat"
        )

    t_vals = np.linspace(0.0, 1.0, n_chordwise_nodes)
    chordwise_weights = chordwise_weights_from_cp_file(cp_distribution_path, t_vals)
    for panel, f_panel in zip(panels, f_aero_wing_vsm_format):
        le_mid = 0.5 * (panel.LE_point_1 + panel.LE_point_2)
        te_mid = 0.5 * (panel.TE_point_1 + panel.TE_point_2)
        vec_chord = te_mid - le_mid

        chord_nodes = le_mid[None, :] + t_vals[:, None] * vec_chord[None, :]
        vsm_wing_nodes_distributed_chordwise.append(chord_nodes)

        f_nodes = chordwise_weights[:, None] * f_panel[None, :]
        vsm_wing_forces_distributed_chordwise.append(f_nodes)

    vsm_wing_nodes_distributed_chordwise = np.vstack(
        vsm_wing_nodes_distributed_chordwise
    )
    vsm_wing_forces_distributed_chordwise = np.vstack(
        vsm_wing_forces_distributed_chordwise
    )

    if is_with_delta_cp_and_weights_plot:
        plot_delta_cp_and_weights(
            Path("data/TUDELFT_V3_KITE/cp_distributions/cp_AOA_8.dat"),
            n_bins=10,
        )

    # mapping the distributed aerodynamic forces to the structural nodes
    sections = build_ordered_sections(struc_nodes, canopy_sections, strut_sections)[0]

    if section_ids is None:
        active_sections = sections
    elif isinstance(section_ids, int):
        active_sections = [sections[section_ids]]
    else:
        active_sections = [sections[sid] for sid in section_ids]

    f_aero_wing = map_aero_forces_to_struct_nodes(
        aero_points=vsm_wing_nodes_distributed_chordwise,
        aero_forces=vsm_wing_forces_distributed_chordwise,
        struc_nodes=struc_nodes,
        sections=active_sections,
    )

    verify_force_moment_conservation(
        aero_points=vsm_wing_nodes_distributed_chordwise,
        aero_forces=vsm_wing_forces_distributed_chordwise,
        struc_nodes=struc_nodes,
        nodal_forces=f_aero_wing,
    )
    if is_with_coupling_plot:
        plot_aerodynamic_forces_chordwise_distributed(
            panel_cps=panel_cp_locations,
            f_aero_chordwise=f_aero_wing_vsm_format,
            vsm_wing_nodes_distributed_chordwise=vsm_wing_nodes_distributed_chordwise,
            vsm_wing_forces_distributed_chordwise=vsm_wing_forces_distributed_chordwise,
            nodes_struc=struc_nodes,
            force_struc=f_aero_wing,
        )
    return f_aero_wing


def initialize_mapping(
    panels: np.ndarray,
    struc_nodes: np.ndarray,
    struc_node_le_indices: np.ndarray,
    struc_node_te_indices: np.ndarray,
) -> np.ndarray:
    """
    For each panel CP, find the two LE and two TE structural‐node indices
    whose y-coordinates bracket the CP’s y. Returns an (n_panels, 4) array
    of [le_lo, le_hi, te_lo, te_hi].

    Args:
        panels (np.ndarray): Array of panel objects with .aerodynamic_center attribute.
        struc_nodes (np.ndarray): Structural node positions (n_nodes,3).
        struc_node_le_indices (np.ndarray): Indices of leading edge nodes.
        struc_node_te_indices (np.ndarray): Indices of trailing edge nodes.

    Returns:
        np.ndarray: Mapping array (n_panels, 4).
    """

    # extract and sort LE candidates by their y
    le_coords = []
    for struc_node_le_idx in struc_node_le_indices:
        le_coords.append(struc_nodes[struc_node_le_idx])

    le_coords = np.array(le_coords)
    le_y = le_coords[:, 1]
    le_order = np.argsort(le_y)
    le_sorted_idx = np.array(struc_node_le_indices)[le_order]
    le_sorted_y = le_y[le_order]

    # same for TE
    te_coords = []
    for struc_node_te_idx in struc_node_te_indices:
        te_coords.append(struc_nodes[struc_node_te_idx])

    te_coords = np.array(te_coords)
    te_y = te_coords[:, 1]
    te_order = np.argsort(te_y)
    te_sorted_idx = np.array(struc_node_te_indices)[te_order]
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


def aero2struc_NN_vsm(
    f_aero_wing_vsm_format: np.ndarray,
    struc_nodes: np.ndarray,
    panel_cps: np.ndarray,
    panel_corner_map: np.ndarray,
    power_for_inverse_weighting: float = 2,
    eps: float = 1e-6,
    is_with_coupling_plot: bool = False,
):
    """
    Distribute each panel's resultant force (at its CoP) onto the four
    structural corner nodes given in panel_corner_map using inverse-distance weighting.

    Args:
        f_aero_wing_vsm_format (np.ndarray): Aerodynamic forces per panel (n_panels,3).
        struc_nodes (np.ndarray): Structural node positions (n_struc,3).
        panel_cps (np.ndarray): Panel control points (n_panels,3).
        panel_corner_map (np.ndarray): Mapping from panels to 4 node indices (n_panels,4).
        p (float): Power for inverse-distance weighting.
        eps (float): Small value to avoid division by zero.
        is_with_coupling_plot (bool): If True, plot the mapping.

    Returns:
        np.ndarray: Forces on structural nodes (n_struc,3).
    """

    n_struc = len(struc_nodes)
    f_aero_wing = np.zeros((n_struc, 3), dtype=float)

    for i, (cp, frc) in enumerate(zip(panel_cps, f_aero_wing_vsm_format)):
        sel_idx = panel_corner_map[i]  # [le_lo, le_hi, te_lo, te_hi]
        sel_coords = struc_nodes[sel_idx]  # (4,3)

        # true inverse-distance weighting across the 4 nodes
        d = np.linalg.norm(sel_coords - cp[None, :], axis=1)
        w = 1.0 / (d**power_for_inverse_weighting + eps)
        w /= np.sum(w)

        f_vals = w[:, None] * frc[None, :]  # (4,3)

        # accumulate
        for local_j, glob_j in enumerate(sel_idx):
            f_aero_wing[glob_j] += f_vals[local_j]

    if is_with_coupling_plot:
        plot_aerodynamic_forces_chordwise_distributed(
            panel_cps=panel_cps,
            f_aero_chordwise=f_aero_wing_vsm_format,
            nodes_struc=struc_nodes,
            force_struc=f_aero_wing,
        )

    return f_aero_wing


# def main(
#     coupling_method: str,
#     f_aero_wing_vsm_format: np.ndarray,
#     struc_nodes: np.ndarray,
#     panel_cp_locations: np.ndarray,
#     aero2struc_mapping: np.ndarray,
#     is_with_coupling_plot: bool,
#     config_aer2struc: dict,
# ):
#     """
#     Main interface for mapping aerodynamic panel forces to structural nodes.

#     Args:
#         coupling_method (str): Coupling method name (e.g., "NN").
#         f_aero_wing_vsm_format (np.ndarray): Aerodynamic forces per panel (n_panels,3).
#         struc_nodes (np.ndarray): Structural node positions (n_struc,3).
#         panel_cp_locations (np.ndarray): Panel control points (n_panels,3).
#         aero2struc_mapping (np.ndarray): Mapping from panels to 4 node indices (n_panels,4).
#         is_with_coupling_plot (bool): If True, plot the mapping.
#         p (float): Power for inverse-distance weighting.
#         eps (float): Small value to avoid division by zero.

#     Returns:
#         np.ndarray: Forces on structural nodes (n_struc,3).
#     """

#     if coupling_method == config_aer2struc["coupling_method"]:
#         return aero2struc_NN_vsm(
#             f_aero_wing_vsm_format,  # (n_panels,3)
#             struc_nodes,  # (n_struc,3)
#             panel_cp_locations,  # (n_panels,3)
#             aero2struc_mapping,  # (n_panels,4)
#             power_for_inverse_weighting=config_aer2struc["power_for_inverse_weighting"],
#             eps=config_aer2struc["eps"],
#             is_with_coupling_plot=is_with_coupling_plot,
#         )
#     else:
#         raise ValueError("Coupling method not recognized; wrong name or typo")
