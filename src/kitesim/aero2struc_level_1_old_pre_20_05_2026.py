import logging
import numpy as np
import matplotlib.pyplot as plt
from kitesim.plotting import plot_aerodynamic_forces_chordwise_distributed


def check_moment_preservation(
    f_aero_panel: np.ndarray,
    panel_cps: np.ndarray,
    f_aero_mapped: np.ndarray,
    struc_nodes: np.ndarray,
    ref_point: np.ndarray = None,
) -> dict:
    """
    Check whether the aero→struc force mapping preserves total force and moment.

    Computes:
        M_aero  = Σ_i  (r_cp_i − r_ref) × F_panel_i       (original panel forces)
        M_struc = Σ_j  (r_node_j − r_ref) × F_mapped_j     (mapped nodal forces)

    A pure inverse-distance-weight mapping preserves force (Σw=1) but generally
    does NOT preserve moment, because the weighted centroid of the corner nodes
    differs from the panel CP.

    Args:
        f_aero_panel:  (n_panels, 3) panel forces at their CPs.
        panel_cps:     (n_panels, 3) panel control-point locations.
        f_aero_mapped: (n_struc, 3)  mapped forces on structural nodes.
        struc_nodes:   (n_struc, 3)  structural node positions.
        ref_point:     (3,)  reference point for moment calculation.
                       Defaults to the origin [0, 0, 0].

    Returns:
        dict with keys:
            F_aero_total:   (3,) total aero force from panels
            F_struc_total:  (3,) total mapped force on struct nodes
            dF:             (3,) force error  (should be ≈0)
            dF_norm:        float  ||dF||
            M_aero:         (3,) total moment from panel forces
            M_struc:        (3,) total moment from mapped nodal forces
            dM:             (3,) moment error
            dM_norm:        float  ||dM||
            dM_rel:         float  ||dM|| / ||M_aero||  (relative moment error)
    """
    if ref_point is None:
        ref_point = np.zeros(3)

    # --- total force ---
    F_aero = np.sum(f_aero_panel, axis=0)
    F_struc = np.sum(f_aero_mapped, axis=0)
    dF = F_struc - F_aero

    # --- total moment about ref_point ---
    M_aero = np.zeros(3)
    for cp, frc in zip(panel_cps, f_aero_panel):
        M_aero += np.cross(cp - ref_point, frc)

    M_struc = np.zeros(3)
    for node, frc in zip(struc_nodes, f_aero_mapped):
        M_struc += np.cross(node - ref_point, frc)

    dM = M_struc - M_aero
    M_aero_norm = np.linalg.norm(M_aero)
    dM_rel = np.linalg.norm(dM) / M_aero_norm if M_aero_norm > 1e-12 else 0.0

    result = {
        "F_aero_total": F_aero,
        "F_struc_total": F_struc,
        "dF": dF,
        "dF_norm": np.linalg.norm(dF),
        "M_aero": M_aero,
        "M_struc": M_struc,
        "dM": dM,
        "dM_norm": np.linalg.norm(dM),
        "dM_rel": dM_rel,
    }

    logging.info(
        f"Moment preservation check (ref={ref_point}):\n"
        f"  Force error  ||dF|| = {result['dF_norm']:.6e} N\n"
        f"  Moment aero  ||M||  = {M_aero_norm:.3f} Nm\n"
        f"  Moment error ||dM|| = {result['dM_norm']:.6e} Nm  "
        f"(relative: {result['dM_rel']:.4%})\n"
        f"  dM components = [{dM[0]:.4f}, {dM[1]:.4f}, {dM[2]:.4f}] Nm"
    )

    return result


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
    Distribute each panel's resultant force (at its CP) onto the four
    structural corner nodes using bilinear interpolation.

    For each panel with corners [le_lo, le_hi, te_lo, te_hi]:
      1. Spanwise: eta from LE y-coordinates
         eta = (y_cp - y_le_lo) / (y_le_hi - y_le_lo), clamped to [0, 1].
      2. Chordwise: project CP onto the LE→TE segment for each span side
         xi_lo = dot(cp - le_lo, te_lo - le_lo) / |te_lo - le_lo|²
         xi_hi = dot(cp - le_hi, te_hi - le_hi) / |te_hi - le_hi|²
      3. Bilinear weights:
         w_le_lo = (1 - eta) * (1 - xi_lo)
         w_te_lo = (1 - eta) * xi_lo
         w_le_hi = eta * (1 - xi_hi)
         w_te_hi = eta * xi_hi

    Force is exactly preserved (weights sum to 1).
    Moment error is proportional to the out-of-plane offset between the
    CP and the bilinear surface of the 4 corner nodes (typically very small).

    Args:
        f_aero_wing_vsm_format (np.ndarray): Aerodynamic forces per panel (n_panels,3).
        struc_nodes (np.ndarray): Structural node positions (n_struc,3).
        panel_cps (np.ndarray): Panel control points (n_panels,3).
        panel_corner_map (np.ndarray): Mapping from panels to 4 node indices (n_panels,4).
        power_for_inverse_weighting (float): Unused, kept for API compatibility.
        eps (float): Small value to avoid division by zero.
        is_with_coupling_plot (bool): If True, plot the mapping.

    Returns:
        np.ndarray: Forces on structural nodes (n_struc,3).
    """

    n_struc = len(struc_nodes)
    f_aero_wing = np.zeros((n_struc, 3), dtype=float)

    for i, (cp, frc) in enumerate(zip(panel_cps, f_aero_wing_vsm_format)):
        le_lo, le_hi, te_lo, te_hi = panel_corner_map[i]

        r_le_lo = struc_nodes[le_lo]
        r_le_hi = struc_nodes[le_hi]
        r_te_lo = struc_nodes[te_lo]
        r_te_hi = struc_nodes[te_hi]

        # --- spanwise weight eta (from LE y-coordinates) ---
        dy_le = r_le_hi[1] - r_le_lo[1]
        if abs(dy_le) < eps:
            eta = 0.5
        else:
            eta = (cp[1] - r_le_lo[1]) / dy_le
        eta = np.clip(eta, 0.0, 1.0)

        # --- chordwise weight xi: project CP onto LE→TE for each span side ---
        # Low-span side (le_lo → te_lo)
        chord_lo = r_te_lo - r_le_lo
        chord_lo_sq = np.dot(chord_lo, chord_lo)
        if chord_lo_sq < eps * eps:
            xi_lo = 0.0
        else:
            xi_lo = np.dot(cp - r_le_lo, chord_lo) / chord_lo_sq
        xi_lo = np.clip(xi_lo, 0.0, 1.0)

        # High-span side (le_hi → te_hi)
        chord_hi = r_te_hi - r_le_hi
        chord_hi_sq = np.dot(chord_hi, chord_hi)
        if chord_hi_sq < eps * eps:
            xi_hi = 0.0
        else:
            xi_hi = np.dot(cp - r_le_hi, chord_hi) / chord_hi_sq
        xi_hi = np.clip(xi_hi, 0.0, 1.0)

        # --- bilinear weights (sum to 1 by construction) ---
        w_le_lo = (1.0 - eta) * (1.0 - xi_lo)
        w_te_lo = (1.0 - eta) * xi_lo
        w_le_hi = eta * (1.0 - xi_hi)
        w_te_hi = eta * xi_hi

        # accumulate
        f_aero_wing[le_lo] += w_le_lo * frc
        f_aero_wing[le_hi] += w_le_hi * frc
        f_aero_wing[te_lo] += w_te_lo * frc
        f_aero_wing[te_hi] += w_te_hi * frc

    if is_with_coupling_plot:
        plot_aerodynamic_forces_chordwise_distributed(
            panel_cps=panel_cps,
            f_aero_chordwise=f_aero_wing_vsm_format,
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


def main(
    coupling_method: str,
    f_aero_wing_vsm_format: np.ndarray,
    struc_nodes: np.ndarray,
    panel_cp_locations: np.ndarray,
    aero2struc_mapping: np.ndarray,
    is_with_coupling_plot: bool,
    config_aer2struc: dict,
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

    Returns:
        np.ndarray: Forces on structural nodes (n_struc,3).
    """

    if coupling_method == config_aer2struc["coupling_method"]:
        return aero2struc_NN_vsm(
            f_aero_wing_vsm_format,  # (n_panels,3)
            struc_nodes,  # (n_struc,3)
            panel_cp_locations,  # (n_panels,3)
            aero2struc_mapping,  # (n_panels,4)
            power_for_inverse_weighting=config_aer2struc["power_for_inverse_weighting"],
            eps=config_aer2struc["eps"],
            is_with_coupling_plot=is_with_coupling_plot,
        )
    else:
        raise ValueError("Coupling method not recognized; wrong name or typo")
