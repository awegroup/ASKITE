import numpy as np
import matplotlib.pyplot as plt
from kitesim.plotting import plot_aerodynamic_forces_chordwise_distributed


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
