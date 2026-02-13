import numpy as np
import matplotlib.pyplot as plt


# TODO: this should be placed in a more general plotting place/module
def plot_aerodynamic_forces_chordwise_distributed(
    panel_cps,
    f_aero_chordwise,
    nodes_struc,
    force_struc=None,
):
    """
    Plot aerodynamic forces distributed chordwise and mapped to structural nodes.

    Args:
        panel_cps (np.ndarray): panel cps (n,3).
        f_aero_chordwise (np.ndarray): Chordwise aerodynamic forces (n,3).
        nodes_struc (np.ndarray): Structural node positions (n_nodes,3).
        force_struc (np.ndarray, optional): Forces on structural nodes (n_nodes,3).

    Returns:
        None. Displays a 3D plot.
    """

    # Create a new figure and set up 3D axes
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    # Scatter plot of chordwise points (blue)
    ax.scatter(
        panel_cps[:, 0],
        panel_cps[:, 1],
        panel_cps[:, 2],
        color="black",
        label="Panel center of pressure",
    )

    # Quiver plot for the forces (red arrows)
    ax.quiver(
        panel_cps[:, 0],
        panel_cps[:, 1],
        panel_cps[:, 2],
        f_aero_chordwise[:, 0],
        f_aero_chordwise[:, 1],
        f_aero_chordwise[:, 2],
        # length=1,
        # normalize=True,
        length=0.01,
        color="black",
        label="Panel force vector",
    )

    if force_struc is None:
        # Scatter plot of structural nodes (wing segment corners)
        ax.scatter(
            nodes_struc[:, 0],
            nodes_struc[:, 1],
            nodes_struc[:, 2],
            color="black",
            label="Wing Segment Corners",
        )

        # Annotate each point with its index
        for idx, point in enumerate(nodes_struc):
            ax.text(point[0], point[1], point[2], f"{idx}", color="black")
    else:
        # Scatter plot of structural nodes (wing segment corners) (green)
        ax.scatter(
            nodes_struc[:, 0],
            nodes_struc[:, 1],
            nodes_struc[:, 2],
            color="blue",
            label="Structural nodes",
        )

        # Quiver plot for the forces on structural nodes (yellow arrows)
        ax.quiver(
            nodes_struc[:, 0],
            nodes_struc[:, 1],
            nodes_struc[:, 2],
            force_struc[:, 0],
            force_struc[:, 1],
            force_struc[:, 2],
            # length=1,
            length=0.01,
            # normalize=True,
            color="red",
            label="Mapped aerodynamic force vector onto structural nodes",
        )

    # Set equal scale for all axes
    points_all = np.concatenate((panel_cps, nodes_struc), axis=0)
    bb = points_all.max(axis=0) - points_all.min(axis=0)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_box_aspect(bb)
    ax.set_title("Aerodynamic Forces and Structural Nodes")
    ax.legend()
    plt.show()


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
