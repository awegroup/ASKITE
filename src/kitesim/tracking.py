import numpy as np


SCALAR_TRACKING_FIELDS = [
    "residual_norm",
    "max_residual",
    "geometry_update_norm",
    "geometry_update_rel_norm",
    "aero_force_update_norm",
    "lift_wing",
    "drag_wing",
    "side_force_wing",
    "C_L_wing",
    "C_D_wing",
    "C_R_wing",
    "glide_ratio_wing",
    "C_L_total_aero",
    "C_D_total_aero",
    "C_R_total_aero",
    "glide_ratio_total_aero",
    "V_a",
    "omega_aitken",
    "regularization_phase",
    "is_structural_converged",
    "is_aero_converged",
    "qs_success",
    "projected_span_m",
    "projected_area_m2",
    "midspan_pitch_deg",
    "mean_twist_deg",
    "max_abs_twist_deg",
]

VECTOR_TRACKING_FIELDS = [
    "f_residual",
    "aero_force_wing_total",
    "aero_force_bridle_total",
    "aero_force_total",
]


def setup_tracking_arrays(n_pts, t_vector):
    """
    Initialize tracking arrays for simulation results.

    Args:
        n_pts (int): Number of nodes/particles.
        t_vector (np.ndarray): Array of time steps.

    Returns:
        dict: Dictionary with preallocated arrays for positions, forces, and tracking metrics.
    """
    nt = len(t_vector)
    tracking_data = {
        "positions": np.zeros((nt, n_pts, 3)),
        "f_ext": np.zeros((nt, n_pts, 3)),
        "f_int": np.zeros((nt, n_pts, 3)),
    }

    for field in VECTOR_TRACKING_FIELDS:
        if field == "f_residual":
            tracking_data[field] = np.zeros((nt, n_pts, 3))
        else:
            tracking_data[field] = np.full((nt, 3), np.nan, dtype=float)

    for field in SCALAR_TRACKING_FIELDS:
        tracking_data[field] = np.full(nt, np.nan, dtype=float)

    return tracking_data


def update_tracking_arrays(
    tracking_data,
    idx,
    struc_nodes,
    f_ext_flat=None,
    f_int_flat=None,
    f_residual_flat=None,
    **metrics,
):
    """
    Update tracking arrays with simulation results for a single time step.

    Args:
        tracking_data (dict): Tracking arrays to update.
        idx (int): Current time step index.
        pos3d (np.ndarray): Current 3D positions (n_nodes, 3).
        f_ext_flat (np.ndarray): Flattened external force vector (n_nodes*3,).
        f_int_flat (np.ndarray): Flattened internal force vector (n_nodes*3,).
        f_residual_flat (np.ndarray): Flattened force residual vector (n_nodes*3,).
        **metrics: Optional scalar/vector tracking fields.

    Returns:
        None. Updates tracking_data in place.
    """
    n_pts = tracking_data["positions"].shape[1]

    if f_ext_flat is None:
        f_ext_flat = np.zeros(n_pts * 3)
    if f_int_flat is None:
        f_int_flat = np.zeros(n_pts * 3)
    if f_residual_flat is None:
        # Backward compatibility: older call sites passed the residual as
        # f_int_flat. New solver call sites pass true f_int and f_residual separately.
        f_residual_flat = f_int_flat

    # 1) Positions
    tracking_data["positions"][idx] = struc_nodes

    # 2) External & internal forces: reshape before storing
    tracking_data["f_ext"][idx] = np.asarray(f_ext_flat).reshape(n_pts, 3)
    tracking_data["f_int"][idx] = np.asarray(f_int_flat).reshape(n_pts, 3)
    tracking_data["f_residual"][idx] = np.asarray(f_residual_flat).reshape(n_pts, 3)

    # 3) Norms
    residual = np.asarray(f_residual_flat, dtype=float).reshape(-1)
    tracking_data["residual_norm"][idx] = np.linalg.norm(residual)
    tracking_data["max_residual"][idx] = np.max(np.abs(residual))

    for name, value in metrics.items():
        if name not in tracking_data:
            continue
        arr = tracking_data[name]
        if arr.ndim == 1:
            tracking_data[name][idx] = float(value)
        else:
            tracking_data[name][idx] = np.asarray(value, dtype=float).reshape(
                arr.shape[1:]
            )
