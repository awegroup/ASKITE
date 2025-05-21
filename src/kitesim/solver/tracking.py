import numpy as np


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
    tracking_data, idx, psystem, points, points_prev, f_external_flat, f_residual_flat
):
    # Unpack 3D storage
    pos3d = tracking_data["positions"]
    ext3d = tracking_data["f_external"]
    res3d = tracking_data["f_residual"]

    n_pts = pos3d.shape[1]

    # 1) Positions
    pos_flat, _ = psystem.x_v_current  # shape (n_pts*3,)
    pos3d[idx] = pos_flat.reshape(n_pts, 3)

    # 2) External & residual forces: reshape before storing
    ext3d[idx] = f_external_flat.reshape(n_pts, 3)
    res3d[idx] = f_residual_flat.reshape(n_pts, 3)

    # 3) Norms
    tracking_data["residual_norm"][idx] = np.linalg.norm(f_residual_flat)
    tracking_data["max_residual"][idx] = np.max(np.abs(f_residual_flat))

    # 4) Deltas
    if points_prev is not None:
        dp = points - points_prev
        tracking_data["pos_change"][idx] = np.linalg.norm(dp)

        if hasattr(psystem, "v_current_2D") and hasattr(psystem, "v_prev_2D"):
            dv = psystem.v_current_2D - psystem.v_prev_2D
            tracking_data["vel_change"][idx] = np.linalg.norm(dv)
