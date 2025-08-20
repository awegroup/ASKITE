import numpy as np


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
    return {
        "positions": np.zeros((nt, n_pts, 3)),
        "f_ext": np.zeros((nt, n_pts, 3)),
        "f_int": np.zeros((nt, n_pts, 3)),
        "residual_norm": np.zeros(nt),
        "max_residual": np.zeros(nt),
    }


def update_tracking_arrays(
    tracking_data, idx, psystem, struc_nodes, nodes_struc_prev, f_ext_flat, f_int_flat
):
    """
    Update tracking arrays with simulation results for a single time step.

    Args:
        tracking_data (dict): Tracking arrays to update.
        idx (int): Current time step index.
        psystem (ParticleSystem): Particle system object.
        struc_nodes (np.ndarray): Current node positions (n_nodes, 3).
        nodes_struc_prev (np.ndarray or None): Previous node positions (n_nodes, 3).
        f_ext_flat (np.ndarray): Flattened external force vector (n_nodes*3,).
        f_int_flat (np.ndarray): Flattened internal force vector (n_nodes*3,).

    Returns:
        None. Updates tracking_data in place.
    """
    # Unpack 3D storage
    pos3d = tracking_data["positions"]
    ext3d = tracking_data["f_ext"]
    int3d = tracking_data["f_int"]

    n_pts = pos3d.shape[1]

    # 1) Positions
    pos_flat, _ = psystem.x_v_current  # shape (n_pts*3,)
    pos3d[idx] = pos_flat.reshape(n_pts, 3)

    # 2) External & internal forces: reshape before storing
    ext3d[idx] = f_ext_flat.reshape(n_pts, 3)
    int3d[idx] = f_int_flat.reshape(n_pts, 3)

    # 3) Norms
    tracking_data["residual_norm"][idx] = np.linalg.norm(f_int_flat)
    tracking_data["max_residual"][idx] = np.max(np.abs(f_int_flat))
