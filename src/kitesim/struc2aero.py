import numpy as np
import copy


# def interpolate_edges(
#     struc_nodes,
#     struc_le_idx_list,
#     struc_te_idx_list,
#     n_aero_panels_per_struc_section,
# ):
#     """
#     Interpolate leading and trailing edge points to create a denser array of points.

#     Args:
#         struc_nodes (np.ndarray): Structural node positions (n_nodes, 3).
#         struc_le_idx_list (list): Indices of leading edge nodes.
#         struc_te_idx_list (list): Indices of trailing edge nodes.
#         n_aero_panels_per_struc_section (int): Number of aerodynamic panels per structural section.

#     Returns:
#         le_dense (np.ndarray): Densely sampled leading edge points.
#         te_dense (np.ndarray): Densely sampled trailing edge points.
#     """

#     # Extract original LE and TE arrays
#     le_arr = np.array([struc_nodes[i] for i in struc_le_idx_list])
#     te_arr = np.array([struc_nodes[i] for i in struc_te_idx_list])

#     le_dense = interpolate_points(le_arr)
#     te_dense = interpolate_points(te_arr)

#     return le_dense, te_dense


def interpolate_points(arr, n_aero_panels_per_struc_section):
    new_arr = []
    for i in range(len(arr) - 1):
        p0, p1 = arr[i], arr[i + 1]
        for j in range(n_aero_panels_per_struc_section):
            t = j / n_aero_panels_per_struc_section
            new_point = (1 - t) * p0 + t * p1
            new_arr.append(new_point)
    new_arr.append(arr[-1])  # Add the final point
    return np.array(new_arr)


def main(
    struc_nodes,
    struc_node_le_indices,
    struc_node_te_indices,
    n_aero_panels_per_struc_section,
):
    """
    Generate arrays of leading and trailing edge points for the aerodynamic solver,
    interpolating if needed.

    Args:
        struc_nodes (np.ndarray): Structural node positions (n_nodes, 3).
        struc_node_le_indices (list): Indices of leading edge nodes.
        struc_node_te_indices (list): Indices of trailing edge nodes.
        n_aero_panels_per_struc_section (int): Number of aerodynamic panels per structural section.

    Returns:
        le_arr (np.ndarray): Interpolated leading edge points.
        te_arr (np.ndarray): Interpolated trailing edge points.
    """

    le_arr = np.array([struc_nodes[i] for i in struc_node_le_indices])
    te_arr = np.array([struc_nodes[i] for i in struc_node_te_indices])

    # TODO: remove hardcoded values
    ## Correct leading edge and trailing edge points for the full wing
    delta_te_arr_LAP_to_full_wing = [
        0.04,
        0.23,
        0.24,
        0.26,
        0.28,
        0.28,
        0.26,
        0.24,
        0.23,
        0.04,
    ]
    for idx in range(len(le_arr)):
        # Compute direction vector from LE to TE
        vec_le_te = te_arr[idx] - le_arr[idx]
        vec_le_te_norm = vec_le_te / np.linalg.norm(vec_le_te)
        # Apply offset along this direction to TE
        te_arr[idx] += delta_te_arr_LAP_to_full_wing[idx] * vec_le_te_norm
        le_arr[idx] -= 0.04 * vec_le_te_norm

    if n_aero_panels_per_struc_section == 1:
        return le_arr, te_arr
    elif n_aero_panels_per_struc_section == 0:
        raise ValueError("n_aero_panels_per_struc_section must be greater than 0.")
    else:
        return interpolate_points(
            le_arr, n_aero_panels_per_struc_section
        ), interpolate_points(te_arr, n_aero_panels_per_struc_section)
