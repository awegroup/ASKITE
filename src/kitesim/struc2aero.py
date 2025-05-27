import numpy as np


def main(
    struc_nodes,
    n_wing_nodes,
    struc_le_idx_list,
    struc_te_idx_list,
    n_aero_panels_per_struc_section,
):
    """
    Generate arrays of leading and trailing edge points for the aerodynamic solver,
    interpolating if needed.

    Args:
        struc_nodes (np.ndarray): Structural node positions (n_nodes, 3).
        struc_le_idx_list (list): Indices of leading edge nodes.
        struc_te_idx_list (list): Indices of trailing edge nodes.
        n_aero_panels_per_struc_section (int): Number of aerodynamic panels per structural section.

    Returns:
        le_arr (np.ndarray): Interpolated leading edge points.
        te_arr (np.ndarray): Interpolated trailing edge points.
    """

    # Correct for smaller chord length structural section than aerodynamic section
    for idx, _ in enumerate(struc_nodes[1:n_wing_nodes]):
        # if LE
        if i in struc_le_idx_list:
            print(f"")
        # if TE
        elif i in struc_te_idx_list:
            print(f"")

    le_arr = struc_nodes[struc_le_idx_list]
    te_arr = struc_nodes[struc_te_idx_list]

    if n_aero_panels_per_struc_section == 0:
        le_arr = struc_nodes[struc_le_idx_list]
        te_arr = struc_nodes[struc_te_idx_list]
    else:

        def interpolate_edges(
            struc_nodes,
            struc_le_idx_list,
            struc_te_idx_list,
            n_aero_panels_per_struc_section,
        ):
            """
            Interpolate leading and trailing edge points to create a denser array of points.

            Args:
                struc_nodes (np.ndarray): Structural node positions (n_nodes, 3).
                struc_le_idx_list (list): Indices of leading edge nodes.
                struc_te_idx_list (list): Indices of trailing edge nodes.
                n_aero_panels_per_struc_section (int): Number of aerodynamic panels per structural section.

            Returns:
                le_dense (np.ndarray): Densely sampled leading edge points.
                te_dense (np.ndarray): Densely sampled trailing edge points.
            """

            # Extract original LE and TE arrays
            le_arr = np.array([struc_nodes[i] for i in struc_le_idx_list])
            te_arr = np.array([struc_nodes[i] for i in struc_te_idx_list])

            def interpolate_points(arr):
                new_arr = []
                for i in range(len(arr) - 1):
                    p0, p1 = arr[i], arr[i + 1]
                    for j in range(n_aero_panels_per_struc_section):
                        t = j / n_aero_panels_per_struc_section
                        new_point = (1 - t) * p0 + t * p1
                        new_arr.append(new_point)
                new_arr.append(arr[-1])  # Add the final point
                return np.array(new_arr)

            le_dense = interpolate_points(le_arr)
            te_dense = interpolate_points(te_arr)

            return le_dense, te_dense

        le_arr, te_arr = interpolate_edges(
            struc_nodes,
            struc_le_idx_list,
            struc_te_idx_list,
            n_aero_panels_per_struc_section,
        )
    return le_arr, te_arr
