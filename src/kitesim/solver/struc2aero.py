import numpy as np


def main(
    struc_nodes, struc_le_idx_list, struc_te_idx_list, n_aero_panels_per_struc_section
):

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
