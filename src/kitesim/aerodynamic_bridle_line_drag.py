import numpy as np


def compute_line_aerodynamic_force(p1, p2, d, va, cd_cable, cf_cable, rho):

    if p1[2] > p2[2]:
        p1, p2 = p2, p1

    length = np.linalg.norm(p2 - p1)
    ej = (p2 - p1) / length
    theta = np.arccos(np.dot(va, ej) / (np.linalg.norm(va) * np.linalg.norm(ej)))

    cd_t = cd_cable * np.sin(theta) ** 3 + np.pi * cf_cable * np.cos(theta) ** 3
    cl_t = (
        cd_cable * np.sin(theta) ** 2 * np.cos(theta)
        - np.pi * cf_cable * np.sin(theta) * np.cos(theta) ** 2
    )
    dir_D = va / np.linalg.norm(va)  # Drag direction
    dir_L = -(ej - np.dot(ej, dir_D) * dir_D)  # Lift direction
    dynamic_pressure_area = 0.5 * rho * np.linalg.norm(va) ** 2 * length * d

    # Calculate lift and drag using the common factor
    lift_j = dynamic_pressure_area * cl_t * dir_L
    drag_j = dynamic_pressure_area * cd_t * dir_D

    return lift_j + drag_j


def main(
    struc_nodes,
    bridle_connectivity_arr,
    bridle_diameters_arr,
    vel_app,
    rho,
    cd_cable,
    cf_cable,
):
    """
    Compute aerodynamic forces on all bridle lines and distribute them to nodes.

    Args:
        struc_nodes: Array of structural node positions, shape (n_nodes, 3)
        bridle_connectivity_arr: Array of bridle line connections, shape (n_bridle_lines, 2)
                                 Each row contains [node_i, node_j] indices
        bridle_diameters_arr: Array of bridle line diameters, shape (n_bridle_lines,)
        vel_app: Apparent wind velocity vector, shape (3,)
        rho: Air density [kg/m³]
        cd_cable: Drag coefficient for cables
        cf_cable: Friction coefficient for cables

    Returns:
        f_aero_bridle: Array of aerodynamic forces on nodes, shape (n_nodes, 3)
    """
    # Initialize force array with zeros for all structural nodes
    f_aero_bridle = np.zeros_like(struc_nodes)

    # Loop through each bridle line segment
    for idx, (ci, cj) in enumerate(bridle_connectivity_arr):
        # Get node positions
        p1 = struc_nodes[ci]
        p2 = struc_nodes[cj]

        # Get diameter for this line segment
        d = bridle_diameters_arr[idx]

        # Compute total aerodynamic force on this line segment
        f_line_total = compute_line_aerodynamic_force(
            p1, p2, d, vel_app, cd_cable, cf_cable, rho
        )

        # Distribute force equally to both nodes (50/50 split)
        f_aero_bridle[ci] += f_line_total / 2.0
        f_aero_bridle[cj] += f_line_total / 2.0

    return f_aero_bridle
