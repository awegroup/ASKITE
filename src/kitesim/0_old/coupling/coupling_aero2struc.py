# %%
import numpy as np
from sklearn.neighbors import NearestNeighbors
import logging


def force2nodes(F, Fpoint, nodes, tangential):

    P1 = line_intersect(nodes[0, :], nodes[1, :], Fpoint, Fpoint + tangential)
    d1 = Fpoint - P1

    M1 = np.cross(d1, F)

    P2 = line_intersect(nodes[3, :], nodes[2, :], Fpoint, Fpoint + tangential)

    d2 = P2 - P1
    Fp2 = np.cross(M1, d2) / np.linalg.norm(d2) ** 2

    Fp1 = F - Fp2

    M3 = np.cross(P1 - nodes[0, :], Fp1)
    d3 = nodes[1, :] - nodes[0, :]
    F3 = np.cross(M3, d3) / np.linalg.norm(d3) ** 2

    node1 = Fp1 - F3
    node2 = F3

    M4 = np.cross(P2 - nodes[2, :], Fp2)
    d4 = nodes[3, :] - nodes[2, :]
    F4 = np.cross(M4, d4) / np.linalg.norm(d4) ** 2

    node4 = F4
    node3 = Fp2 - F4

    Fnode = np.array([node1, node2, node3, node4])

    return Fnode


def moment2nodes(M, Mpoint, nodes, tangential):
    d = tangential * 0.05

    dF = np.cross(M, d)
    dF = dF / np.linalg.norm(dF)
    Fmag = np.linalg.norm(M) / np.linalg.norm(np.cross(dF, d))
    F = dF * Fmag

    P1 = Mpoint + d

    Fnode1 = force2nodes(F, P1, nodes, tangential)
    Fnode2 = force2nodes(-F, Mpoint, nodes, tangential)

    Fnode = np.array(Fnode1) + np.array(Fnode2)

    return Fnode


def line_intersect(p1, p2, p3, p4):

    p13 = np.empty(3)
    p43 = np.empty(3)
    p21 = np.empty(3)
    pa = np.empty(3)
    pb = np.empty(3)

    p13[0] = p1[0] - p3[0]
    p13[1] = p1[1] - p3[1]
    p13[2] = p1[2] - p3[2]

    p43[0] = p4[0] - p3[0]
    p43[1] = p4[1] - p3[1]
    p43[2] = p4[2] - p3[2]

    p21[0] = p2[0] - p1[0]
    p21[1] = p2[1] - p1[1]
    p21[2] = p2[2] - p1[2]

    d1343 = p13[0] * p43[0] + p13[1] * p43[1] + p13[2] * p43[2]
    d4321 = p43[0] * p21[0] + p43[1] * p21[1] + p43[2] * p21[2]
    d1321 = p13[0] * p21[0] + p13[1] * p21[1] + p13[2] * p21[2]
    d4343 = p43[0] * p43[0] + p43[1] * p43[1] + p43[2] * p43[2]
    d2121 = p21[0] * p21[0] + p21[1] * p21[1] + p21[2] * p21[2]

    denom = d2121 * d4343 - d4321 * d4321

    numer = d1343 * d4321 - d1321 * d4343

    mua = numer / denom
    mub = (d1343 + d4321 * mua) / d4343

    pa[0] = p1[0] + mua * p21[0]
    pa[1] = p1[1] + mua * p21[1]
    pa[2] = p1[2] + mua * p21[2]
    pb[0] = p3[0] + mub * p43[0]
    pb[1] = p3[1] + mub * p43[1]
    pb[2] = p3[2] + mub * p43[2]

    return pa


def aero2struc(pos, ci, cj, plates, F, M, ringvec, controlpoints):

    lift_force = np.zeros(pos.shape)
    N_struct = len(plates)
    N_split = int(len(controlpoints) / N_struct)
    for i in np.arange(0, len(controlpoints)):  # looping through each panel
        sec = (N_struct - 1) - int((i + 1) / N_split - 0.01)
        # Fi = (F[i][0] +F[i][1])*np.linalg.norm(ringvec[i]['r0'])
        # Fi = (F[i])*np.linalg.norm(ringvec[i]['r0'])
        # Mi = M[i]*np.linalg.norm(ringvec[i]['r0'])
        Fi = F[i]
        Mi = M[i]
        Mi = Mi * controlpoints[i]["airf_coord"][:, 2]

        if sec > 4:
            Pnodes = np.array(
                [
                    pos[plates[sec][0], :],
                    pos[plates[sec][1], :],
                    pos[plates[sec][2], :],
                    pos[plates[sec][3], :],
                ]
            )
        else:
            Pnodes = np.array(
                [
                    pos[plates[sec][1], :],
                    pos[plates[sec][0], :],
                    pos[plates[sec][3], :],
                    pos[plates[sec][2], :],
                ]
            )
        Fnode = force2nodes(
            Fi,
            controlpoints[i]["coordinates_aoa"],
            Pnodes,
            controlpoints[i]["tangential"],
        )
        # print(sum(Fnode)-Fi)
        # M1 = np.cross(Pnodes[0]-controlpoints[i]['coordinates_aoa'], Fnode[0,:])
        # M2 = np.cross(Pnodes[1]-controlpoints[i]['coordinates_aoa'], Fnode[1,:])
        # M3 = np.cross(Pnodes[2]-controlpoints[i]['coordinates_aoa'], Fnode[2,:])
        # M4 = np.cross(Pnodes[3]-controlpoints[i]['coordinates_aoa'], Fnode[3,:])
        # MT = M1+M2+M3+M4
        # print(MT)
        Fnode += moment2nodes(
            Mi,
            controlpoints[i]["coordinates_aoa"],
            Pnodes,
            controlpoints[i]["tangential"],
        )

        if sec > 4:

            lift_force[plates[sec][0], :] += Fnode[0]
            lift_force[plates[sec][1], :] += Fnode[1]
            lift_force[plates[sec][2], :] += Fnode[2]
            lift_force[plates[sec][3], :] += Fnode[3]
        else:
            lift_force[plates[sec][1], :] += Fnode[0]
            lift_force[plates[sec][0], :] += Fnode[1]
            lift_force[plates[sec][3], :] += Fnode[2]
            lift_force[plates[sec][2], :] += Fnode[3]

    return lift_force


# %% Above is URIs old method, that also considers the momentum balance.


# Distributing VSM-loads chordwise
def calculate_midpoints(wingpanels):
    midpoints = []
    for panel in wingpanels:
        p1, p2, p3, p4 = panel["p1"], panel["p2"], panel["p3"], panel["p4"]
        # Calculate midpoint of the trailing edge
        point_mid_te = p3 + 0.5 * (p4 - p3)
        # Calculate midpoint of the leading edge
        point_mid_le = p1 + 0.5 * (p2 - p1)
        midpoints.append({"point_mid_te": point_mid_te, "point_mid_le": point_mid_le})
    return midpoints


def interpolate_chordwise_points(midpoints, n_chordwise_aero_nodes):
    chordwise_points = np.empty(
        (len(midpoints), n_chordwise_aero_nodes, 3)
    )  # Initialize array to hold chordwise points
    for panel_idx, panel in enumerate(midpoints):
        point_mid_le = panel["point_mid_le"]
        point_mid_te = panel["point_mid_te"]

        # Linearly interpolate between point_mid_le and point_mid_te
        for i in range(n_chordwise_aero_nodes):
            # Calculate interpolation factor
            t = i / (n_chordwise_aero_nodes - 1)
            # Linear interpolation
            interpolated_point = (1 - t) * point_mid_le + t * point_mid_te
            chordwise_points[panel_idx, i] = interpolated_point
    return chordwise_points


def generate_distribution(n_chordwise_aero_nodes):
    x = np.linspace(0, 1, n_chordwise_aero_nodes)
    y = np.zeros_like(x)
    peak_index = int(len(x) * 0.25)
    y[:peak_index] = x[:peak_index] / (0.25)
    y[peak_index:] = (1 - x[peak_index:]) / (0.75)
    return x, y / np.sum(y)


# def aero2struc_NN_vsm(
#     n_chordwise_aero_nodes,
#     wing_aero,
#     force_aero_wing_VSM,
#     points_wing_segment_corners_aero_orderded,
#     points_structural_nodes,
#     plate_point_indices,
#     is_with_coupling_plot,
#     le_idx,
#     te_idx,
#     results_aero,
# ):

#     ## Step 1: Define a chordwise distribution
#     # 1.1) Calculate midpoints of the wingpanels
#     # midpoints = calculate_midpoints(wingpanels)
#     # THIS SHOULD NOT BE MIDPOINTS, BUT SHOULD BE THE POINTS ON THE AC TO CP LINE
#     # (AS THIS LINE IS NOT NECESSARILY IN THE MIDDLE OF THE PANEL)
#     midpoints = []
#     for panel in wing_aero.panels:
#         # define a vector from the aerodynamic center to the control point
#         vec_from_ac_to_cp_half_chord_length = (
#             panel.aerodynamic_center - panel.control_point
#         )

#         # define a leading-edge and trailing-edge point
#         point_mid_le = (
#             panel.aerodynamic_center + 0.5 * vec_from_ac_to_cp_half_chord_length
#         )
#         point_mid_te = panel.control_point - 0.5 * vec_from_ac_to_cp_half_chord_length
#         midpoints.append({"point_mid_te": point_mid_te, "point_mid_le": point_mid_le})

#     # 1.2) Interpolate chordwise points, creating a chordwise aero-mesh component
#     chordwise_points = interpolate_chordwise_points(midpoints, n_chordwise_aero_nodes)
#     # 1.3) Generate a distribution
#     chordwise_distribution, percentage_distribution = generate_distribution(
#         n_chordwise_aero_nodes
#     )

#     ## Step 2: Distribute the aerodynamic forces chordwise
#     force_aero_wing_VSM_distributed_chordwise = np.zeros(
#         (len(force_aero_wing_VSM), n_chordwise_aero_nodes, 3)
#     )
#     for i, points in enumerate(chordwise_points):
#         for j, point in enumerate(points):
#             force_aero_wing_VSM_distributed_chordwise[i, j] = (
#                 force_aero_wing_VSM[i] * percentage_distribution[j]
#             )

#     if is_with_coupling_plot:
#         plot_aerodynamic_forces_chordwise_distributed(
#             points_aero_chordwise=chordwise_points.reshape(-1, 3),
#             force_aero_chordwise=force_aero_wing_VSM_distributed_chordwise.reshape(
#                 -1, 3
#             ),
#             points_struc=points_wing_segment_corners_aero_orderded,
#         )

#     # force_aero_wing_segment_corners_struc_mesh = np.zeros_like(points_structural_nodes)
#     # # looping over each panel
#     # for i, (force_arr, chordwise_point_arr) in enumerate(
#     #     zip(force_aero_wing_VSM_distributed_chordwise, chordwise_points)
#     # ):
#     #     LAPs_of_this_panel = []
#     #     for index in plate_point_indices[i]:
#     #         LAPs_of_this_panel.append(points_structural_nodes[index])

#     #     f_LAPs_of_this_panel = map_idw(
#     #         chordwise_point_arr,
#     #         force_arr,
#     #         np.array(LAPs_of_this_panel).reshape(-1, 3),
#     #     )
#     #     # add these forces to the structural mesh
#     #     for idx, index in enumerate(plate_point_indices[i]):
#     #         force_aero_wing_segment_corners_struc_mesh[index] += f_LAPs_of_this_panel[
#     #             idx
#     #         ]

#     # return force_aero_wing_segment_corners_struc_mesh

#     from sklearn.neighbors import NearestNeighbors

#     force_struc = np.zeros_like(points_structural_nodes)

#     # Build two small KD-trees once per panel
#     for panel_i, (force_arr, chord_pts) in enumerate(
#         zip(force_aero_wing_VSM_distributed_chordwise, chordwise_points)
#     ):

#         # 1) find the two nearest LE nodes for *every* chordwise pt at once
#         le_coords = points_structural_nodes[le_idx]  # (n_le,3)
#         le_tree = NearestNeighbors(n_neighbors=2).fit(le_coords)
#         _, le_inds = le_tree.kneighbors(chord_pts)  # (n_chordwise, 2)
#         sel_le_glob = np.array(le_idx)[le_inds]  # (n_chordwise, 2)

#         # 2) same for TE
#         te_coords = points_structural_nodes[te_idx]
#         te_tree = NearestNeighbors(n_neighbors=2).fit(te_coords)
#         _, te_inds = te_tree.kneighbors(chord_pts)
#         sel_te_glob = np.array(te_idx)[te_inds]  # (n_chordwise, 2)

#         # 3) concatenate once
#         #    shape (n_chordwise, 4), each row = [le1,le2,te1,te2]
#         sel_idx = np.hstack((sel_le_glob, sel_te_glob))

#         # 4) for each chordwise pt, scatter to its 4 chosen nodes
#         for pt, frc, idx4 in zip(chord_pts, force_arr, sel_idx):
#             sel_coords = points_structural_nodes[idx4]  # (4,3)
#             f4 = map_idw(
#                 aero_nodes=pt[np.newaxis, :],  # (1,3)
#                 aero_forces=frc[np.newaxis, :],  # (1,3)
#                 struct_nodes=sel_coords,  # (4,3)
#             )  # returns (4,3)
#             for local_j, global_j in enumerate(idx4):
#                 force_struc[global_j] += f4[local_j]

#     if is_with_coupling_plot:

#         plot_aerodynamic_forces_chordwise_distributed(
#             points_aero_chordwise=chordwise_points.reshape(-1, 3),
#             force_aero_chordwise=force_aero_wing_VSM_distributed_chordwise.reshape(
#                 -1, 3
#             ),
#             points_struc=points_structural_nodes,
#             force_struc=force_struc,
#         )

#     # at the end:
#     return force_struc


# %% above is part of the old method, that used chordwise distributed nodes to handle the mapping.


# TODO: this should be placed in a more general plotting place/module
def plot_aerodynamic_forces_chordwise_distributed(
    points_aero_chordwise,
    force_aero_chordwise,
    points_struc,
    force_struc=None,
):
    import matplotlib.pyplot as plt

    # Create a new figure and set up 3D axes
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    # Scatter plot of chordwise points (blue)
    ax.scatter(
        points_aero_chordwise[:, 0],
        points_aero_chordwise[:, 1],
        points_aero_chordwise[:, 2],
        color="black",
        label="Chordwise Points",
    )

    # Quiver plot for the forces (red arrows)
    ax.quiver(
        points_aero_chordwise[:, 0],
        points_aero_chordwise[:, 1],
        points_aero_chordwise[:, 2],
        force_aero_chordwise[:, 0],
        force_aero_chordwise[:, 1],
        force_aero_chordwise[:, 2],
        # length=1,
        # normalize=True,
        length=0.01,
        color="black",
        label="Force Vectors",
    )

    if force_struc is None:
        # Scatter plot of structural nodes (wing segment corners) (green)
        ax.scatter(
            points_struc[:, 0],
            points_struc[:, 1],
            points_struc[:, 2],
            color="black",
            label="Wing Segment Corners",
        )

        # Annotate each point with its index
        for idx, point in enumerate(points_struc):
            ax.text(point[0], point[1], point[2], f"{idx}", color="black")
    else:
        # Scatter plot of structural nodes (wing segment corners) (green)
        ax.scatter(
            points_struc[:, 0],
            points_struc[:, 1],
            points_struc[:, 2],
            color="blue",
            label="Wing Segment Corners",
        )

        # Quiver plot for the forces on structural nodes (yellow arrows)
        ax.quiver(
            points_struc[:, 0],
            points_struc[:, 1],
            points_struc[:, 2],
            force_struc[:, 0],
            force_struc[:, 1],
            force_struc[:, 2],
            # length=1,
            length=0.01,
            # normalize=True,
            color="red",
            label="Force Vectors on Structural Nodes",
        )

    # Set equal scale for all axes
    points_all = np.concatenate((points_aero_chordwise, points_struc), axis=0)
    bb = points_all.max(axis=0) - points_all.min(axis=0)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_box_aspect(bb)
    ax.set_title("Aerodynamic Forces and Structural Nodes")
    ax.legend()
    plt.show()


def aero2struc_NN_vsm(
    force_aero_wing_VSM: np.ndarray,  # (n_panels,3)
    points_structural_nodes: np.ndarray,  # (n_struc,3)
    panel_cps: np.ndarray,  # (n_panels,3)
    panel_corner_map: np.ndarray,  # (n_panels,4)
    p: float = 2,
    eps: float = 1e-6,
    is_with_coupling_plot: bool = False,
):
    """
    Distribute each panel's resultant force (at its CoP) onto the four
    structural corner nodes given in panel_corner_map.
    """
    n_struc = len(points_structural_nodes)
    force_struc = np.zeros((n_struc, 3), dtype=float)

    for i, (cp, frc) in enumerate(zip(panel_cps, force_aero_wing_VSM)):
        sel_idx = panel_corner_map[i]  # [le_lo, le_hi, te_lo, te_hi]
        sel_coords = points_structural_nodes[sel_idx]  # (4,3)

        # true inverse-distance weighting across the 4 nodes
        d = np.linalg.norm(sel_coords - cp[None, :], axis=1)
        w = 1.0 / (d**p + eps)
        w /= np.sum(w)

        f_vals = w[:, None] * frc[None, :]  # (4,3)

        # accumulate
        for local_j, glob_j in enumerate(sel_idx):
            force_struc[glob_j] += f_vals[local_j]

    if is_with_coupling_plot:
        plot_aerodynamic_forces_chordwise_distributed(
            points_aero_chordwise=panel_cps,
            force_aero_chordwise=force_aero_wing_VSM,
            points_struc=points_structural_nodes,
            force_struc=force_struc,
        )

    return force_struc
