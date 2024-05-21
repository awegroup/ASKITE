# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 10:18:54 2023

@author: ocayon

Functions to plot the VSM-PSM in various styles
"""

# %% Importing functions

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.proj3d import proj_transform
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.proj3d import proj_transform
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import plotly.graph_objects as go
import plotly.offline as offline
import matplotlib

# own functions
from src.post_processing import post_processing_utils
from src.structural import structural_model


# %% Useful functions
def _arrow3D(ax, x, y, z, dx, dy, dz, *args, **kwargs):
    """Add an 3d arrow to an `Axes3D` instance."""

    arrow = Arrow3D(x, y, z, dx, dy, dz, *args, **kwargs)
    ax.add_artist(arrow)


class Arrow3D(FancyArrowPatch):
    def __init__(self, x, y, z, dx, dy, dz, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._xyz = (x, y, z)
        self._dxdydz = (dx, dy, dz)

    def draw(self, renderer):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        super().draw(renderer)

    def do_3d_projection(self, renderer=None):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))

        return np.min(zs)


# %% Plotting functions


def plot_VSM_PSM(
    wingpanels, controlpoints, rings, F, coord_L, pos, ax, col_kite, elev, plates
):
    # Aerodynamic part
    width = 1
    mksize = 1
    N_struct = 9
    N_split = int(len(wingpanels) / N_struct)
    secp = 10
    for panel in wingpanels:
        coord = np.array([panel["p1"], panel["p2"], panel["p3"]])
        ax.plot(
            coord[:, 0],
            coord[:, 1],
            coord[:, 2],
            "#000000",
            linewidth=width,
            linestyle="--",
        )
        coord = np.array([panel["p3"], panel["p4"]])
        ax.plot(coord[:, 0], coord[:, 1], coord[:, 2], "#000000", linewidth=width)
    for i in range(len(wingpanels)):
        sec = (N_struct - 1) - int((i + 1) / N_split - 0.01)
        if sec != secp:
            coord = np.array([wingpanels[i]["p1"], wingpanels[i]["p4"]])
            ax.plot(
                coord[:, 0], coord[:, 1], coord[:, 2], "#000000", linewidth=width * 5
            )
        coord = np.array([wingpanels[i]["p1"], wingpanels[i]["p2"]])
        ax.plot(coord[:, 0], coord[:, 1], coord[:, 2], "#000000", linewidth=width * 5)
        secp = sec
    coord = np.array([wingpanels[i]["p2"], wingpanels[i]["p3"]])
    ax.plot(coord[:, 0], coord[:, 1], coord[:, 2], "#000000", linewidth=width * 5)
    # for cp in controlpoints:
    # ax.plot(cp['coordinates'][0],cp['coordinates'][1],cp['coordinates'][2],'orange', marker = '.', markersize = mksize)
    # ax.plot(cp['coordinates_aoa'][0],cp['coordinates_aoa'][1],cp['coordinates_aoa'][2],'#00008B', marker = '.', markersize = mksize)
    for ring in rings:
        for filament in ring:
            if filament["id"] == "trailing_inf1" or filament["id"] == "trailing_inf2":
                coord = np.array([filament["x1"], filament["x1"] + filament["dir"] * 4])
                # ax.plot(coord[:,0], coord[:,1], coord[:,2], '#0D23C2',linewidth = width,linestyle = '--', alpha = 0.3)
            elif filament["id"] == "bound":
                coord = np.array([filament["x1"], filament["x2"]])
                ax.plot(
                    coord[:, 0],
                    coord[:, 1],
                    coord[:, 2],
                    "#1E90FF",
                    linewidth=width,
                    alpha=0.6,
                )
            else:
                coord = np.array([filament["x1"], filament["x2"]])
                # ax.plot(coord[:,0], coord[:,1], coord[:,2],'#0D23C2',linewidth = width,linestyle = '--', alpha=0.6)

    setattr(Axes3D, "arrow3D", _arrow3D)
    for i in range(len(F)):
        a = coord_L[i]
        b = (F[i][0]) / np.linalg.norm(F[int(len(F) / 2)][0]) * 3
        c = (F[i][1]) / np.linalg.norm(F[int(len(F) / 2)][1]) * 1
        ax.arrow3D(
            a[0],
            a[1],
            a[2],
            b[0],
            b[1],
            b[2],
            mutation_scale=5,
            linewidth=width,
            arrowstyle="-|>",
            fc="#2E8B57",
            ec="#2E8B57",
        )
        ax.arrow3D(
            a[0],
            a[1],
            a[2],
            c[0],
            c[1],
            c[2],
            mutation_scale=5,
            linewidth=width,
            arrowstyle="-|>",
            fc="#800000",
            ec="#800000",
        )

    # getting the connectivity
    # ci_bridle,cj_bridle = structural_mesher.get_bridle_connectivity()
    # plate_point_indices = structural_mesher.get_plate_point_indices()
    # ci_wing,cj_wing     = structural_mesher.get_wing_connectivity(plate_point_indices)
    # tube_indices        = structural_mesher.get_tube_indices()

    TE = [2, 8, 14, 20, 26, 32, 38, 44, 50]
    width = 0.7
    msize = 5

    # Fill canopy
    for i in range(len(plates)):
        x = np.array(
            [
                pos[plates[i][0], 0],
                pos[plates[i][1], 0],
                pos[plates[i][2], 0],
                pos[plates[i][3], 0],
            ]
        )
        y = np.array(
            [
                pos[plates[i][0], 1],
                pos[plates[i][1], 1],
                pos[plates[i][2], 1],
                pos[plates[i][3], 1],
            ]
        )
        z = np.array(
            [
                pos[plates[i][0], 2],
                pos[plates[i][1], 2],
                pos[plates[i][2], 2],
                pos[plates[i][3], 2],
            ]
        )

        # Create a polygonal surface between the lines
        vertices = [(xi, yi, zi) for xi, yi, zi in zip(x, y, z)]
        # vertices += [(xi, yi, zi) for xi, yi, zi in zip(x[::-1], z[::-1], np.zeros_like(x))]
        poly = Poly3DCollection([vertices], alpha=0.1, facecolor="black")
        ax.add_collection3d(poly)

    # ip = 1
    # for i in range(0, len(ci_kite)):
    #     if i in tube_idx:
    #         ax.plot([pos[ci_kite[i], 0], pos[cj_kite[i], 0]], [pos[ci_kite[i], 1], pos[cj_kite[i], 1]], [pos[ci_kite[i], 2], pos[cj_kite[i], 2]],color=col_kite,linewidth = width*6)
    #     elif i in TE:
    #         ip += 1
    #         ax.plot([pos[ci_kite[i], 0], pos[cj_kite[i], 0]], [pos[ci_kite[i], 1], pos[cj_kite[i], 1]], [pos[ci_kite[i], 2], pos[cj_kite[i], 2]],color=col_kite,linewidth = width)

    for i in range(0, len(ci)):
        if i > 21:
            ax.plot(
                [pos[ci[i], 0], pos[cj[i], 0]],
                [pos[ci[i], 1], pos[cj[i], 1]],
                [pos[ci[i], 2], pos[cj[i], 2]],
                color="#000000",
                linewidth=width,
                marker=".",
                markersize=msize,
            )
        else:
            ax.plot(
                [pos[ci[i], 0], pos[cj[i], 0]],
                [pos[ci[i], 1], pos[cj[i], 1]],
                [pos[ci[i], 2], pos[cj[i], 2]],
                color="#000000",
                linewidth=width,
                marker=".",
                markersize=msize,
            )

    # Plot KCU
    # ax.plot(pos[21, 0],pos[21, 1],pos[21, 2],color='#000000',linewidth = width,marker = 'v',markersize = msize*3)
    # Plot pulleys
    # ax.plot(pos[24, 0],pos[24, 1],pos[24, 2],color='orange',linewidth = width,marker = 'P',markersize = msize*1.5)
    # ax.plot(pos[28, 0],pos[28, 1],pos[28, 2],color='orange',linewidth = width,marker = 'P',markersize = msize*1.5)
    # Plot tapes
    # id1 = 21
    # id2 = 27
    # ax.plot([pos[id1, 0],pos[id2, 0]],[pos[id1, 1],pos[id2, 1]],[pos[id1, 2],pos[id2, 2]],color = '#333333',linewidth = width*3)
    # id2 = 22
    # ax.plot([pos[id1, 0],pos[id2, 0]],[pos[id1, 1],pos[id2, 1]],[pos[id1, 2],pos[id2, 2]],color = '#333333',linewidth = width*3)
    # id2 = 23
    # ax.plot([pos[id1, 0],pos[id2, 0]],[pos[id1, 1],pos[id2, 1]],[pos[id1, 2],pos[id2, 2]],color = '#333333',linewidth = width*3)
    # Plot KCU
    # ax.plot(pos[21, 0],pos[21, 1],pos[21, 2],color='#000000',linewidth = width,marker = 'v',markersize = msize*3)
    # Plot point mass particles
    for i in range(0, len(ci)):
        ax.plot(
            pos[ci[i], 0],
            pos[ci[i], 1],
            pos[ci[i], 2],
            color="#A605CD",
            linewidth=width,
            marker=".",
            markersize=msize,
        )
        ax.plot(
            pos[cj[i], 0],
            pos[cj[i], 1],
            pos[cj[i], 2],
            color="#A605CD",
            linewidth=width,
            marker=".",
            markersize=msize,
        )

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()
    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5 * max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

    legend_elements = [
        Line2D([0], [0], color=col_kite, lw=width * 5),
        Line2D([0], [0], color=col_kite, lw=width, ls="--"),
        Line2D([0], [0], color="#2E8B57", lw=width),
        Line2D([0], [0], color="#800000", lw=width),
        Line2D([0], [0], color="#1E90FF", lw=width),
        Line2D(
            [0], [0], marker=".", color="w", markerfacecolor="#A605CD", markersize=10
        ),
    ]

    ax.legend(
        legend_elements,
        [
            "Inflatable tubes",
            "VSM discretization",
            "Lift force",
            "Drag force",
            "Lifting line",
            "Point mass particles",
        ],
        frameon=False,
        loc="center right",
    )


# %% Plot in external window


def get_xyz_list_MinMax_range(N_D, points_list):  # use internally for making nice plots
    x_lst, y_lst, z_lst = [], [], []

    if N_D == 4:  # if it has 3 dimensions, i.e. list of list of lists of points
        flat_point_list = []
        for list_3 in points_list:
            for list in list_3:
                for point in list:
                    # Extract the x, y, and z coordinates from the array
                    x_lst.append(point[0])
                    y_lst.append(point[1])
                    z_lst.append(point[2])
                    flat_point_list.append(point[0])
                    flat_point_list.append(point[1])
                    flat_point_list.append(point[2])

        # Get the min and max range of the points, for aspect ratio
        min_range = min(flat_point_list)
        max_range = max(flat_point_list)

    elif N_D == 3:  # if it has 3 dimensions, i.e. list of lists of points
        flat_point_list = []
        for list in points_list:
            for point in list:
                # Extract the x, y, and z coordinates from the array
                x_lst.append(point[0])
                y_lst.append(point[1])
                z_lst.append(point[2])
                flat_point_list.append([point[0], point[1], point[2]])

        # Get the min and max range of the points, for aspect ratio
        min_range = min(np.concatenate(flat_point_list))
        max_range = max(np.concatenate(flat_point_list))

    elif N_D == 2:  # it has 2 dimensions, i.e. a list of point
        flat_point_list = []
        for point in points_list:
            # Extract the x, y, and z coordinates from the array
            x_lst.append(point[0])
            y_lst.append(point[1])
            z_lst.append(point[2])
            flat_point_list.append([point[0], point[1], point[2]])
            min_range = min(np.concatenate(flat_point_list))
            max_range = max(np.concatenate(flat_point_list))

    return x_lst, y_lst, z_lst, min_range, max_range


## defining a color bar generator work around
def gen_color_bar(line_trace):
    """
    Generates a trace which shows a colorbar based on a line plot.

    Relevant issue: https://github.com/plotly/plotly.py/issues/1085
    """
    return go.Scatter3d(
        x=line_trace.x,
        y=line_trace.y,
        z=line_trace.z,
        mode="markers",
        marker=go.scatter3d.Marker(
            color=line_trace.line.color,
            colorscale=line_trace.line.to_plotly_json()[
                "colorscale"
            ],  # Waiting on https://github.com/plotly/plotly.py/issues/1087
            showscale=line_trace.line.showscale,
            opacity=0.00000000000001,  # Make invisible, visible=False disables color bar
        ),
        hoverinfo="none",
        showlegend=False,
    )


def plot_kite(
    point_list,
    line_list=[[[], [], [], "black", 1, []]],
    title="3D Plot",
    surface_list=[[False, [], [], "lightgrey", 0.8]],
):
    """plots a 3D plot of the points offline (.html) using plotly
    input   : point_list,line_list,surface_list

        point_list = [[N_D,points_list,color,size],...]
            N_D               : number of dimensions (2,3,4)
            points_list       : list of points
            color             : color of the points
            size              : size of the points

        line_list = [[points,conn_i,conn_j,color,width,elongation_values],...]
            points            : list of points
            conn_i            : list of indices of the points to connect
            conn_j            : list of indices of the points to connect
            color             : color of the lines
            width             : width of the lines
            elongation_values : list of elongation values of the lines

        surface_list = [[boolean,points,plate_point_indices,colour,opacity],...]
            boolean            : True if the surface is to be plotted
            points             : list of points
            plate_point_indices: list of indices of the points to connect
            colour             : colour of the surface
            opacity            : opacity of the surface

    output  : 3D plot"""

    ### Points, gathering data and appending to data_list
    data_list = []
    min_range_list, max_range_list = [], []
    for i in range(len(point_list)):  # looping through each point_list
        N_D_i, points_list_i, color_i, size_i = point_list[i]
        x_lst_i, y_lst_i, z_lst_i, min_range_i, max_range_i = get_xyz_list_MinMax_range(
            N_D_i, points_list_i
        )

        # appending the min/max ranges
        min_range_list.append(min_range_i)
        max_range_list.append(max_range_i)

        data_list.append(
            go.Scatter3d(
                x=x_lst_i,
                y=y_lst_i,
                z=z_lst_i,
                mode="markers",
                marker=dict(color=color_i, size=size_i),
                showlegend=False,
            )
        )

    ### Lines, gathering data and appending to data_list

    ## Looping through each line_list
    for idx in range(len(line_list)):
        color_scale = "bluered"
        (
            points_idx,
            conn_i_idx,
            conn_j_idx,
            line_color_idx,
            line_width_idx,
            elongation_values,
        ) = line_list[idx]

        ## If no elongation values are given, do the usual print
        if len(elongation_values) == 0:
            for i in range(len(conn_i_idx)):
                data_list.append(
                    go.Scatter3d(
                        x=[points_idx[conn_i_idx[i]][0], points_idx[conn_j_idx[i]][0]],
                        y=[points_idx[conn_i_idx[i]][1], points_idx[conn_j_idx[i]][1]],
                        z=[points_idx[conn_i_idx[i]][2], points_idx[conn_j_idx[i]][2]],
                        mode="lines",
                        line=dict(width=line_width_idx, color=line_color_idx),
                        showlegend=False,
                    )
                )

        ## If elongation values are given, make the colorbar
        else:
            normalized_elongation_values = [
                (x - min(elongation_values))
                / (max(elongation_values) - min(elongation_values))
                for x in elongation_values
            ]

            def create_color_scale(values):
                colors = []
                for value in values:
                    # Interpolate between blue (0, 0, 1) and red (1, 0, 0)
                    r = 1 - value
                    b = value
                    g = 0.0
                    colors.append((r, g, b))

                return colors

            normalized_elongation_values_colors = create_color_scale(
                normalized_elongation_values
            )

            colors_based_on_percentage, labels, width_based_on_percentage = [], [], []
            for i, percentage in enumerate(elongation_values):
                if percentage > 5:
                    colors_based_on_percentage.append("red")
                    labels.append(f" 5% < dL  (very elongated)")
                    width_based_on_percentage.append(5)
                elif percentage > 2 and percentage < 5:
                    colors_based_on_percentage.append("orange")
                    labels.append(f" 2% < dL < 5%  (elongated)")
                    width_based_on_percentage.append(4)
                elif percentage > 1 and percentage < 2:
                    colors_based_on_percentage.append("yellow")
                    labels.append(f" 1% < dL < 2%  (elongated)")
                    width_based_on_percentage.append(3.5)
                elif percentage > -1 and percentage < 1:
                    colors_based_on_percentage.append("green")
                    labels.append(f"-1% < dL < 1% (acceptable)")
                    width_based_on_percentage.append(3)

                elif percentage > -2 and percentage < -1:
                    colors_based_on_percentage.append("turquoise")
                    labels.append(f"-2% < dL <-1%      (slack)")
                    width_based_on_percentage.append(3.5)
                elif percentage > -5 and percentage < -2:
                    colors_based_on_percentage.append("blue")
                    labels.append(f"-5% < dL <-2%      (slack)")
                    width_based_on_percentage.append(4)
                elif percentage < -5:
                    colors_based_on_percentage.append("black")
                    labels.append(f"      dL <-5% (very slack)")
                    width_based_on_percentage.append(5)

            # getting the points to plot all lines at once
            x_list, y_list, z_list = (
                np.zeros((len(conn_i_idx), 2)),
                np.zeros((len(conn_i_idx), 2)),
                np.zeros((len(conn_i_idx), 2)),
            )
            used_name_list = []
            for i in range(len(conn_i_idx)):
                x = [points_idx[conn_i_idx[i]][0], points_idx[conn_j_idx[i]][0]]
                y = [points_idx[conn_i_idx[i]][1], points_idx[conn_j_idx[i]][1]]
                z = [points_idx[conn_i_idx[i]][2], points_idx[conn_j_idx[i]][2]]
                x_list[i], y_list[i], z_list[i] = x, y, z

                if labels[i] not in used_name_list:
                    used_name_list.append(labels[i])
                    name_i = labels[i]
                    showlegend_i = True
                else:
                    name_i = None
                    showlegend_i = False

                data_list.append(
                    go.Scatter3d(
                        x=x,
                        y=y,
                        z=z,
                        mode="lines",
                        name=name_i,
                        line=go.scatter3d.Line(
                            color=colors_based_on_percentage[i],
                            width=width_based_on_percentage[i],
                            # showlegend = True,
                            # opacity = 1
                            # coloraxis=labels[i]
                            # colorbar = dict(
                            #     title = labels[i]
                            # ),
                            # colorscale=colors_based_on_percentage,
                            # showscale=True,
                            # symbol=labels[i]
                        ),
                        showlegend=showlegend_i,
                    )
                )

            # # getting the colorbar
            # line_trace_list = go.Scatter3d( x=x_list, y=y_list, z=z_list,
            #                                 mode='lines',
            #                                 line=go.scatter3d.Line(
            #                                     color=colors_based_on_percentage[i],
            #                                     colorscale=color_scale,
            #                                     showscale=True
            #                                     ),
            #                                 showlegend=False)

            # data_list.append(gen_color_bar(line_trace_list))

            # cmin_value = 1.2*np.min(elongation_values)
            # cmax_value = 1.2*np.max(elongation_values)

            # for i in range(len(conn_i_idx)):
            #     if i == 0: # only add the colorbar once
            #         fig.add_trace(go.Scatter3d(
            #             x=[points_idx[conn_i_idx[i]][0], points_idx[conn_j_idx[i]][0]],
            #             y=[points_idx[conn_i_idx[i]][1], points_idx[conn_j_idx[i]][1]],
            #             z=[points_idx[conn_i_idx[i]][2], points_idx[conn_j_idx[i]][2]],
            #             mode='lines',
            #             line=dict(
            #                 width=line_width_idx,
            #                 color=elongation_values[i],
            #                 colorscale=color_scale,
            #                 cmin = cmin_value,
            #                 cmax = cmax_value,
            #                 colorbar=dict(
            #                         thickness=15,
            #                         title = 'Elongation [mm]',
            #                         titlefont=dict(family='PT Sans Narrow', size=20, color='black'),
            #                         titleside = 'top',
            #                         # tickvals=[cmin_value, cmax_value],
            #                         # ticktext=[f'{cmin_value:.2f} (Slack)', f'{cmax_value:.2f} (Elongated)'],
            #                         outlinewidth=0)),
            #             showlegend= False))
            #     else:
            #         fig.add_trace(go.Scatter3d(
            #             x=[points_idx[conn_i_idx[i]][0], points_idx[conn_j_idx[i]][0]],
            #             y=[points_idx[conn_i_idx[i]][1], points_idx[conn_j_idx[i]][1]],
            #             z=[points_idx[conn_i_idx[i]][2], points_idx[conn_j_idx[i]][2]],
            #             mode='lines',
            #             line=dict(
            #                 width=line_width_idx,
            #                 color=elongation_values[i],
            #                 colorscale=color_scale),
            #             showlegend= False))

    ### Surfaces, gathering data and appending to data_list -as polygones

    ## looping through the surface_list = [surface_1,surface2]
    for surface in surface_list:
        # if surface[0], surface = [boolean,points,plate_point_indices,color,opacity]
        if surface[0] == True:
            # retrieving data
            (
                boolean,
                points,
                plate_point_indices,
                surface_color,
                surface_opacity,
            ) = surface

            for plate in plate_point_indices:  # loop through each plate ("segment")
                # define the x,y,z coordinates of each plate corner
                if len(plate) == 4:
                    x_surface = [
                        points[plate[0]][0],
                        points[plate[1]][0],
                        points[plate[2]][0],
                        points[plate[3]][0],
                    ]
                    y_surface = [
                        points[plate[0]][1],
                        points[plate[1]][1],
                        points[plate[2]][1],
                        points[plate[3]][1],
                    ]
                    z_surface = [
                        points[plate[0]][2],
                        points[plate[1]][2],
                        points[plate[2]][2],
                        points[plate[3]][2],
                    ]
                elif len(plate) == 3:
                    x_surface = [
                        points[plate[0]][0],
                        points[plate[1]][0],
                        points[plate[2]][0],
                    ]
                    y_surface = [
                        points[plate[0]][1],
                        points[plate[1]][1],
                        points[plate[2]][1],
                    ]
                    z_surface = [
                        points[plate[0]][2],
                        points[plate[1]][2],
                        points[plate[2]][2],
                    ]

                data_list.append(
                    go.Mesh3d(
                        x=x_surface,
                        y=y_surface,
                        z=z_surface,
                        color=surface_color,
                        opacity=surface_opacity,
                    )
                )

    ## defining the plot_range
    min_range = min(min_range_list) - 0.1 * max(max_range_list)
    max_range = max(max_range_list) + 0.1 * max(max_range_list)

    # # adding a custom legend
    # custom_legend = {'line':'A', 'circle-open': 'B', 'circle':'C'}

    # data_list.append(go.Scatter(
    #             x=[None],
    #             y=[None],
    #             mode="markers",
    #             name="A",
    #             marker=dict(size=7, color="blue", symbol='line'),
    #         ))

    ### Create a figure, set aspect ratio and plot data
    fig = go.FigureWidget(
        ## Layout part first
        layout=go.Layout(
            scene=dict(
                aspectmode="cube",
                bgcolor="#FFFFFF",
                xaxis=dict(
                    range=[min_range, max_range],
                    showgrid=False,
                    showticklabels=False,
                    zeroline=False,
                    showaxeslabels=False,
                    showbackground=False,
                    visible=False,
                ),
                yaxis=dict(
                    range=[min_range, max_range],
                    showgrid=False,
                    showticklabels=False,
                    zeroline=False,
                    showaxeslabels=False,
                    showbackground=False,
                    visible=False,
                ),
                zaxis=dict(
                    range=[min_range, max_range],
                    showgrid=False,
                    showticklabels=False,
                    zeroline=False,
                    showaxeslabels=False,
                    showbackground=False,
                    visible=False,
                ),
            ),
            title=dict(
                text=title,
                x=0.5,
                y=0.95,
                font=dict(family="PT Sans Narrow", size=25, color="black"),
            ),
        ),
        # Legend=dict(title = 'legend'),
        # Entering the data part
        data=data_list,
    )

    # Show the plot
    offline.plot(fig, filename="filename.html", auto_open=True)
    # fig.savefig("figure.svg")
    # fig.show()

    return ()


def plot_initial_geometry(config, points_between_dict):
    plot_points_initial = [2, config.kite.points_ini, "black", 2]
    plot_points_pulley = [
        2,
        config.kite.points_ini[config.kite.pulley.point_indices],
        "red",
        3,
    ]
    points_between_IN_dict = [
        config.kite.points_ini[int(key)] for key in points_between_dict.keys()
    ]
    points_edges_dict_values_0 = [
        config.kite.points_ini[int(value[0])] for value in points_between_dict.values()
    ]
    points_edges_dict_values_1 = [
        config.kite.points_ini[int(value[1])] for value in points_between_dict.values()
    ]

    # gathering all the plot_points
    plot_points_list = [plot_points_initial, plot_points_pulley]

    if points_edges_dict_values_0 != []:
        plot_points_list.append([2, points_edges_dict_values_0, "purple", 10])
    if points_edges_dict_values_1 != []:
        plot_points_list.append([2, points_edges_dict_values_1, "yellow", 10])
    if points_between_IN_dict != []:
        plot_points_list.append([2, points_between_IN_dict, "blue", 8])

    ##TODO: below is the code for testing-if you are accurately
    # capturing the rotational resistances particles
    ##TODO: remove this code at some point
    # ##TODO: should actually be able to split up the LE and Struts
    # # the indices are sorted from left-to-right
    # # the indices_mid are even sorted from TE to LE
    # # and all have the same length
    # indices_le = [int(value[0]) for value in points_between_dict.values()]
    # indices_mid = [int(key) for key in points_between_dict.keys()]
    # indices_te = [int(value[1]) for value in points_between_dict.values()]

    # le_rotational_resistance_dict = {}
    # # making unique copies
    # indices_le_unique = np.unique(np.copy(indices_le))
    # indices_le_unique = sorted(
    #     indices_le_unique, key=lambda idx: -config.kite.points_ini[idx][1]
    # )
    # indices_mid_unique = np.unique(np.copy(indices_mid))
    # indices_te_unique = np.unique(np.copy(indices_te))
    # indices_te_unique = sorted(
    #     indices_te_unique, key=lambda idx: -config.kite.points_ini[idx][1]
    # )

    # # first point
    # le_rotational_resistance_dict[str(indices_le_unique[0])] = [
    #     indices_le_unique[1],
    #     indices_te_unique[0],
    # ]
    # # middle points
    # for idx in range(1, len(indices_le_unique) - 1):
    #     le_rotational_resistance_dict[str(indices_le_unique[idx])] = [
    #         indices_le_unique[idx - 1],
    #         indices_le_unique[idx + 1],
    #     ]
    # # last point
    # le_rotational_resistance_dict[str(indices_le_unique[-1])] = [
    #     indices_le_unique[-2],
    #     indices_te_unique[-1],
    # ]

    # # TODO: continue with the middle-points
    # indices_le_sorted = sorted(
    #     indices_le, key=lambda idx: -config.kite.points_ini[idx][1]
    # )
    # indices_te_sorted = sorted(
    #     indices_te, key=lambda idx: -config.kite.points_ini[idx][1]
    # )
    # strut_rotational_resistance_dict = {}

    # ## loop through all the le points, knowing that when there is no more te point, the strut is finished
    # # the indices_mid are aranged from left-to-right and from TE-to-LE
    # for i, index in enumerate(indices_le_sorted[:-1]):
    #     # TODO: implement some logic, something like below
    #     # if has another point back and front
    #     if (
    #         indices_le_sorted[i] == indices_le_sorted[i + 1]
    #         and indices_le_sorted[i] == indices_le_sorted[i - 1]
    #     ):
    #         strut_rotational_resistance_dict[str(indices_mid[i])] = [
    #             indices_mid[i - 1],
    #             indices_mid[i + 1],
    #         ]
    #     elif indices_le_sorted[i] == indices_le_sorted[i + 1]:
    #         strut_rotational_resistance_dict[str(indices_mid[i])] = [
    #             indices_te_sorted[i],
    #             indices_mid[i + 1],
    #         ]
    #     elif indices_le_sorted[i] == indices_le_sorted[i - 1]:
    #         strut_rotational_resistance_dict[str(indices_mid[i])] = [
    #             indices_mid[i - 1],
    #             indices_le_sorted[i],
    #         ]
    #     elif (
    #         indices_le_sorted[i] != indices_le_sorted[i + 1]
    #         and indices_le_sorted[i] != indices_le_sorted[i - 1]
    #     ):
    #         strut_rotational_resistance_dict[str(indices_mid[i])] = [
    #             indices_te_sorted[i],
    #             indices_le_sorted[i],
    #         ]

    # ## Printing the rotational resistance pairing dict.
    # # test_dict = le_rotational_resistance_dict
    # test_dict = strut_rotational_resistance_dict
    # key = str(indices_mid[38])

    # # test_dict = le_rotational_resistance_dict
    # # key = str(indices_le_unique[5])
    # # middle point
    # plot_points_list.append(
    #     [
    #         2,
    #         [
    #             config.kite.points_ini[int(key)],
    #             config.kite.points_ini[int(key)],
    #         ],
    #         "orange",
    #         15,
    #     ]
    # )
    # # connecter points
    # plot_points_list.append(
    #     [
    #         2,
    #         config.kite.points_ini[test_dict[key]],
    #         "deeppink",
    #         15,
    #     ]
    # )

    # test_dict = le_rotational_resistance_dict
    # key = str(indices_le_unique[13])
    # # middle point
    # plot_points_list.append(
    #     [
    #         2,
    #         [
    #             config.kite.points_ini[int(key)],
    #             config.kite.points_ini[int(key)],
    #         ],
    #         "orange",
    #         15,
    #     ]
    # )
    # # connecter points
    # plot_points_list.append(
    #     [
    #         2,
    #         config.kite.points_ini[test_dict[key]],
    #         "deeppink",
    #         15,
    #     ]
    # )

    plot_lines_pulley = [
        config.kite.points_ini,
        config.kite.pulley.ci,
        config.kite.pulley.cj,
        "red",
        4,
        [],
    ]
    plot_lines_bridle = [
        config.kite.points_ini,
        config.kite.connectivity.bridle_ci,
        config.kite.connectivity.bridle_cj,
        "black",
        1.5,
        [],
    ]
    plot_lines_wing = [
        config.kite.points_ini,
        config.kite.connectivity.wing_ci,
        config.kite.connectivity.wing_cj,
        "seagreen",
        2,
        [],
    ]
    plot_lines_tube = [
        config.kite.points_ini,
        config.kite.connectivity.wing_ci[config.kite.connectivity.tube_line_indices],
        config.kite.connectivity.wing_cj[config.kite.connectivity.tube_line_indices],
        "seagreen",
        10,
        [],
    ]

    plot_kite(
        plot_points_list,
        [plot_lines_pulley, plot_lines_bridle, plot_lines_wing, plot_lines_tube],
        title=config.kite_name,
    )
    return


# %% Aerodynamic Plots


def plot_geometry_aero_old(
    wingpanels, controlpoints, rings, F, coord_L, ax, plot, N_segments_struc
):
    width = 1
    mksize = 5
    # N_segments_struc = 13
    N_split = int(len(wingpanels) / N_segments_struc)
    secp = N_segments_struc + 1
    if plot == "True":
        for panel in wingpanels:
            coord = np.array([panel["p1"], panel["p2"], panel["p3"], panel["p4"]])
            ax.plot(coord[:, 0], coord[:, 1], coord[:, 2], "#000000", linewidth=width)
        for i in range(len(wingpanels)):
            sec = (N_segments_struc - 1) - int((i + 1) / N_split - 0.01)
            if sec != secp:
                coord = np.array([wingpanels[i]["p1"], wingpanels[i]["p4"]])
                ax.plot(
                    coord[:, 0],
                    coord[:, 1],
                    coord[:, 2],
                    "#000000",
                    linewidth=width * 5,
                )
            coord = np.array([wingpanels[i]["p1"], wingpanels[i]["p2"]])
            ax.plot(
                coord[:, 0], coord[:, 1], coord[:, 2], "#000000", linewidth=width * 5
            )
            secp = sec
        coord = np.array([wingpanels[i]["p2"], wingpanels[i]["p3"]])
        ax.plot(coord[:, 0], coord[:, 1], coord[:, 2], "#000000", linewidth=width * 5)
        for cp in controlpoints:
            ax.plot(
                cp["coordinates"][0],
                cp["coordinates"][1],
                cp["coordinates"][2],
                "orange",
                marker=".",
                markersize=mksize,
            )
            ax.plot(
                cp["coordinates_aoa"][0],
                cp["coordinates_aoa"][1],
                cp["coordinates_aoa"][2],
                "#23B52180",
                marker=".",
                markersize=mksize,
            )
        for ring in rings:
            for filament in ring:
                if (
                    filament["id"] == "trailing_inf1"
                    or filament["id"] == "trailing_inf2"
                ):
                    coord = np.array(
                        [filament["x1"], filament["x1"] + filament["dir"] * 4]
                    )
                    ax.plot(
                        coord[:, 0],
                        coord[:, 1],
                        coord[:, 2],
                        "#0D23C233",
                        linewidth=width,
                        linestyle="--",
                    )
                else:
                    coord = np.array([filament["x1"], filament["x2"]])
                    ax.plot(
                        coord[:, 0],
                        coord[:, 1],
                        coord[:, 2],
                        "#0D23C233",
                        linewidth=width,
                        linestyle="--",
                    )

        setattr(Axes3D, "arrow3D", _arrow3D)
        for i in range(len(F)):
            a = coord_L[i]
            b = (
                (F[i][0] + F[i][1])
                / np.linalg.norm(F[int(len(F) / 2)][0] + F[int(len(F) / 2)][1])
                * 2
            )
            ax.arrow3D(
                a[0],
                a[1],
                a[2],
                b[0],
                b[1],
                b[2],
                mutation_scale=5,
                linewidth=width,
                arrowstyle="-|>",
                fc="#23B521",
                ec="#23B521",
            )

        legend_elements = [
            Line2D([0], [0], color="#000000", lw=width * 5),
            Line2D([0], [0], color="#000000", lw=width),
            Line2D([0], [0], color="#0D23C280", linestyle="--", lw=width),
            # Line2D([0], [0], color='#23B521', lw=width),
            Line2D(
                [0],
                [0],
                marker=".",
                color="w",
                markerfacecolor="#23B521",
                markersize=10,
            ),
            Line2D(
                [0], [0], marker=".", color="w", markerfacecolor="orange", markersize=10
            ),
        ]

        ax.legend(
            legend_elements,
            [
                "Inflatable tubes",
                "Aerodynamic discretization",
                "Horseshoe vortices",
                "Lifting Line points ($1/4$c)",
                "Control points ($3/4$c)",
            ],
            frameon=False,
            loc="upper right",
            bbox_to_anchor=(1.25, 1),
        )


def plot_particle_system_model(points, elongation_values, fig, ax, config):
    ## initializing
    width_bridles = 1.5
    width_wings = 1.0

    # colourmap = matplotlib.colors.ListedColormap(
    #     ["blue", "black", "black", "darkorange"]
    # )
    # colourmap.set_over("red")
    # colourmap.set_under("navy")
    # bounds = [-1, -0.5, 0.5, 1]

    ## dark mode
    colourmap = matplotlib.cm.get_cmap("viridis")
    bounds = [-2.5, -2.0, -1.5, -0.5, 0.5, 1.5, 2.0, 2.5]
    bounds = [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]
    norm = matplotlib.colors.BoundaryNorm(bounds, colourmap.N)

    for idx, (i, j) in enumerate(
        zip(config.kite.connectivity.bridle_ci, config.kite.connectivity.bridle_cj)
    ):
        if idx in config.kite.pulley.line_indices:
            lwidth = 1.8 * width_bridles
            lstyle = "solid"  # (0, (5, 1))  # densely dashed
        else:
            lwidth = 1.0 * width_bridles
            lstyle = "solid"

        p1 = points[i]
        p2 = points[j]
        ax.plot(
            [p1[0], p2[0]],
            [p1[1], p2[1]],
            [p1[2], p2[2]],
            linestyle=lstyle,
            c=colourmap(norm(elongation_values[idx])),
            lw=lwidth,
        )

    for idx, (i, j) in enumerate(
        zip(config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj)
    ):
        if idx in config.kite.connectivity.tube_line_indices:
            lwidth = 2.5
        else:
            lwidth = 1.0
        p1 = points[i]
        p2 = points[j]
        ax.plot(
            [p1[0], p2[0]],
            [p1[1], p2[1]],
            [p1[2], p2[2]],
            linestyle="-",
            c=colourmap(norm(elongation_values[idx])),
            lw=lwidth,
        )

    ## configuring colour bar
    ax_pos = ax.get_position()

    cax = fig.add_axes(
        [
            ax_pos.x0 + ax_pos.width * 0.3,
            ax_pos.y0,  # ax_pos.y0 - 0.05,
            ax_pos.width * 0.5,
            ax_pos.height * 0.02,
        ]
    )

    fig.colorbar(
        matplotlib.cm.ScalarMappable(norm=norm, cmap=colourmap),
        cax=cax,
        boundaries=[-10] + bounds + [10],
        extend="both",
        extendfrac="auto",
        ticks=bounds,
        spacing="uniform",
        orientation="horizontal",
        label="Element elongation %",
    )

    # plotting the knots and pulleys
    for i, point in enumerate(points):
        if i in config.kite.pulley.point_indices:
            colour = "black"
            msize = 4.0
            mtype = "o"  # Diamond
        else:
            colour = "black"
            msize = 3.0
            mtype = "o"
        ax.plot(point[0], point[1], point[2], colour, marker=mtype, markersize=msize)


def plot_geometry_aero(
    wingpanels,
    controlpoints,
    rings,
    F_rel,
    coord_L,
    is_with_aero_plot,
    ax,
    N_segments_struc,
):
    # colors
    # col_3_4_cp_points = "orange"  #'orange
    # col_1_4_ll_points = "#23B52180"  #'#23B52180'
    # col_horsehoe = "#23B52180"  #'0D23C233'
    # col_vector = "limegreen"  #'#23B521'
    # col_aero_mesh = "orchid"
    # col_surfaces = "orchid"
    col_3_4_cp_points = "orange"  #'orange
    col_1_4_ll_points = "#23B52180"  #'#23B52180'
    colorrrr = "deeppink"
    col_horsehoe = colorrrr  #'0D23C233'
    col_vector = colorrrr  #'#23B521'
    col_aero_mesh = colorrrr
    col_surfaces = colorrrr

    width_multiplier_aero_mesh = 1.4
    width = 1.2
    mksize = 7

    opacity_surfaces = 0.15
    # N_segments_struc = 13

    ## Plotting forces
    setattr(Axes3D, "arrow3D", _arrow3D)
    for i in range(len(F_rel)):
        a = coord_L[i]
        F = F_rel
        b = (
            (F[i][0] + F[i][1])
            / np.linalg.norm(F[int(len(F) / 2)][0] + F[int(len(F) / 2)][1])
            * 2
        )
        # b = (force_aero_wing_refined[i])/np.linalg.norm(force_aero_wing_refined[int(len(force_aero_wing_refined)/2)])*2
        ax.arrow3D(
            a[0],
            a[1],
            a[2],
            b[0],
            b[1],
            b[2],
            mutation_scale=5,
            linewidth=width * 2.0,
            arrowstyle="-|>",
            fc=col_vector,
            ec=col_vector,
        )

    N_split = int(len(wingpanels) / N_segments_struc)
    secp = N_segments_struc + 1
    if is_with_aero_plot:
        for panel in wingpanels:
            vertices = [panel["p1"], panel["p2"], panel["p3"], panel["p4"]]
            ## plotting the vertices
            for vertex in vertices:
                ax.plot(
                    vertex[0],
                    vertex[1],
                    vertex[2],
                    col_aero_mesh,
                    marker=".",
                    markersize=mksize,
                )
            ## colouring the surfaces
            ax.add_collection3d(
                Poly3DCollection(
                    [vertices],
                    facecolors=col_surfaces,
                    alpha=opacity_surfaces,
                    edgecolors="none",
                )
            )

            coord = np.array([panel["p1"], panel["p2"], panel["p3"], panel["p4"]])
            ax.plot(
                coord[:, 0],
                coord[:, 1],
                coord[:, 2],
                col_aero_mesh,
                linewidth=width * width_multiplier_aero_mesh,
            )
        for i in range(len(wingpanels)):
            sec = (N_segments_struc - 1) - int((i + 1) / N_split - 0.01)
            if sec != secp:
                coord = np.array([wingpanels[i]["p1"], wingpanels[i]["p4"]])
                ax.plot(
                    coord[:, 0],
                    coord[:, 1],
                    coord[:, 2],
                    col_aero_mesh,
                    linewidth=width * width_multiplier_aero_mesh,
                )
            coord = np.array([wingpanels[i]["p1"], wingpanels[i]["p2"]])
            ax.plot(
                coord[:, 0],
                coord[:, 1],
                coord[:, 2],
                col_aero_mesh,
                linewidth=width * width_multiplier_aero_mesh,
            )
            secp = sec
        coord = np.array([wingpanels[i]["p2"], wingpanels[i]["p3"]])
        # ax.plot(
        #     coord[:, 0], coord[:, 1], coord[:, 2], col_aero_mesh, linewidth=width * 5
        # )
        for cp in controlpoints:
            ax.plot(
                cp["coordinates"][0],
                cp["coordinates"][1],
                cp["coordinates"][2],
                col_3_4_cp_points,
                marker=".",
                markersize=mksize,
            )
            ax.plot(
                cp["coordinates_aoa"][0],
                cp["coordinates_aoa"][1],
                cp["coordinates_aoa"][2],
                col_1_4_ll_points,
                marker=".",
                markersize=mksize,
            )
        for ring in rings:
            for filament in ring:
                if (
                    filament["id"] == "trailing_inf1"
                    or filament["id"] == "trailing_inf2"
                ):
                    coord = np.array(
                        [filament["x1"], filament["x1"] + filament["dir"] * 4]
                    )
                    ax.plot(
                        coord[:, 0],
                        coord[:, 1],
                        coord[:, 2],
                        col_horsehoe,
                        linewidth=width,
                        linestyle="--",
                    )
                else:
                    coord = np.array([filament["x1"], filament["x2"]])
                    ax.plot(
                        coord[:, 0],
                        coord[:, 1],
                        coord[:, 2],
                        col_horsehoe,
                        linewidth=width,
                        linestyle="--",
                    )

        legend_elements = [
            # Line2D([0], [0], color=col_aero_mesh, lw=width * 3),
            Line2D([0], [0], color=col_aero_mesh, lw=width),
            Line2D([0], [0], color=col_horsehoe, linestyle="--", lw=width),
            # Line2D([0], [0], color='#23B521', lw=width),
            Line2D(
                [0],
                [0],
                marker=".",
                color="w",
                markerfacecolor=col_1_4_ll_points,
                markersize=10,
            ),
            Line2D(
                [0],
                [0],
                marker=".",
                color="w",
                markerfacecolor=col_3_4_cp_points,
                markersize=10,
            ),
        ]
        # TODO: can also add other legend elements than only aero for
        ## configuring the legend
        legend_names = [
            "Aerodynamic discretization",
            "Horseshoe vortices",
            "Lifting Line points ($1/4$c)",
            "Control points ($3/4$c)",
        ]

        # ax.legend(
        #     legend_elements,
        #     legend_names,
        #     frameon=False,
        #     # loc="upper left",
        #     loc="upper center"
        #     # bbox_to_anchor=(1.25, 1),
        # )
        ## Plotting Legend as well

    return


# def plot_panel(wingpanels, controlpoints, rings, F, coord_L, ax, plot):
#     width = 1
#     mksize = 5
#     N_struct = 9
#     N_split = int(len(wingpanels) / N_struct)
#     secp = 10
#     if plot == "True":
#         for panel in wingpanels:
#             coord = np.array([panel["p1"], panel["p2"], panel["p3"], panel["p4"]])
#             ax.plot(coord[:, 0], coord[:, 1], coord[:, 2], "#000000", linewidth=width)
#         for i in range(len(wingpanels)):
#             sec = (N_struct - 1) - int((i + 1) / N_split - 0.01)
#             if sec != secp:
#                 coord = np.array([wingpanels[i]["p1"], wingpanels[i]["p4"]])
#                 ax.plot(
#                     coord[:, 0],
#                     coord[:, 1],
#                     coord[:, 2],
#                     "#000000",
#                     linewidth=width * 5,
#                 )
#             coord = np.array([wingpanels[i]["p1"], wingpanels[i]["p2"]])
#             ax.plot(
#                 coord[:, 0], coord[:, 1], coord[:, 2], "#000000", linewidth=width * 5
#             )
#             secp = sec
#         coord = np.array([wingpanels[i]["p2"], wingpanels[i]["p3"]])
#         ax.plot(coord[:, 0], coord[:, 1], coord[:, 2], "#000000", linewidth=width * 5)
#         for cp in controlpoints:
#             ax.plot(
#                 cp["coordinates"][0],
#                 cp["coordinates"][1],
#                 cp["coordinates"][2],
#                 "orange",
#                 marker=".",
#                 markersize=mksize,
#             )
#             ax.plot(
#                 cp["coordinates_aoa"][0],
#                 cp["coordinates_aoa"][1],
#                 cp["coordinates_aoa"][2],
#                 "#23B52180",
#                 marker=".",
#                 markersize=mksize,
#             )
#         for ring in rings:
#             for filament in ring:
#                 if (
#                     filament["id"] == "trailing_inf1"
#                     or filament["id"] == "trailing_inf2"
#                 ):
#                     coord = np.array(
#                         [filament["x1"], filament["x1"] + filament["dir"] * 4]
#                     )
#                     ax.plot(
#                         coord[:, 0],
#                         coord[:, 1],
#                         coord[:, 2],
#                         "#0D23C233",
#                         linewidth=width,
#                         linestyle="--",
#                     )
#                 else:
#                     coord = np.array([filament["x1"], filament["x2"]])
#                     ax.plot(
#                         coord[:, 0],
#                         coord[:, 1],
#                         coord[:, 2],
#                         "#0D23C233",
#                         linewidth=width,
#                         linestyle="--",
#                     )

#         setattr(Axes3D, "arrow3D", _arrow3D)
#         for i in range(len(F)):
#             a = coord_L[i]
#             b = (
#                 (F[i][0] + F[i][1])
#                 / np.linalg.norm(F[int(len(F) / 2)][0] + F[int(len(F) / 2)][1])
#                 * 2
#             )
#             ax.arrow3D(
#                 a[0],
#                 a[1],
#                 a[2],
#                 b[0],
#                 b[1],
#                 b[2],
#                 mutation_scale=5,
#                 linewidth=width,
#                 arrowstyle="-|>",
#                 fc="#23B521",
#                 ec="#23B521",
#             )

#         legend_elements = [
#             Line2D([0], [0], color="#000000", lw=width * 5),
#             Line2D([0], [0], color="#000000", lw=width),
#             Line2D([0], [0], color="#0D23C280", linestyle="--", lw=width),
#             Line2D([0], [0], color="#23B521", lw=width),
#             Line2D(
#                 [0],
#                 [0],
#                 marker=".",
#                 color="w",
#                 markerfacecolor="#23B521",
#                 markersize=10,
#             ),
#             Line2D(
#                 [0], [0], marker=".", color="w", markerfacecolor="orange", markersize=10
#             ),
#         ]

#         ax.legend(
#             legend_elements,
#             [
#                 "Inflatable tubes",
#                 "Aerodynamic discretization",
#                 "Horseshoe vortices",
#                 "Local Aerodynamic Forces",
#                 "Lifting Line points ($1/4$c)",
#                 "Control points ($3/4$c)",
#             ],
#             frameon=False,
#             loc="center",
#         )


def _arrow3D(ax, x, y, z, dx, dy, dz, *args, **kwargs):
    """Add an 3d arrow to an `Axes3D` instance."""

    arrow = Arrow3D(x, y, z, dx, dy, dz, *args, **kwargs)
    ax.add_artist(arrow)


class Arrow3D(FancyArrowPatch):
    def __init__(self, x, y, z, dx, dy, dz, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._xyz = (x, y, z)
        self._dxdydz = (dx, dy, dz)

    def draw(self, renderer):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        super().draw(renderer)

    def do_3d_projection(self, renderer=None):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))

        return np.min(zs)


def plot_wind_vector(
    VEL_WIND, VEL_KITE, vel_app, controlpoints, ax, v_wind_color="blue"
):
    # Calculate the maximum dimension of the plot
    max_dimension = 0
    for cp in controlpoints:
        max_dimension = max(max_dimension, np.linalg.norm(cp["coordinates"]))

    # Calculate scaling factor for wind vectors
    max_vel = max(
        np.linalg.norm(VEL_WIND), np.linalg.norm(VEL_KITE), np.linalg.norm(vel_app)
    )
    vel_scaling = 0.3 * max_dimension / max_vel
    width = 2

    ### Plotting flow vectors
    if np.linalg.norm(VEL_KITE) == 0:
        label = "Va (=Vw)"
        setattr(Axes3D, "arrow3D", _arrow3D)
        scaled_vector = -1 * VEL_WIND * vel_scaling  # Scale the vector
        arrow_start = (0, 0, 0)
        arrow_end = scaled_vector
        ax.arrow3D(
            arrow_start[0],
            arrow_start[1],
            arrow_start[2],
            arrow_end[0],
            arrow_end[1],
            arrow_end[2],
            mutation_scale=10,
            linewidth=width,
            arrowstyle="<|-",
            fc=v_wind_color,
            ec=v_wind_color,
        )
        arrow_end = arrow_end - 0.5
        ax.text(
            arrow_end[0],
            arrow_end[1],
            arrow_end[2],
            label,
            color=v_wind_color,
            fontsize=17,
        )

    else:  # if Vk nonzero, we plot all three vector
        ## plot inflow vectors
        vectors = [VEL_WIND, -VEL_KITE, vel_app]
        labels = ["Vw", "-Vk", "Va"]

        setattr(Axes3D, "arrow3D", _arrow3D)
        for vector, label in zip(vectors, labels):
            scaled_vector = -1 * vector * vel_scaling  # Scale the vector
            arrow_start = (0, 0, 0)
            arrow_end = scaled_vector
            ax.arrow3D(
                arrow_start[0],
                arrow_start[1],
                arrow_start[2],
                arrow_end[0],
                arrow_end[1],
                arrow_end[2],
                mutation_scale=10,
                linewidth=width,
                arrowstyle="<|-",
                fc=v_wind_color,
                ec=v_wind_color,
            )
            arrow_end = arrow_end - 0.5
            ax.text(
                arrow_end[0],
                arrow_end[1],
                arrow_end[2],
                label,
                color=v_wind_color,
                fontsize=17,
            )


def plot_axis(ax, axis_color="black"):
    ### Simulation Axis
    # Create arrow patches
    vectors = [np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])]
    labels = ["X", "Y", "Z"]
    width = 2.5

    setattr(Axes3D, "arrow3D", _arrow3D)
    for vector, label in zip(vectors, labels):
        arrow_start = (0, 0, 0)
        arrow_end = 1.2 * vector
        ax.arrow3D(
            arrow_start[0],
            arrow_start[1],
            arrow_start[2],
            arrow_end[0],
            arrow_end[1],
            arrow_end[2],
            mutation_scale=10,
            linewidth=width,
            arrowstyle="-|>",
            fc=axis_color,
            ec=axis_color,
        )
        arrow_end = arrow_end * 1.05
        ax.text(
            arrow_end[0],
            arrow_end[1],
            arrow_end[2],
            label,
            color=axis_color,
            fontsize=18,
        )

    # ### Wind Axis
    # labels = ["Xw", "Yw", "Zw"]
    # # Create arrow patches
    # for i, value in enumerate(VEL_WIND):
    #     if VEL_WIND[i] != 0:
    #         if i == 0:
    #             arrow_ends = [
    #                 np.array([-1, 0, 0]),
    #                 np.array([-2, 1, 0]),
    #                 np.array([-2, 0, 1]),
    #             ]
    #         elif i == 2:
    #             arrow_ends = [
    #                 np.array([-1, 0, 0]),
    #                 np.array([-2, 1, 0]),
    #                 np.array([-2, 0, 1]),
    #             ]

    # arrow_color = "blue"
    # width = 2

    # setattr(Axes3D, "arrow3D", _arrow3D)
    # for arrow_end, label in zip(arrow_ends, labels):
    #     arrow_start = (-2, 0, 0)
    #     print("arrow_end", arrow_end)
    #     ax.arrow3D(
    #         arrow_start[0],
    #         arrow_start[1],
    #         arrow_start[2],
    #         arrow_end[0],
    #         arrow_end[1],
    #         arrow_end[2],
    #         mutation_scale=10,
    #         linewidth=width,
    #         arrowstyle="-|>",
    #         fc=arrow_color,
    #         ec=arrow_color,
    #     )
    #     arrow_end = arrow_end * 1.05
    #     ax.text(
    #         arrow_end[0],
    #         arrow_end[1],
    #         arrow_end[2],
    #         label,
    #         color=arrow_color,
    #         fontsize=15,
    #     )

    # Adding x, y, z axes vectors
    # axis_factor = 0.1
    # axis_scaling = max_dimension*axis_factor

    # x_displacement = -0 * axis_scaling
    # vector_size = 1.5 * axis_scaling
    # font_size = 27
    # ax.quiver(x_displacement, 0, 0, vector_size, 0, 0, color='black', arrow_length_ratio=0.2, linewidth=5)
    # ax.text(x_displacement+vector_size*1.2, 0, 0, 'x', color='black', fontsize=font_size)
    # ax.quiver(x_displacement, 0, 0, 0, vector_size, 0, color='black', arrow_length_ratio=0.2, linewidth=5)
    # ax.text(x_displacement, vector_size*1.2, 0, 'y', color='black', fontsize=font_size)
    # ax.quiver(x_displacement, 0, 0, 0, 0, vector_size, color='black', arrow_length_ratio=0.2, linewidth=5)
    # ax.text(x_displacement, 0, vector_size*1.2, 'z', color='black', fontsize=font_size)


def plot_aero(
    points,
    elongation_values,
    vel_app,
    wingpanels,
    controlpoints,
    rings,
    coord_L,
    F_rel,
    config,
    elev=10,
    azim=230,
    it_number=00,
):
    # for creating a responsive plot
    # matplotlib widget
    VEL_WIND = config.vel_wind
    VEL_KITE = config.vel_kite
    N_segments = config.kite.n_segments

    # style = "dark_background"
    style = "default"
    # style = "seaborn-whitegrid"
    # style = "seaborn-darkgrid"

    plt.style.use(style)
    # height = 0.7 * max(points[:, 2])
    width = 1.0 * max(config.kite.points_ini[:, 1])
    width = np.round(width, 0)
    fig = plt.figure(figsize=(width, width))
    # fig = plt.figure()
    # ax = fig.add_subplot(1, 1, 1, projection="3d")
    ax = plt.axes(projection="3d")

    if style == "dark_background":
        axis_color = "white"
    else:
        axis_color = "black"

    plot_axis(ax, axis_color)
    plot_wind_vector(VEL_WIND, VEL_KITE, vel_app, controlpoints, ax, "dodgerblue")
    # TODO: legend is off
    plot_geometry_aero(
        wingpanels,
        controlpoints,
        rings,
        F_rel,
        coord_L,
        config.is_with_aero_geometry,
        ax,
        N_segments,
    )

    plot_particle_system_model(points, elongation_values, fig, ax, config)

    ax.axis("off")
    ax.set_xlim((-width, width))
    ax.set_ylim((-width, width))
    ax.set_zlim((0, 2 * width))
    ax.view_init(elev=elev, azim=azim)
    # bbox = fig.bbox_inches.from_bounds(1.5, 3, 7, 4.5)

    # TODO: add an is_with_svg_plot option
    # plt.savefig(f"results/aero_structural_plot{config.kite_name}_idx22.svg", format="svg")
    # TODO: this here is to create an animation, from each simulated frame
    plt.savefig(f"{config.output_path}/animation/plot_iteration_{it_number}.png")
    if it_number != 00:  # TODO: closing when making animation
        plt.close()

    # plt.savefig(
    #     f"{config.output_path}/torque_paper/images/turning.svg", format="svg"
    # )

    return


# %% Printing out solution


def get_slack_values(
    pos, points_ini, ci, cj, ci_kite, cj_kite, springL, springL_kite, TE_extension, u_p
):  ##TODO: imported by functions_print
    slack_index = []
    for i in range(0, len(ci)):
        if i > 21:
            # LE bridles
            slack_index.append("#00000080")
        else:
            # TE bridles
            slack_index.append("#00000080")

    dL_kite_stetch_lst, dL_kite_slack_lst = np.zeros(len(ci_kite)), np.zeros(
        len(ci_kite)
    )  # Initializing, with additional entry for when it remains 0
    for i in range(0, len(ci_kite)):  # i loops over each (bridle or kite) line
        sep_vec = (
            pos[ci_kite[i]] - pos[cj_kite[i]]
        )  # Vector(ci --> cj) separating the points, indicating direction
        sep = np.linalg.norm(
            sep_vec
        )  # Absolute magnitude of the vector (works both ways), indicating
        dL_perc = 100 * ((sep - springL_kite[i]) / sep)
        if dL_perc > 0:
            dL_kite_stetch_lst[i] = dL_perc  # Append when bridle is longer
        else:
            dL_kite_slack_lst[i] = dL_perc

    dL_bridle_slack_lst, dL_bridle_stretch_lst = np.zeros(len(ci_kite)), np.zeros(
        len(ci_kite)
    )  # Initializing, with additional entry for when it remains 0
    for i in range(0, len(ci)):  # i loops over each (bridle or kite) line
        if i == 4 or i == 5:  # An attempt at making the pulley line equal colour
            sep_vec_1 = (
                pos[ci[i]] - pos[cj[i]]
            )  # Vector(ci --> cj) separating the points, indicating direction
            sep_1 = np.linalg.norm(
                sep_vec_1
            )  # Absolute magnitude of the vector (works both ways), indicating strength
            unit_vector_1 = sep_vec_1 / sep_1  # Define the unit_vector (ci --> cj)

            i_n = i * 2 - 2  # 4 --> 6, 6 --> 8
            sep_vec_2 = (
                pos[ci[i_n]] - pos[cj[i_n]]
            )  # Vector(ci --> cj) separating the points, indicating direction
            sep_2 = np.linalg.norm(
                sep_vec_2
            )  # Absolute magnitude of the vector (works both ways), indicating strength
            unit_vector_2 = sep_vec_2 / sep_2  # Define the unit_vector (ci --> cj)

            dL = ((sep_1 + sep_2) - (springL[i] + springL[i_n])) / (
                springL[i] + springL[i_n]
            )  # SpringL is defined on a range(0,len(ci)) loop (works both ways)
            dL_perc = 100 * dL

        elif i != 6 and i != 8:
            sep_vec = (
                pos[ci[i]] - pos[cj[i]]
            )  # Vector(ci --> cj) separating the points, indicating direction
            sep = np.linalg.norm(
                sep_vec
            )  # Absolute magnitude of the vector (works both ways), indicating strength
            dL_perc = 100 * ((sep - springL[i]) / springL[i])
        elif i == 6:
            dL_perc = dL_bridle_stretch_lst[3]
        elif i == 8:
            dL_perc = dL_bridle_stretch_lst[5]

        if dL_perc > 0:
            dL_bridle_stretch_lst[i] = dL_perc  # Append when bridle is longer
        else:
            dL_bridle_slack_lst[i] = dL_perc
            if dL_perc > -2.5:
                slack_index[i] = "#1CEFCC"
            elif dL_perc <= -2.5 and dL_perc > -5:
                slack_index[i] = "#FFBD19"
            elif dL_perc < -5:
                slack_index[i] = "#EF1CEC"

    return (
        pos,
        slack_index,
        dL_bridle_stretch_lst,
        dL_bridle_slack_lst,
        dL_kite_stetch_lst,
        dL_kite_slack_lst,
    )


# %% structural plots?


##TODO: used
def distance(A, B):
    return np.linalg.norm(B - A)


##TODO: used
def plot_kiteforces(pos, ci, cj, ci_kite, cj_kite, lift_force, ax, elev):
    """Make axes of 3D plot have equal scale so that spheres appear as spheres,
      is is one possible solution to Matplotlib's
      ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

      Input
        ax: a matplotlib axis, e.g., as output from plt.gca().
    cubes as cubes, etc..  Th
    """
    setattr(Axes3D, "arrow3D", _arrow3D)
    width = 0.5
    pos = pos / 10
    for i in range(0, len(ci_kite)):
        ax.plot(
            [pos[ci_kite[i], 0], pos[cj_kite[i], 0]],
            [pos[ci_kite[i], 1], pos[cj_kite[i], 1]],
            [pos[ci_kite[i], 2], pos[cj_kite[i], 2]],
            color="black",
            linewidth=width,
        )  # ,legend='kite')
    width = 1
    # for i in range(0,len(lift_force)):
    #     ax.plot([pos[i, 0], pos[i, 0]+lift_force[i,0]], [pos[i, 1], pos[i, 1]+lift_force[i,1]], [pos[i, 2], pos[i, 2]+lift_force[i,2]],color='green',linewidth = width)#,legend='kite')

    for i in range(0, len(lift_force)):
        a = pos[i]
        b = lift_force[i]
        ax.arrow3D(
            a[0],
            a[1],
            a[2],
            b[0],
            b[1],
            b[2],
            mutation_scale=5,
            linewidth=0.5,
            arrowstyle="-|>",
            fc="green",
            ec="green",
        )

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()
    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)
    ax.view_init(elev=elev, azim=45)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5 * max([x_range, y_range, z_range])

    # ax.view_init(elev=2000, azim=None)
    # if dist is not None:
    #     ax.dist = dist
    # if elev is not None:
    #     ax.elev = elev

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])
    return


##TODO: used
def plot_kite_deformation(
    pos, ci, cj, ci_kite, cj_kite, plates, col_kite, ax, line_boolean, elev, tube_idx
):
    """Make axes of 3D plot have equal scale so that spheres appear as spheres,
      is is one possible solution to Matplotlib's
      ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

      Input
        ax: a matplotlib axis, e.g., as output from plt.gca().
    cubes as cubes, etc..  Th
    """
    TE = [2, 8, 14, 20, 26, 32, 38, 44, 50]
    width = 0.7
    msize = 5

    # Fill canopy
    for i in range(len(plates)):
        x = np.array(
            [
                pos[plates[i][0], 0],
                pos[plates[i][1], 0],
                pos[plates[i][2], 0],
                pos[plates[i][3], 0],
            ]
        )
        y = np.array(
            [
                pos[plates[i][0], 1],
                pos[plates[i][1], 1],
                pos[plates[i][2], 1],
                pos[plates[i][3], 1],
            ]
        )
        z = np.array(
            [
                pos[plates[i][0], 2],
                pos[plates[i][1], 2],
                pos[plates[i][2], 2],
                pos[plates[i][3], 2],
            ]
        )

        # Create a polygonal surface between the lines
        vertices = [(xi, yi, zi) for xi, yi, zi in zip(x, y, z)]
        # vertices += [(xi, yi, zi) for xi, yi, zi in zip(x[::-1], z[::-1], np.zeros_like(x))]
        poly = Poly3DCollection([vertices], alpha=0.1, facecolor=col_kite)
        ax.add_collection3d(poly)

    # ax.scatter(pos[:,0], pos[:,1], pos[:,2], color=color)

    for i in range(0, len(ci_kite)):
        if i in tube_idx:
            ax.plot(
                [pos[ci_kite[i], 0], pos[cj_kite[i], 0]],
                [pos[ci_kite[i], 1], pos[cj_kite[i], 1]],
                [pos[ci_kite[i], 2], pos[cj_kite[i], 2]],
                color=col_kite,
                linewidth=width * 5,
                solid_capstyle="round",
            )
        elif i in TE:
            ax.plot(
                [pos[ci_kite[i], 0], pos[cj_kite[i], 0]],
                [pos[ci_kite[i], 1], pos[cj_kite[i], 1]],
                [pos[ci_kite[i], 2], pos[cj_kite[i], 2]],
                color=col_kite,
                linewidth=width,
                solid_capstyle="round",
            )

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()
    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = max([0.5 * x_range, 0.5 * y_range, 0.3 * z_range])

    # ax.view_init(elev=2000, azim=None)
    # if dist is not None:
    #     ax.dist = dist
    # if elev is not None:
    #     ax.elev = elev

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])
    return


##TODO: used
def plot_kite_matplotlib(
    pos,
    ci,
    cj,
    ci_kite,
    cj_kite,
    plates,
    slack_index,
    col_kite,
    ax,
    line_boolean,
    elev,
    tube_idx,
):
    """Make axes of 3D plot have equal scale so that spheres appear as spheres,
      is is one possible solution to Matplotlib's
      ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

      Input
        ax: a matplotlib axis, e.g., as output from plt.gca().
    cubes as cubes, etc..  Th
    """
    TE = [2, 8, 14, 20, 26, 32, 38, 44, 50]
    width = 0.7
    msize = 5

    # Fill canopy
    for i in range(len(plates)):
        x = np.array(
            [
                pos[plates[i][0], 0],
                pos[plates[i][1], 0],
                pos[plates[i][2], 0],
                pos[plates[i][3], 0],
            ]
        )
        y = np.array(
            [
                pos[plates[i][0], 1],
                pos[plates[i][1], 1],
                pos[plates[i][2], 1],
                pos[plates[i][3], 1],
            ]
        )
        z = np.array(
            [
                pos[plates[i][0], 2],
                pos[plates[i][1], 2],
                pos[plates[i][2], 2],
                pos[plates[i][3], 2],
            ]
        )

        # Create a polygonal surface between the lines
        vertices = [(xi, yi, zi) for xi, yi, zi in zip(x, y, z)]
        # vertices += [(xi, yi, zi) for xi, yi, zi in zip(x[::-1], z[::-1], np.zeros_like(x))]
        poly = Poly3DCollection([vertices], alpha=0.1, facecolor="black")
        ax.add_collection3d(poly)
    plt.rcParams.update({"font.size": 11})
    plt.rcParams.update(
        {
            "text.usetex": True,
            "font.family": "serif",
            "font.serif": ["Computer Modern Roman"],
        }
    )
    # ax.scatter(pos[:,0], pos[:,1], pos[:,2], color=color)
    if line_boolean == True:
        for i in range(0, len(ci)):
            color = slack_index[i]
            ax.plot(
                [pos[ci[i], 0], pos[cj[i], 0]],
                [pos[ci[i], 1], pos[cj[i], 1]],
                [pos[ci[i], 2], pos[cj[i], 2]],
                color=color,
                linewidth=width,
            )

    for i in range(0, len(ci_kite)):
        if i in tube_idx:
            ax.plot(
                [pos[ci_kite[i], 0], pos[cj_kite[i], 0]],
                [pos[ci_kite[i], 1], pos[cj_kite[i], 1]],
                [pos[ci_kite[i], 2], pos[cj_kite[i], 2]],
                color=col_kite,
                linewidth=width * 5,
                solid_capstyle="round",
            )
        elif i in TE:
            ax.plot(
                [pos[ci_kite[i], 0], pos[cj_kite[i], 0]],
                [pos[ci_kite[i], 1], pos[cj_kite[i], 1]],
                [pos[ci_kite[i], 2], pos[cj_kite[i], 2]],
                color=col_kite,
                linewidth=width,
                solid_capstyle="round",
            )

    legend_elements = [
        Line2D([0], [0], color="#00000080", lw=0.5),
        Line2D([0], [0], color="#1CEFCC", lw=0.5),
        Line2D([0], [0], color="#FFBD19", lw=0.5),
        # Line2D([0], [0], color='#EF1CEC', lw=0.5)
    ]
    ax.legend(
        legend_elements,
        ["no slack", "slack $< 2.5\%$", "$2.5\% <$ slack $< 5\%$", "slack $> 5\%$"],
        frameon=False,
        loc="center right",
    )

    # Plot tapes
    id1 = 21
    id2 = 27
    ax.plot(
        [pos[id1, 0], pos[id2, 0]],
        [pos[id1, 1], pos[id2, 1]],
        [pos[id1, 2], pos[id2, 2]],
        color="#000000",
        linewidth=width * 3,
    )
    id2 = 22
    ax.plot(
        [pos[id1, 0], pos[id2, 0]],
        [pos[id1, 1], pos[id2, 1]],
        [pos[id1, 2], pos[id2, 2]],
        color="#000000",
        linewidth=width * 3,
    )
    id2 = 23
    ax.plot(
        [pos[id1, 0], pos[id2, 0]],
        [pos[id1, 1], pos[id2, 1]],
        [pos[id1, 2], pos[id2, 2]],
        color="#000000",
        linewidth=width * 3,
    )

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()
    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5 * max([x_range, y_range, z_range])

    # ax.view_init(elev=2000, azim=None)
    # if dist is not None:
    #     ax.dist = dist
    # if elev is not None:
    #     ax.elev = elev

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])
    return


##TODO: used
def plot_kite_matplotlib_pretty(
    pos, ci, cj, ax, line_boolean, col_kite, elev, tube_idx, ci_kite, cj_kite, plates
):
    """Make axes of 3D plot have equal scale so that spheres appear as spheres,
      is is one possible solution to Matplotlib's
      ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

      Input
        ax: a matplotlib axis, e.g., as output from plt.gca().
    cubes as cubes, etc..  Th
    """

    TE = [2, 8, 14, 20, 26, 32, 38, 44, 50]
    width = 0.7
    msize = 5

    # Fill canopy
    for i in range(len(plates)):
        x = np.array(
            [
                pos[plates[i][0], 0],
                pos[plates[i][1], 0],
                pos[plates[i][2], 0],
                pos[plates[i][3], 0],
            ]
        )
        y = np.array(
            [
                pos[plates[i][0], 1],
                pos[plates[i][1], 1],
                pos[plates[i][2], 1],
                pos[plates[i][3], 1],
            ]
        )
        z = np.array(
            [
                pos[plates[i][0], 2],
                pos[plates[i][1], 2],
                pos[plates[i][2], 2],
                pos[plates[i][3], 2],
            ]
        )

        # Create a polygonal surface between the lines
        vertices = [(xi, yi, zi) for xi, yi, zi in zip(x, y, z)]
        # vertices += [(xi, yi, zi) for xi, yi, zi in zip(x[::-1], z[::-1], np.zeros_like(x))]
        poly = Poly3DCollection([vertices], alpha=0.1, facecolor="black")
        ax.add_collection3d(poly)

    ip = 1
    for i in range(0, len(ci_kite)):
        if i in tube_idx:
            ax.plot(
                [pos[ci_kite[i], 0], pos[cj_kite[i], 0]],
                [pos[ci_kite[i], 1], pos[cj_kite[i], 1]],
                [pos[ci_kite[i], 2], pos[cj_kite[i], 2]],
                color=col_kite,
                linewidth=width * 6,
            )
        elif i in TE:
            ip += 1
            ax.plot(
                [pos[ci_kite[i], 0], pos[cj_kite[i], 0]],
                [pos[ci_kite[i], 1], pos[cj_kite[i], 1]],
                [pos[ci_kite[i], 2], pos[cj_kite[i], 2]],
                color=col_kite,
                linewidth=width,
            )

    for i in range(0, len(ci)):
        if i > 21:
            ax.plot(
                [pos[ci[i], 0], pos[cj[i], 0]],
                [pos[ci[i], 1], pos[cj[i], 1]],
                [pos[ci[i], 2], pos[cj[i], 2]],
                color="#23B521",
                linewidth=width,
                marker=".",
                markersize=msize,
            )
        else:
            ax.plot(
                [pos[ci[i], 0], pos[cj[i], 0]],
                [pos[ci[i], 1], pos[cj[i], 1]],
                [pos[ci[i], 2], pos[cj[i], 2]],
                color="#0D23C2",
                linewidth=width,
                marker=".",
                markersize=msize,
            )

    # Plot KCU
    ax.plot(
        pos[21, 0],
        pos[21, 1],
        pos[21, 2],
        color="#000000",
        linewidth=width,
        marker="v",
        markersize=msize * 3,
    )
    # Plot pulleys
    ax.plot(
        pos[24, 0],
        pos[24, 1],
        pos[24, 2],
        color="orange",
        linewidth=width,
        marker="P",
        markersize=msize * 1.5,
    )
    ax.plot(
        pos[28, 0],
        pos[28, 1],
        pos[28, 2],
        color="orange",
        linewidth=width,
        marker="P",
        markersize=msize * 1.5,
    )
    # Plot tapes
    id1 = 21
    id2 = 27
    ax.plot(
        [pos[id1, 0], pos[id2, 0]],
        [pos[id1, 1], pos[id2, 1]],
        [pos[id1, 2], pos[id2, 2]],
        color="#0D23C2",
        linewidth=width * 3,
    )
    id2 = 22
    ax.plot(
        [pos[id1, 0], pos[id2, 0]],
        [pos[id1, 1], pos[id2, 1]],
        [pos[id1, 2], pos[id2, 2]],
        color="#0D23C2",
        linewidth=width * 3,
    )
    id2 = 23
    ax.plot(
        [pos[id1, 0], pos[id2, 0]],
        [pos[id1, 1], pos[id2, 1]],
        [pos[id1, 2], pos[id2, 2]],
        color="#0D23C2",
        linewidth=width * 3,
    )
    # Plot KCU
    ax.plot(
        pos[21, 0],
        pos[21, 1],
        pos[21, 2],
        color="#000000",
        linewidth=width,
        marker="v",
        markersize=msize * 3,
    )
    # Plot point mass particles
    for i in range(0, len(ci)):
        ax.plot(
            pos[ci[i], 0],
            pos[ci[i], 1],
            pos[ci[i], 2],
            color="#A605CD",
            linewidth=width,
            marker=".",
            markersize=msize,
        )
        ax.plot(
            pos[cj[i], 0],
            pos[cj[i], 1],
            pos[cj[i], 2],
            color="#A605CD",
            linewidth=width,
            marker=".",
            markersize=msize,
        )

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()
    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5 * max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])
    from matplotlib.lines import Line2D

    legend_elements = [
        Line2D([0], [0], color=col_kite, lw=width * 5),
        Line2D([0], [0], color=col_kite, lw=width),
        Line2D([0], [0], color="#23B521", lw=width),
        Line2D([0], [0], color="#0D23C2", lw=width),
        Line2D([0], [0], color="#0D23C2", lw=width * 3),
        Line2D([0], [0], marker="v", color="w", markerfacecolor="black", markersize=15),
        Line2D(
            [0], [0], marker=".", color="w", markerfacecolor="#A605CD", markersize=10
        ),
        Line2D(
            [0], [0], marker="P", color="w", markerfacecolor="orange", markersize=10
        ),
    ]

    ax.legend(
        legend_elements,
        [
            "Inflatable tubes",
            "Trailing edge",
            "Power lines",
            "Steering lines",
            "Steering tapes",
            "KCU",
            "Point mass particles",
            "Pulleys",
        ],
        frameon=False,
        loc="center right",
    )

    plt.tight_layout()
    return


def _arrow3D(ax, x, y, z, dx, dy, dz, *args, **kwargs):
    """Add an 3d arrow to an `Axes3D` instance."""

    arrow = Arrow3D(x, y, z, dx, dy, dz, *args, **kwargs)
    ax.add_artist(arrow)


class Arrow3D(FancyArrowPatch):
    def __init__(self, x, y, z, dx, dy, dz, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._xyz = (x, y, z)
        self._dxdydz = (dx, dy, dz)

    def draw(self, renderer):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        super().draw(renderer)

    def do_3d_projection(self, renderer=None):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))

        return np.min(zs)


def plot_alex_master_thesis(position, last_index, connectivity_matrix, params):
    plt.close("all")
    n = params["n"]

    # plotting & graph configuration
    # Data from layout after 1 iteration step
    X = []
    Y = []
    Z = []
    for i in range(n):
        X.append(position[f"x{i + 1}"].iloc[0])
        Y.append(position[f"y{i + 1}"].iloc[0])
        Z.append(position[f"z{i + 1}"].iloc[0])

    fig = plt.figure(figsize=plt.figaspect(0.5))
    ax = fig.add_subplot(1, 2, 2, projection="3d")
    ax2 = fig.add_subplot(1, 2, 1, projection="3d")

    # ensuring the axis are scaled properly
    ax.set_xlim(-5, 5)
    ax.set_ylim(-5, 5)
    ax.set_zlim(0, 10)
    ax2.set_xlim(-5, 5)
    ax2.set_ylim(-5, 5)
    ax2.set_zlim(0, 10)

    # data from final timestep
    X_f = []
    Y_f = []
    Z_f = []
    for i in range(n):
        X_f.append(position[f"x{i + 1}"].iloc[last_index])
        Y_f.append(position[f"y{i + 1}"].iloc[last_index])
        Z_f.append(position[f"z{i + 1}"].iloc[last_index])

    # plot inital layout
    ax.scatter(X, Y, Z, c="red")
    for indices in connectivity_matrix:
        ax.plot(
            [X[indices[0]], X[indices[1]]],
            [Y[indices[0]], Y[indices[1]]],
            [Z[indices[0]], Z[indices[1]]],
            color="black",
        )

    # plot final found shape
    ax2.scatter(X_f, Y_f, Z_f, c="red")
    for indices in connectivity_matrix:
        ax2.plot(
            [X_f[indices[0]], X_f[indices[1]]],
            [Y_f[indices[0]], Y_f[indices[1]]],
            [Z_f[indices[0]], Z_f[indices[1]]],
            color="black",
        )

    plt.show()
    return


# #%% stack part 1

# import numpy as np
# import plotly.graph_objects as go

# # define function
# def get_3D_plot_stack(x,y,z,connections,elongation_values):

#     data_list = [] #initializing empty list

#     ### Points, gathering data and appending to data_list
#     for i,(xi,yi,zi) in enumerate(zip(x,y,z)): # looping through each point_list
#         data_list.append(go.Scatter3d(x=[xi], y=[yi], z=[zi],
#                                       mode='markers',
#                                       marker=dict(color='black',size=2),
#                                       showlegend= False))

#     ### Lines, gathering data and appending to data_list
#     x_conn, y_conn, z_conn = np.empty((len(connections),2)), np.empty((len(connections),2)), np.empty((len(connections),2))
#     for i,(conn,elongation_i) in enumerate(zip(connections,elongation_values)):
#         xi_conn = [x[conn[0]], x[conn[1]]]
#         yi_conn = [y[conn[0]], y[conn[1]]]
#         zi_conn = [z[conn[0]], z[conn[1]]]

#         # storing data
#         x_conn[i], y_conn[i], z_conn[i] = xi_conn, yi_conn, zi_conn

#         data_list.append(go.Scatter3d(
#                                     x=xi_conn, y=yi_conn, z=zi_conn,
#                                     mode='lines',
#                                     line=go.scatter3d.Line(
#                                         width = 4,
#                                         color=elongation_i,
#                                         colorscale='Viridis',
#                                         showscale=True,
#                                     ),
#                                 showlegend=False
#                             ))

#     ### Create figure
#     fig = go.FigureWidget(data=data_list)

#     fig.show()
#     return()

# x,y,z = np.random.random_sample((3,10)) # random points
# connections = np.array([[0,1],[9,2],[2,3],[5,7],[6,8],[2,8],[1,2],[4,5]]) # line connections
# elongation_values = np.random.random_sample((len(connections))) # random colors

# get_3D_plot_stack(x,y,z,connections,elongation_values)


# #%% stack part 2

# import numpy as np
# import plotly.graph_objects as go

# def get_3D_plot_stack(x,y,z,connections,elongation_values):

#     data_list = [] #initializing an empty list

#     ### Points, gathering data and appending to data_list
#     for i,(xi,yi,zi) in enumerate(zip(x,y,z)): # looping through each point_list
#         data_list.append(go.Scatter3d(x=[xi], y=[yi], z=[zi],
#                                       mode='markers',
#                                       marker=dict(color='black',size=2),
#                                       showlegend= False))

#     ### Lines, gathering data and appending to data_list
#     x_conn, y_conn, z_conn = np.empty((len(connections),2)), np.empty((len(connections),2)), np.empty((len(connections),2))
#     for i,(conn,elongation_i) in enumerate(zip(connections,elongation_values)):
#         xi_conn = [x[conn[0]], x[conn[1]]]
#         yi_conn = [y[conn[0]], y[conn[1]]]
#         zi_conn = [z[conn[0]], z[conn[1]]]

#         # storing data
#         x_conn[i], y_conn[i], z_conn[i] = xi_conn, yi_conn, zi_conn

#         data_list.append(go.Scatter3d(
#                                     x=xi_conn, y=yi_conn, z=zi_conn,
#                                     mode='lines',
#                                     line=go.scatter3d.Line(
#                                         width = 4,
#                                         color=elongation_i,
#                                         colorscale='Viridis',
#                                         showscale=False,
#                                     ),
#                                 showlegend=False
#                             ))

#     # getting the colorbar once
#     line_trace_all = go.Scatter3d( x=x_conn, y=y_conn, z=z_conn,
#                                     mode='lines',
#                                     line=go.scatter3d.Line(
#                                         color=elongation_values,
#                                         colorscale='Viridis',
#                                         showscale=True),
#                                     showlegend=False)

#     data_list.append(line_trace_all)

#     ### Create figure
#     fig = go.FigureWidget(data=data_list)

#     fig.show()
#     return()

# x,y,z = np.random.random_sample((3,10))
# connections = np.array([[0,1],[9,2],[2,3],[5,7],[6,8],[2,8],[1,2],[4,5]]) #random
# elongation_values = np.random.random_sample((len(connections)))

# get_3D_plot_stack(x,y,z,connections,elongation_values)

# #%% for stack overflow

# import numpy as np
# import plotly.graph_objects as go

# def gen_color_bar(line_trace): #from GitHub
#     """
#     Generates a trace which shows a colorbar based on a line plot.

#     Relevant issue: https://github.com/plotly/plotly.py/issues/1085
#     """
#     return go.Scatter3d(
#         x=line_trace.x, y=line_trace.y, z=line_trace.z,
#         mode="markers",
#         marker=go.scatter3d.Marker(
#             color=line_trace.line.color,
#             colorscale=line_trace.line.to_plotly_json()["colorscale"], # Waiting on https://github.com/plotly/plotly.py/issues/1087
#             showscale=line_trace.line.showscale,
#             opacity=0.00000000000001 # Make invisible, visible=False disables color bar
#         ),
#         hoverinfo="none",
#         showlegend=False
#     )

# def get_3D_plot_stack(x,y,z,connections,elongation_values):

#     data_list = [] #initializing an empty list

#     ### Points, gathering data and appending to data_list
#     for i,(xi,yi,zi) in enumerate(zip(x,y,z)): # looping through each point_list
#         data_list.append(go.Scatter3d(x=[xi], y=[yi], z=[zi],
#                                       mode='markers',
#                                       marker=dict(color='black',size=2),
#                                       showlegend= False))

#     ### Lines, gathering data and appending to data_list
#     x_conn, y_conn, z_conn = np.empty((len(connections),2)), np.empty((len(connections),2)), np.empty((len(connections),2))
#     for i,(conn,elongation_i) in enumerate(zip(connections,elongation_values)):
#         xi_conn = [x[conn[0]], x[conn[1]]]
#         yi_conn = [y[conn[0]], y[conn[1]]]
#         zi_conn = [z[conn[0]], z[conn[1]]]

#         # storing data
#         x_conn[i], y_conn[i], z_conn[i] = xi_conn, yi_conn, zi_conn

#         data_list.append(go.Scatter3d(
#                                     x=xi_conn, y=yi_conn, z=zi_conn,
#                                     mode='lines',
#                                     line=go.scatter3d.Line(
#                                         width = 4,
#                                         color=elongation_i,
#                                         colorscale='Viridis',
#                                         showscale=False,
#                                     ),
#                                 showlegend=False
#                             ))

#     # getting the colorbar
#     line_trace_list = go.Scatter3d( x=x_conn, y=y_conn, z=z_conn,
#                                     mode='lines',
#                                     line=go.scatter3d.Line(
#                                         color=elongation_values,
#                                         colorscale='Viridis',
#                                         showscale=True),
#                                     showlegend=False)

#     data_list.append(line_trace_list)

#     ### Create figure
#     fig = go.FigureWidget(data=data_list)

#     fig.show()
#     return()

# x,y,z = np.random.random_sample((3,10))
# connections = np.array([[0,1],[9,2],[2,3],[5,7],[6,8],[2,8],[1,2],[4,5]]) #random
# elongation_values = np.random.random_sample((len(connections)))

# get_3D_plot_stack(x,y,z,connections,elongation_values)
