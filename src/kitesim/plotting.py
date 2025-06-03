import numpy as np
import matplotlib.pyplot as plt

import matplotlib.pyplot as plt
from typing import Optional
from cycler import cycler
from matplotlib.widgets import Slider

PALETTE = {
    "Black": "#000000",
    "Orange": "#E69F00",
    "Sky Blue": "#56B4E9",
    "Bluish Green": "#009E73",
    "Yellow": "#F0E442",
    "Blue": "#0072B2",
    "Vermillion": "#D55E00",
    "Reddish Purple": "#CC79A7",
}


def set_plot_style():
    """
    Set the default style for plots using LaTeX and custom color palette.

    Tips:
    - If you specify colors, they will still be used.
    - If you want to change the axis margins:
        1. try with ax.xlim and ax.ylim
        2. try by changing the 'axes.autolimit_mode' parameter to data
    - more?
    """

    # Define the color palette as a list of colors
    color_cycle = [
        PALETTE["Black"],
        PALETTE["Orange"],
        PALETTE["Sky Blue"],
        PALETTE["Bluish Green"],
        PALETTE["Yellow"],
        PALETTE["Blue"],
        PALETTE["Vermillion"],
        PALETTE["Reddish Purple"],
    ]

    # Apply Seaborn style and custom settings
    # plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams.update(
        {
            "text.usetex": True,
            "font.family": "serif",
            "font.serif": ["Computer Modern Roman"],
            ## Axes settings
            "axes.titlesize": 15,
            "axes.labelsize": 13,
            "axes.linewidth": 1.0,
            "axes.edgecolor": "#C5C5C5",
            "axes.labelcolor": "black",
            "axes.autolimit_mode": "round_numbers",
            "axes.xmargin": 0,  # Remove extra margin
            "axes.ymargin": 0,  # Remove extra margin
            ## Grid settings
            "axes.grid": True,
            "axes.grid.axis": "both",
            "grid.alpha": 0.5,
            "grid.color": "#C5C5C5",
            "grid.linestyle": "-",
            "grid.linewidth": 1.0,
            ## Line settings
            "lines.linewidth": 1,
            "lines.markersize": 6,
            # "lines.color": "grey",,
            "figure.titlesize": 15,
            "pgf.texsystem": "pdflatex",  # Use pdflatex
            "pgf.rcfonts": False,
            "figure.figsize": (15, 5),  # Default figure size
            "axes.prop_cycle": cycler(
                "color", color_cycle
            ),  # Set the custom color cycle
            ## tick settings
            "xtick.color": "#C5C5C5",
            "ytick.color": "#C5C5C5",
            "xtick.labelcolor": "black",
            "ytick.labelcolor": "black",
            "xtick.labelsize": 13,
            "ytick.labelsize": 13,
            "xtick.top": True,  # Show ticks on both sides
            "xtick.bottom": True,
            "ytick.left": True,
            "ytick.right": True,
            "xtick.direction": "in",  # Direction for x-axis ticks
            "ytick.direction": "in",  # Direction for y-axis ticks
            ## legend settings
            "legend.fontsize": 15,
        }
    )


def main(
    struc_nodes,
    kite_connectivity,
    rest_lengths,
    struc_nodes_initial=None,
    f_ext=None,
    title="PSM State",
    body_aero=None,
    chord_length=2.6,  # new argument for scaling force vectors
    is_with_node_indices=False,
    lightred_color="#FF5F5F",
):
    """
    Plot the current (and optionally initial) structure state in 3D.

    Args:
        struc_nodes (np.ndarray): Current node positions (n_nodes, 3).
        kite_connectivity (array-like): List/array of [i, j, ...] giving spring connections.
        struc_nodes_initial (np.ndarray, optional): Initial node positions (n_nodes, 3).
        f_ext (np.ndarray or None): Optional external forces, shape (n_nodes, 3) or flat.
        title (str): Figure title.
        chord_length (float): Maximum length for force vectors (for scaling).

    Returns:
        None. Displays a 3D plot.
    """
    kite_connectivity = np.array(kite_connectivity)  # Ensure numeric array

    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    # Try to maximize the plot window (works for Qt backends)
    try:
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
    except Exception:
        # Fallback: set a large figure size
        fig.set_size_inches(16, 9, forward=True)

    # Plot initial state if provided (single color for reference)
    if struc_nodes_initial is not None:
        ax.scatter(
            *(struc_nodes_initial.T), color="black", marker="o", s=10, label="Initial"
        )
        # Draw initial lines in gray
        for idx, (i, j, *rest) in enumerate(kite_connectivity):
            p1, p2 = struc_nodes_initial[i], struc_nodes_initial[j]
            ax.plot(
                [p1[0], p2[0]],
                [p1[1], p2[1]],
                [p1[2], p2[2]],
                color="gray",
                linewidth=1,
                alpha=0.5,
            )

    # Plot current state
    ax.scatter(*(struc_nodes.T), color="black", marker="o", s=10, label="Particles")
    # Compute normalized elongation for each connection
    cmap = plt.get_cmap("viridis")  # Changed from "coolwarm" to "viridis"
    norm_values = []
    for idx, (i, j, *rest) in enumerate(kite_connectivity):
        p1, p2 = struc_nodes[i], struc_nodes[j]
        curr_length = np.linalg.norm(p2 - p1)
        # Change: multiply by 100 for percent change
        norm_val = 100 * (curr_length - rest_lengths[idx]) / rest_lengths[idx]
        norm_values.append(norm_val)
    # Normalize for color mapping
    vmin, vmax = min(norm_values), max(norm_values)
    for idx, (i, j, *rest) in enumerate(kite_connectivity):
        p1, p2 = struc_nodes[i], struc_nodes[j]
        color = cmap((norm_values[idx] - vmin) / (vmax - vmin) if vmax > vmin else 0.5)
        ax.plot(
            [p1[0], p2[0]],
            [p1[1], p2[1]],
            [p1[2], p2[2]],
            color=color,
            linewidth=2,  # thicker lines
        )

    # Optionally plot external forces
    if f_ext is not None:
        arr = np.array(f_ext)
        if arr.ndim == 1:
            arr = arr.reshape(-1, 3)
        # Scale all vectors so the longest has length chord_length/5 (larger arrows)
        norms = np.linalg.norm(arr, axis=1)
        max_norm = np.max(norms) if np.max(norms) > 0 else 1.0
        scale = (chord_length / 2) / max_norm
        arr_scaled = arr * scale
        ax.quiver(
            struc_nodes[:, 0],
            struc_nodes[:, 1],
            struc_nodes[:, 2],
            arr_scaled[:, 0],
            arr_scaled[:, 1],
            arr_scaled[:, 2],
            length=1,
            normalize=False,
            color=lightred_color,
        )

    if is_with_node_indices:
        # Annotate node indices
        for i, pos in enumerate(struc_nodes):
            ax.text(pos[0], pos[1], pos[2], str(i), color="black", fontsize=8)
        if struc_nodes_initial is not None:
            for i, pos in enumerate(struc_nodes_initial):
                ax.text(pos[0], pos[1], pos[2], str(i), color="blue", fontsize=8)

    # If aero mesh nodes are provided, plot them
    if body_aero is not None:
        aero_mesh_nodes = []
        for panel in body_aero.panels:
            for cp in panel.corner_points:
                aero_mesh_nodes.append(cp)
        # make unique
        aero_mesh_nodes = np.unique(np.array(aero_mesh_nodes), axis=0)
        ax.scatter(
            *(aero_mesh_nodes.T),
            color=lightred_color,
            marker="^",
            s=10,
            label="Aero Mesh",
        )

    # Set aspect ratio to equal
    all_pts = (
        struc_nodes
        if struc_nodes_initial is None
        else np.vstack((struc_nodes, struc_nodes_initial))
    )
    bb = all_pts.max(axis=0) - all_pts.min(axis=0)
    ax.set_box_aspect(bb)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title(title)
    # Move legend outside plot area to the right
    ax.legend(loc="center right", bbox_to_anchor=(1.05, 0.5))
    # Always add colorbar for current state
    if len(norm_values) > 0:
        import matplotlib as mpl

        sm = mpl.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, shrink=0.7, pad=0.1)
        # Change: update colorbar label to percent
        cbar.set_label("Normalized rest length change (\%)")
        # Add contracted (bottom) and elongated (top) labels
        cbar.ax.text(
            1.5,
            -0.1,
            "contracted",
            va="center",
            ha="left",
            fontsize=11,
            rotation=0,
            color="black",
            transform=cbar.ax.transAxes,
        )
        cbar.ax.text(
            1.5,
            1.1,
            "elongated",
            va="center",
            ha="left",
            fontsize=11,
            rotation=0,
            color="black",
            transform=cbar.ax.transAxes,
        )
    plt.show()


def interactive_plot(
    tracking_data,
    kite_connectivity,
    rest_lengths,
    f_ext=None,
    title="PSM State",
    chord_length=2.6,
    lightred_color="#FF5F5F",
    elev=15,
    azim=-120,
    t_per_step=0.5,
):
    """
    Interactive plot for tracked positions over iterations.

    Args:
        tracking_data (dict): Output from tracking.setup_tracking_arrays.
        kite_connectivity (array-like): Connectivity for plotting.
        rest_lengths (np.ndarray): Rest lengths for each connection.
        f_ext (np.ndarray or None): External forces (nt, n_nodes, 3) or None.
        title (str): Plot title.
        chord_length (float): For force vector scaling.
        lightred_color (str): Color for forces/aero mesh.
        elev (float): Initial elevation angle for 3D plot.
        azim (float): Initial azimuth angle for 3D plot.

    Returns:
        None. Displays an interactive plot.
    """
    import matplotlib.pyplot as plt

    positions = tracking_data["positions"]  # shape (nt, n_nodes, 3)
    nt = positions.shape[0]
    kite_conn = np.array(kite_connectivity)
    rest_lengths = np.array(rest_lengths)

    if f_ext is not None and f_ext.ndim == 3:
        f_ext_all = f_ext
    else:
        f_ext_all = None

    # --- Compute maximum elongation over all time steps and connections ---
    max_abs_elong = 0.0
    for t in range(nt):
        struc_nodes = positions[t]
        elong = []
        for idx, (i, j, *_) in enumerate(kite_conn):
            p1, p2 = struc_nodes[i], struc_nodes[j]
            curr_length = np.linalg.norm(p2 - p1)
            percent_change = 100 * (curr_length - rest_lengths[idx]) / rest_lengths[idx]
            elong.append(percent_change)
        max_abs_elong = max(max_abs_elong, np.max(np.abs(elong)))
    normalized_elongation_limit = max_abs_elong * 1.01  # add 1% margin

    fig = plt.figure(figsize=(12, 7))
    from matplotlib.gridspec import GridSpec

    gs = GridSpec(
        2, 2, width_ratios=[20, 1], height_ratios=[20, 1], hspace=0.25, wspace=0.15
    )
    ax = fig.add_subplot(gs[0, 0], projection="3d")
    ax.view_init(elev=elev, azim=azim)

    vmin, vmax = -normalized_elongation_limit, normalized_elongation_limit
    cmap = plt.get_cmap("viridis")
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    idx0 = 0

    def draw_state(idx):
        ax.cla()
        ax.view_init(elev=elev, azim=azim)
        struc_nodes = positions[idx]
        ax.scatter(*(struc_nodes.T), color="black", marker="o", s=10, label="Particles")
        norm_values = []
        for idx2, (i, j, *rest) in enumerate(kite_conn):
            p1, p2 = struc_nodes[i], struc_nodes[j]
            curr_length = np.linalg.norm(p2 - p1)
            norm_val = 100 * (curr_length - rest_lengths[idx2]) / rest_lengths[idx2]
            norm_values.append(norm_val)
        for idx2, (i, j, *rest) in enumerate(kite_conn):
            p1, p2 = struc_nodes[i], struc_nodes[j]
            color = cmap(norm(norm_values[idx2]))
            ax.plot(
                [p1[0], p2[0]],
                [p1[1], p2[1]],
                [p1[2], p2[2]],
                color=color,
                linewidth=2,
            )
        if f_ext_all is not None:
            arr = np.array(f_ext_all[idx])
            if arr.ndim == 1:
                arr = arr.reshape(-1, 3)
            norms = np.linalg.norm(arr, axis=1)
            max_norm = np.max(norms) if np.max(norms) > 0 else 1.0
            scale = (chord_length / 2) / max_norm
            arr_scaled = arr * scale
            ax.quiver(
                struc_nodes[:, 0],
                struc_nodes[:, 1],
                struc_nodes[:, 2],
                arr_scaled[:, 0],
                arr_scaled[:, 1],
                arr_scaled[:, 2],
                length=1,
                normalize=False,
                color=lightred_color,
            )
        all_pts = struc_nodes
        bb = all_pts.max(axis=0) - all_pts.min(axis=0)
        ax.set_box_aspect(bb)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.set_title(f"{title} (step {idx+1}/{nt})")

    draw_state(idx0)

    import matplotlib as mpl

    ax_cbar = fig.add_subplot(gs[0, 1])
    cbar = plt.colorbar(sm, cax=ax_cbar, orientation="vertical")
    cbar.set_label("Normalized rest length change (\%)")
    cbar.ax.text(
        1.5,
        -0.1,
        "contracted",
        va="center",
        ha="left",
        fontsize=11,
        rotation=0,
        color="black",
        transform=cbar.ax.transAxes,
    )
    cbar.ax.text(
        1.5,
        1.1,
        "elongated",
        va="center",
        ha="left",
        fontsize=11,
        rotation=0,
        color="black",
        transform=cbar.ax.transAxes,
    )

    from matplotlib.widgets import Slider, Button

    ax_slider = fig.add_subplot(gs[1, 0])
    plt.subplots_adjust(bottom=0.15)
    slider = Slider(
        ax_slider,
        "Step",
        1,
        nt,
        valinit=1,
        valstep=1,
        valfmt="%d",
    )
    ax_slider.set_facecolor("white")

    # Add play button to the right of the slider
    button_width = 0.08
    button_height = 0.04
    button_left = ax_slider.get_position().x1 + 0.01
    button_bottom = ax_slider.get_position().y0 + 0.01
    ax_play = fig.add_axes([button_left, button_bottom, button_width, button_height])
    play_button = Button(ax_play, "Play")

    import time

    play_state = {"playing": False}

    def play_animation(event):
        if play_state["playing"]:
            play_state["playing"] = False
            play_button.label.set_text("Play")
            return
        play_state["playing"] = True
        play_button.label.set_text("Pause")
        idx = int(slider.val) - 1
        while play_state["playing"] and idx < nt - 1 and plt.fignum_exists(fig.number):
            idx += 1
            slider.set_val(idx + 1)
            plt.pause(t_per_step)
        play_state["playing"] = False
        play_button.label.set_text("Play")

    play_button.on_clicked(play_animation)

    def update(val):
        idx = int(slider.val) - 1
        draw_state(idx)
        fig.canvas.draw_idle()

    slider.on_changed(update)
    plt.show()
