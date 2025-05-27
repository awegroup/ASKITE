import numpy as np
import matplotlib.pyplot as plt


def main(
    struc_nodes,
    kite_connectivity,
    rest_lengths,
    struc_nodes_initial=None,
    f_ext=None,
    title="PSM State",
):
    """
    Plot the current (and optionally initial) structure state in 3D.

    Args:
        struc_nodes (np.ndarray): Current node positions (n_nodes, 3).
        kite_connectivity (array-like): List/array of [i, j, ...] giving spring connections.
        struc_nodes_initial (np.ndarray, optional): Initial node positions (n_nodes, 3).
        f_ext (np.ndarray or None): Optional external forces, shape (n_nodes, 3) or flat.
        title (str): Figure title.

    Returns:
        None. Displays a 3D plot.
    """
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    # Plot initial state if provided
    if struc_nodes_initial is not None:
        ax.scatter(
            *(struc_nodes_initial.T), color="black", marker="o", s=10, label="Initial"
        )
        # Draw initial lines
        cmap = plt.get_cmap("coolwarm")
        norm_values = []
        for idx, (i, j, *rest) in enumerate(kite_connectivity):
            p1, p2 = struc_nodes_initial[i], struc_nodes_initial[j]
            # compute delta rest length
            rest_length = np.linalg.norm(p2 - p1)
            # normalize value for color mapping
            norm_val = (rest_length - rest_lengths[idx]) / rest_lengths[idx]
            norm_values.append(norm_val)
            color = cmap(norm_val)
            ax.plot(
                [p1[0], p2[0]],
                [p1[1], p2[1]],
                [p1[2], p2[2]],
                color=color,
                linewidth=1,
                alpha=0.5,
            )

    # Plot current state
    ax.scatter(*(struc_nodes.T), color="black", marker="o", s=10, label="Current")
    for i, j, *rest in kite_connectivity:
        p1, p2 = struc_nodes[i], struc_nodes[j]
        ax.plot(
            [p1[0], p2[0]],
            [p1[1], p2[1]],
            [p1[2], p2[2]],
            color="black",
            linewidth=1,
        )

    # Optionally plot external forces
    if f_ext is not None:
        arr = np.array(f_ext)
        if arr.ndim == 1:
            arr = arr.reshape(-1, 3)
        for pos, frc in zip(struc_nodes, arr):
            ax.quiver(*pos, *frc, length=1, normalize=True, color="red")

    # Annotate node indices
    for i, pos in enumerate(struc_nodes):
        ax.text(pos[0], pos[1], pos[2], str(i), color="black", fontsize=8)
    if struc_nodes_initial is not None:
        for i, pos in enumerate(struc_nodes_initial):
            ax.text(pos[0], pos[1], pos[2], str(i), color="blue", fontsize=8)

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
    ax.legend()
    # Add colorbar if initial state was plotted
    if struc_nodes_initial is not None and len(norm_values) > 0:
        import matplotlib as mpl

        sm = mpl.cm.ScalarMappable(
            cmap=cmap, norm=plt.Normalize(vmin=min(norm_values), vmax=max(norm_values))
        )
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, shrink=0.7, pad=0.1)
        cbar.set_label("Normalized rest length change")
    plt.show()
    ax.legend()
    plt.show()
