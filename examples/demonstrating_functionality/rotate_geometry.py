from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

# Allow running this script directly from the repository root without
# requiring an editable install or PYTHONPATH setup.
PROJECT_DIR = Path(__file__).resolve().parents[2]
SRC_DIR = PROJECT_DIR / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from kitesim import read_struc_geometry_yaml_level_1
from kitesim.utils import load_yaml, rotate_geometry


def _resolve_level_1_geometry_path(project_dir):
    """Resolve TUDELFT_V3_KITE level-1 geometry file path."""
    candidate_paths = [
        project_dir
        / "data"
        / "TUDELFT_V3_KITE"
        / "struc_geometry_level_1_manual_JULIA.yaml",
        project_dir / "data" / "TUDELFT_V3_KITE" / "struc_geometry_level_1_manual.yaml",
    ]
    for path in candidate_paths:
        if path.exists():
            return path
    raise FileNotFoundError(
        "Could not find a level-1 geometry YAML. Checked:\n"
        + "\n".join(str(p) for p in candidate_paths)
    )


def _set_axes_equal_3d(ax, points):
    """Set equal scale on all three 3D axes based on point cloud bounds."""
    mins = points.min(axis=0)
    maxs = points.max(axis=0)
    center = 0.5 * (mins + maxs)
    span = np.max(maxs - mins)
    half = 0.5 * span if span > 0 else 1.0

    ax.set_xlim(center[0] - half, center[0] + half)
    ax.set_ylim(center[1] - half, center[1] + half)
    ax.set_zlim(center[2] - half, center[2] + half)
    ax.set_proj_type("ortho")
    if hasattr(ax, "set_box_aspect"):
        ax.set_box_aspect((1.0, 1.0, 1.0))


def _plot_structure(
    ax,
    nodes,
    connectivity,
    color_line,
    color_node,
    alpha=1.0,
    plot_pivot=False,
):
    """Plot structural nodes and connectivity in a 3D axis."""
    for ci, cj in connectivity:
        i = int(ci)
        j = int(cj)
        p1 = nodes[i]
        p2 = nodes[j]
        ax.plot(
            [p1[0], p2[0]],
            [p1[1], p2[1]],
            [p1[2], p2[2]],
            color=color_line,
            linewidth=1.0,
            alpha=alpha,
        )

    ax.scatter(
        nodes[:, 0],
        nodes[:, 1],
        nodes[:, 2],
        s=10,
        color=color_node,
        alpha=min(1.0, alpha + 0.1),
    )
    if plot_pivot:
        ax.scatter(
            [nodes[0, 0]],
            [nodes[0, 1]],
            [nodes[0, 2]],
            s=35,
            color="tab:red",
        )

    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_zlabel("z [m]")


def main():
    project_dir = PROJECT_DIR
    struc_geometry_path = _resolve_level_1_geometry_path(project_dir)
    struc_geometry = load_yaml(struc_geometry_path)

    geometry_data = read_struc_geometry_yaml_level_1.main(struc_geometry)
    struc_nodes = np.array(geometry_data[0], dtype=float)
    kite_connectivity_arr = np.array(geometry_data[7], dtype=int)

    step_definitions = [
        ("Step 1: +30 deg about x", [30.0, 0.0, 0.0], "x"),
        ("Step 2: +30 deg about y", [0.0, 30.0, 0.0], "y"),
        ("Step 3: +30 deg about z", [0.0, 0.0, 30.0], "z"),
    ]

    individual_pairs = []
    for _, angles_deg, axis_name in step_definitions:
        nodes_after = rotate_geometry(struc_nodes, angle_deg=angles_deg)
        individual_pairs.append((axis_name, struc_nodes.copy(), nodes_after))

    fig, axes = plt.subplots(1, 3, figsize=(18, 6), subplot_kw={"projection": "3d"})
    all_points = np.vstack([struc_nodes] + [pair[2] for pair in individual_pairs])

    for col, (axis_name, nodes_before, nodes_after) in enumerate(individual_pairs):
        ax = axes[col]
        _plot_structure(
            ax,
            nodes_before,
            kite_connectivity_arr,
            color_line="0.6",
            color_node="0.35",
            alpha=0.35,
        )
        _plot_structure(
            ax,
            nodes_after,
            kite_connectivity_arr,
            color_line="tab:blue",
            color_node="black",
            alpha=0.95,
            plot_pivot=True,
        )
        ax.set_title(f"Individual: +30 deg about {axis_name}")
        ax.plot([], [], [], color="0.6", alpha=0.6, linewidth=2, label="initial")
        ax.plot([], [], [], color="tab:blue", linewidth=2, label="rotated")
        ax.plot([], [], [], color="tab:red", linewidth=0, marker="o", label="pivot (node 0)")
        _set_axes_equal_3d(ax, all_points)
        ax.view_init(elev=24, azim=-58)
        ax.text2D(0.03, 0.95, f"axis: {axis_name}", transform=ax.transAxes, fontsize=9)
        ax.legend(loc="upper right", fontsize=8)

    fig.suptitle(
        f"TUDELFT_V3_KITE level-1 geometry: individual rotations ({struc_geometry_path.name})",
        fontsize=12,
    )
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
