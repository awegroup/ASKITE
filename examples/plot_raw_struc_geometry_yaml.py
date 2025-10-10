import numpy as np
import matplotlib.pyplot as plt
import yaml
from pathlib import Path


def plot_struct_geometry_all_in_surfplan_yaml(yaml_path, show_plot=True):
    yaml_path = Path(yaml_path)
    if not yaml_path.exists():
        raise FileNotFoundError(f"Structural geometry file not found: {yaml_path}")

    with open(yaml_path, "r") as f:
        config = yaml.safe_load(f)

    wing_particles_cfg = config.get("wing_particles", {})
    wing_particles_data = wing_particles_cfg.get("data", [])
    strut_tubes_cfg = config.get("strut_tubes", {})
    strut_tubes_data = strut_tubes_cfg.get("data", [])
    bridle_particles_cfg = config.get("bridle_particles", {})
    bridle_particles_data = bridle_particles_cfg.get("data", [])
    bridle_connections_cfg = config.get("bridle_connections", {})
    bridle_connections_data = bridle_connections_cfg.get("data", [])
    bridle_point_node = config.get("bridle_point_node", [0, 0, 0])

    if not wing_particles_data and not bridle_particles_data:
        print("No structural data available to plot.")
        return

    # Build coordinate maps
    wing_coords = {
        int(row[0]): np.array([float(row[1]), float(row[2]), float(row[3])])
        for row in wing_particles_data
    }
    bridle_coords = {
        int(row[0]): np.array([float(row[1]), float(row[2]), float(row[3])])
        for row in bridle_particles_data
    }

    # Add bridle_point node (node 0 - KCU attachment point)
    bridle_point_coord = np.array(
        [
            float(bridle_point_node[0]),
            float(bridle_point_node[1]),
            float(bridle_point_node[2]),
        ]
    )
    bridle_coords[0] = bridle_point_coord

    all_coords = dict(wing_coords)
    all_coords.update(bridle_coords)

    strut_node_ids = set()
    for entry in strut_tubes_data:
        if len(entry) >= 3:
            strut_node_ids.add(int(entry[1]))
            strut_node_ids.add(int(entry[2]))

    # Prepare figure
    fig = plt.figure(figsize=(15, 10))
    ax = fig.add_subplot(111, projection="3d")
    ax.set_title(yaml_path.stem)
    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.set_zlabel("Z [m]")
    ax.grid(True, alpha=0.3)
    ax.view_init(elev=20, azim=-120)

    base_marker = 20

    # Plot wing nodes
    wing_points_strut = []
    wing_sizes_strut = []

    wing_points_other = []
    wing_sizes_other = []

    for node_id, coord in wing_coords.items():
        is_strut_node = node_id in strut_node_ids
        is_le_node = node_id % 2 == 1  # LE nodes were added first for each rib
        marker_size = base_marker * (2 if is_le_node else 1)

        if is_strut_node:
            wing_points_strut.append(coord)
            wing_sizes_strut.append(marker_size)
        else:
            wing_points_other.append(coord)
            wing_sizes_other.append(marker_size)

    if wing_points_other:
        wing_points_other = np.array(wing_points_other)
        ax.scatter(
            wing_points_other[:, 0],
            wing_points_other[:, 1],
            wing_points_other[:, 2],
            c="dimgray",
            s=wing_sizes_other,
            alpha=0.8,
            label="Wing Node",
        )
    if wing_points_strut:
        wing_points_strut = np.array(wing_points_strut)
        ax.scatter(
            wing_points_strut[:, 0],
            wing_points_strut[:, 1],
            wing_points_strut[:, 2],
            c="blue",
            s=wing_sizes_strut,
            alpha=0.9,
            label="Strut Node",
        )

    # Annotate wing nodes
    for node_id, coord in wing_coords.items():
        color = "blue" if node_id in strut_node_ids else "dimgray"
        ax.text(
            coord[0],
            coord[1],
            coord[2],
            f"{node_id}",
            fontsize=8,
            color=color,
        )

    # Plot chord lines (LE to TE pairs in order)
    wing_ids_sorted = sorted(wing_coords.keys())
    for i in range(0, len(wing_ids_sorted), 2):
        if i + 1 >= len(wing_ids_sorted):
            continue
        le_id = wing_ids_sorted[i]
        te_id = wing_ids_sorted[i + 1]
        if le_id in wing_coords and te_id in wing_coords:
            le_coord = wing_coords[le_id]
            te_coord = wing_coords[te_id]
            ax.plot(
                [le_coord[0], te_coord[0]],
                [le_coord[1], te_coord[1]],
                [le_coord[2], te_coord[2]],
                c="black",
                linewidth=1.0,
                alpha=0.8,
            )

    # Plot bridle nodes
    if bridle_coords:
        # Separate node 0 (KCU) from other bridle nodes
        bridle_points_regular = []
        for node_id, coord in bridle_coords.items():
            if node_id != 0:
                bridle_points_regular.append(coord)

        # Plot regular bridle nodes
        if bridle_points_regular:
            bridle_points_regular = np.array(bridle_points_regular)
            ax.scatter(
                bridle_points_regular[:, 0],
                bridle_points_regular[:, 1],
                bridle_points_regular[:, 2],
                c="darkorange",
                s=base_marker,
                alpha=0.8,
                label="Bridle Node",
            )

        # Plot node 0 (KCU) with distinctive marker
        if 0 in bridle_coords:
            kcu_coord = bridle_coords[0]
            ax.scatter(
                [kcu_coord[0]],
                [kcu_coord[1]],
                [kcu_coord[2]],
                c="red",
                s=base_marker * 3,
                marker="^",
                alpha=1.0,
                label="KCU (Node 0)",
                edgecolors="black",
                linewidths=2,
            )

        # Annotate all bridle nodes
        for node_id, coord in bridle_coords.items():
            color = "red" if node_id == 0 else "darkorange"
            fontsize = 10 if node_id == 0 else 8
            fontweight = "bold" if node_id == 0 else "normal"
            ax.text(
                coord[0],
                coord[1],
                coord[2],
                f"{node_id}",
                fontsize=fontsize,
                color=color,
                weight=fontweight,
            )

    # Plot bridle connections
    pulley_plotted = False  # Track if we've added the pulley label
    regular_plotted = False  # Track if we've added the regular bridle label

    for conn in bridle_connections_data:
        if len(conn) < 3:
            continue

        # Check if this is a pulley connection (4 elements: name, ci, ck, cj)
        is_pulley = len(conn) >= 4

        if is_pulley:
            # Pulley connection: ci -> ck -> cj (two line segments)
            ci = int(conn[1])
            ck = int(conn[2])  # Pulley node
            cj = int(conn[3])

            if ci in all_coords and ck in all_coords and cj in all_coords:
                p1 = all_coords[ci]
                pk = all_coords[ck]
                p2 = all_coords[cj]

                # Plot first segment (ci -> ck)
                ax.plot(
                    [p1[0], pk[0]],
                    [p1[1], pk[1]],
                    [p1[2], pk[2]],
                    c="purple",
                    linewidth=1.5,
                    alpha=0.8,
                    label="Pulley Connection" if not pulley_plotted else "",
                )
                # Plot second segment (ck -> cj)
                ax.plot(
                    [pk[0], p2[0]],
                    [pk[1], p2[1]],
                    [pk[2], p2[2]],
                    c="purple",
                    linewidth=1.5,
                    alpha=0.8,
                )
                pulley_plotted = True
        else:
            # Regular connection: ci -> cj (one line segment)
            ci = int(conn[1])
            cj = int(conn[2])

            if ci in all_coords and cj in all_coords:
                p1 = all_coords[ci]
                p2 = all_coords[cj]
                ax.plot(
                    [p1[0], p2[0]],
                    [p1[1], p2[1]],
                    [p1[2], p2[2]],
                    c="darkorange",
                    linewidth=1.0,
                    alpha=0.7,
                    label="Bridle Connection" if not regular_plotted else "",
                )
                regular_plotted = True

    # Equal aspect ratio
    if all_coords:
        coord_array = np.array(list(all_coords.values()))
        xyz_min = coord_array.min(axis=0)
        xyz_max = coord_array.max(axis=0)
        centers = (xyz_max + xyz_min) / 2.0
        max_range = (xyz_max - xyz_min).max()
        half_range = max_range * 0.6 if max_range > 0 else 1.0
        ax.set_xlim(centers[0] - half_range, centers[0] + half_range)
        ax.set_ylim(centers[1] - half_range, centers[1] + half_range)
        ax.set_zlim(centers[2] - half_range, centers[2] + half_range)

    ax.legend(loc="upper right")

    if show_plot:
        plt.show()
    else:
        plt.close(fig)


if __name__ == "__main__":
    # Example usage
    PROJECT_DIR = Path(__file__).resolve().parents[1]

    example_yaml_path = (
        Path(PROJECT_DIR)
        / "data"
        / "TUDELFT_V3_KITE"
        / "struc_geometry_level_2_manual.yaml"
    )
    plot_struct_geometry_all_in_surfplan_yaml(example_yaml_path)
