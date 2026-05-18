from pathlib import Path
import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import numpy as np
import yaml


SCRIPT_DIR = Path(__file__).resolve().parent
KITE_DIR = SCRIPT_DIR.parent

AERO_GEOMETRY_PATH = (
    KITE_DIR
    / "2D_airfoils_polars_plots_BEST"
    / "aero_geometry_CFD_CAD_derived.yaml"
)
AIRFOIL_PROFILE_DIR = (
    KITE_DIR / "2D_airfoils_polars_plots_BEST" / "2D_airfoils_sliced_from_CFD_CAD"
)
PSM_GEOMETRY_PATH = (
    SCRIPT_DIR
    / "results"
    / "struc_geometry_PSM_reduced_stretched_mapped_to_FEM_stretched.yaml"
)
SURFPLAN_AERO_GEOMETRY_PATH = Path(
    "/home/jellepoland/ownCloud/phd/code/SurfplanAdapter/processed_data/"
    "TUDELFT_V3_KITE/aero_geometry.yaml"
)


def load_yaml(path):
    with Path(path).open("r") as stream:
        return yaml.safe_load(stream)


def make_header_map(section):
    return {name: idx for idx, name in enumerate(section["headers"])}


def normalize(vector):
    norm = np.linalg.norm(vector)
    if norm == 0.0:
        raise ValueError("Cannot normalize a zero-length vector.")
    return vector / norm


def get_inferred_vup(chord_unit):
    up_hint = np.array([0.0, 0.0, 1.0])
    inferred_vup = up_hint - np.dot(up_hint, chord_unit) * chord_unit
    if np.linalg.norm(inferred_vup) < 1e-12:
        up_hint = np.array([0.0, 1.0, 0.0])
        inferred_vup = up_hint - np.dot(up_hint, chord_unit) * chord_unit
    return inferred_vup


def get_section_frame(le_point, te_point, vup_vector=None):
    chord_vector = te_point - le_point
    chord_length = np.linalg.norm(chord_vector)
    if chord_length == 0.0:
        raise ValueError("Wing section has coincident LE and TE points.")

    x_local = chord_vector / chord_length

    if vup_vector is None:
        vup_vector = get_inferred_vup(x_local)

    # Match the SurfplanAdapter plotting convention: VUP is used directly as
    # local airfoil-up. The spanwise normal follows from chord x VUP.
    y_local = normalize(vup_vector)
    z_local = normalize(np.cross(x_local, y_local))
    return chord_length, x_local, y_local, z_local


def load_airfoil_profile(profile_path):
    coords = []
    with Path(profile_path).open("r") as stream:
        for line in stream:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.replace(",", " ").split()
            if len(parts) < 2:
                continue

            try:
                coords.append([float(parts[0]), float(parts[1])])
            except ValueError:
                continue

    if not coords:
        raise ValueError(f"No airfoil coordinates found in {profile_path}")
    return np.array(coords, dtype=float)


def transform_profile_to_world(profile_coords, le_point, te_point, vup_vector=None):
    chord_length, x_local, y_local, z_local = get_section_frame(
        le_point, te_point, vup_vector
    )

    profile_x = profile_coords[:, 0] * chord_length
    profile_y = profile_coords[:, 1] * chord_length
    profile_z = np.zeros_like(profile_x)

    return (
        le_point
        + profile_x[:, None] * x_local
        + profile_y[:, None] * y_local
        + profile_z[:, None] * z_local
    )


def load_surfplan_vup_sections(surfplan_aero_geometry_path):
    surfplan_aero_geometry_path = Path(surfplan_aero_geometry_path)
    if not surfplan_aero_geometry_path.exists():
        return []

    config = load_yaml(surfplan_aero_geometry_path)
    wing_sections = config.get("wing_sections", {})
    if not wing_sections:
        return []

    headers = make_header_map(wing_sections)
    required_headers = {
        "LE_y",
        "TE_y",
        "VUP_x",
        "VUP_y",
        "VUP_z",
    }
    if not required_headers.issubset(headers):
        return []

    vup_sections = []
    for row in wing_sections["data"]:
        le_y = float(row[headers["LE_y"]])
        te_y = float(row[headers["TE_y"]])
        y_mid = 0.5 * (le_y + te_y)
        vup = np.array(
            [
                float(row[headers["VUP_x"]]),
                float(row[headers["VUP_y"]]),
                float(row[headers["VUP_z"]]),
            ]
        )
        if np.linalg.norm(vup) == 0.0:
            continue
        vup_sections.append({"y_mid": y_mid, "abs_y_mid": abs(y_mid), "vup": vup})

    return vup_sections


def map_surfplan_vup_to_section(le_point, te_point, surfplan_vup_sections):
    if not surfplan_vup_sections:
        return None

    target_y = 0.5 * (le_point[1] + te_point[1])
    abs_target_y = abs(target_y)

    if abs_target_y < 1e-9:
        positive_sections = [
            section for section in surfplan_vup_sections if section["y_mid"] > 0.0
        ]
        negative_sections = [
            section for section in surfplan_vup_sections if section["y_mid"] < 0.0
        ]
        if positive_sections and negative_sections:
            positive_match = min(
                positive_sections,
                key=lambda section: abs(section["abs_y_mid"] - abs_target_y),
            )
            negative_match = min(
                negative_sections,
                key=lambda section: abs(section["abs_y_mid"] - abs_target_y),
            )
            center_vup = 0.5 * (positive_match["vup"] + negative_match["vup"])
            center_vup[1] = 0.0
            return center_vup

    target_side = np.sign(target_y)
    same_side_sections = [
        section
        for section in surfplan_vup_sections
        if np.sign(section["y_mid"]) == target_side
    ]
    candidates = same_side_sections if same_side_sections else surfplan_vup_sections
    return min(candidates, key=lambda section: abs(section["y_mid"] - target_y))["vup"]


def load_psm_wing_nodes_and_connections(psm_geometry_path):
    config = load_yaml(psm_geometry_path)

    particle_headers = make_header_map(config["wing_particles"])
    connection_headers = make_header_map(config["wing_connections"])

    wing_nodes = {}
    for row in config["wing_particles"]["data"]:
        node_id = int(row[particle_headers["id"]])
        wing_nodes[node_id] = np.array(
            [
                float(row[particle_headers["x"]]),
                float(row[particle_headers["y"]]),
                float(row[particle_headers["z"]]),
            ]
        )

    wing_connections = []
    for row in config["wing_connections"]["data"]:
        wing_connections.append(
            (
                str(row[connection_headers["name"]]),
                int(row[connection_headers["ci"]]),
                int(row[connection_headers["cj"]]),
            )
        )

    return wing_nodes, wing_connections


def set_equal_axes(ax, points):
    points = np.asarray(points, dtype=float)
    mins = points.min(axis=0)
    maxs = points.max(axis=0)
    center = 0.5 * (mins + maxs)
    max_range = np.max(maxs - mins)
    half_range = 0.55 * max_range

    ax.set_xlim(center[0] - half_range, center[0] + half_range)
    ax.set_ylim(center[1] - half_range, center[1] + half_range)
    ax.set_zlim(center[2] - half_range, center[2] + half_range)
    ax.set_box_aspect((1, 1, 1))


def plot_aero_geometry_and_psm_wing_nodes(
    aero_geometry_path=AERO_GEOMETRY_PATH,
    profile_dir=AIRFOIL_PROFILE_DIR,
    psm_geometry_path=PSM_GEOMETRY_PATH,
    surfplan_aero_geometry_path=SURFPLAN_AERO_GEOMETRY_PATH,
):
    aero_config = load_yaml(aero_geometry_path)
    wing_sections = aero_config["wing_sections"]
    section_headers = make_header_map(wing_sections)
    surfplan_vup_sections = load_surfplan_vup_sections(surfplan_aero_geometry_path)
    if surfplan_vup_sections:
        print(f"Using mapped SurfplanAdapter VUPs from {surfplan_aero_geometry_path}")
    else:
        print("No SurfplanAdapter VUPs found; inferring section-up from global z.")

    fig = plt.figure(figsize=(15, 10))
    ax = fig.add_subplot(111, projection="3d")

    all_points = []
    plotted_airfoil_label = False
    plotted_le_label = False
    plotted_te_label = False
    plotted_chord_label = False

    for section in wing_sections["data"]:
        airfoil_id = int(section[section_headers["airfoil_id"]])
        le_point = np.array(
            [
                float(section[section_headers["LE_x"]]),
                float(section[section_headers["LE_y"]]),
                float(section[section_headers["LE_z"]]),
            ]
        )
        te_point = np.array(
            [
                float(section[section_headers["TE_x"]]),
                float(section[section_headers["TE_y"]]),
                float(section[section_headers["TE_z"]]),
            ]
        )

        profile_path = Path(profile_dir) / f"{airfoil_id}.dat"
        profile_coords = load_airfoil_profile(profile_path)
        vup_vector = None
        if {"VUP_x", "VUP_y", "VUP_z"}.issubset(section_headers):
            vup_vector = np.array(
                [
                    float(section[section_headers["VUP_x"]]),
                    float(section[section_headers["VUP_y"]]),
                    float(section[section_headers["VUP_z"]]),
                ]
            )
        elif surfplan_vup_sections:
            vup_vector = map_surfplan_vup_to_section(
                le_point, te_point, surfplan_vup_sections
            )

        world_profile = transform_profile_to_world(
            profile_coords, le_point, te_point, vup_vector
        )

        ax.plot(
            world_profile[:, 0],
            world_profile[:, 1],
            world_profile[:, 2],
            color="black",
            linewidth=1.0,
            alpha=0.7,
            label="CFD/CAD airfoil sections" if not plotted_airfoil_label else None,
        )
        plotted_airfoil_label = True

        ax.scatter(
            [le_point[0]],
            [le_point[1]],
            [le_point[2]],
            c="tab:blue",
            s=12,
            alpha=0.85,
            label="Aero leading edge" if not plotted_le_label else None,
        )
        plotted_le_label = True

        ax.scatter(
            [te_point[0]],
            [te_point[1]],
            [te_point[2]],
            c="tab:red",
            s=12,
            alpha=0.85,
            label="Aero trailing edge" if not plotted_te_label else None,
        )
        plotted_te_label = True

        ax.plot(
            [le_point[0], te_point[0]],
            [le_point[1], te_point[1]],
            [le_point[2], te_point[2]],
            color="0.25",
            linestyle="--",
            linewidth=0.6,
            alpha=0.45,
            label="Aero chord lines" if not plotted_chord_label else None,
        )
        plotted_chord_label = True

        all_points.extend(world_profile)
        all_points.extend([le_point, te_point])

    wing_nodes, wing_connections = load_psm_wing_nodes_and_connections(psm_geometry_path)
    plotted_connection_label = False
    for _, ci, cj in wing_connections:
        if ci not in wing_nodes or cj not in wing_nodes:
            continue
        p_i = wing_nodes[ci]
        p_j = wing_nodes[cj]
        ax.plot(
            [p_i[0], p_j[0]],
            [p_i[1], p_j[1]],
            [p_i[2], p_j[2]],
            color="tab:orange",
            linewidth=0.8,
            alpha=0.35,
            label="PSM wing connections" if not plotted_connection_label else None,
        )
        plotted_connection_label = True

    psm_points = np.array([wing_nodes[node_id] for node_id in sorted(wing_nodes)])
    ax.scatter(
        psm_points[:, 0],
        psm_points[:, 1],
        psm_points[:, 2],
        c="tab:orange",
        marker="^",
        s=42,
        edgecolors="black",
        linewidths=0.4,
        label="PSM wing nodes",
    )

    for node_id in sorted(wing_nodes):
        point = wing_nodes[node_id]
        ax.text(
            point[0],
            point[1],
            point[2],
            str(node_id),
            color="tab:orange",
            fontsize=8,
        )

    all_points.extend(psm_points)
    set_equal_axes(ax, all_points)

    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.set_zlabel("Z [m]")
    ax.set_title("CFD/CAD Aero Geometry and Mapped PSM Wing Nodes")
    ax.view_init(elev=20, azim=-120)
    ax.grid(True, alpha=0.3)
    ax.legend()
    plt.tight_layout()
    plt.show()


def main():
    plot_aero_geometry_and_psm_wing_nodes()


if __name__ == "__main__":
    main()
