from __future__ import annotations

import copy
import importlib.util
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import numpy as np
import yaml
from matplotlib.lines import Line2D

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

PROJECT_DIR = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_DIR / "src"
DATA_DIR = PROJECT_DIR / "data" / "TUDELFT_V3_KITE"
RESULTS_DIR = PROJECT_DIR / "results" / "TUDELFT_V3_KITE"

if SRC_DIR.exists() and str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

try:
    from kitesim import read_struc_geometry_yaml_level_2
except ImportError:
    read_struc_geometry_yaml_level_2 = None

LEVEL_1_YAML_PATH = DATA_DIR / "struc_geometry_level_1_manual.yaml"
LEVEL_2_YAML_PATH = DATA_DIR / "struc_geometry_level_2_manual.yaml"

OUTPUT_PDF = RESULTS_DIR / "compare_raw_struc_geometry_yaml_level_1_and_2.pdf"
OUTPUT_SVG = RESULTS_DIR / "compare_raw_struc_geometry_yaml_level_1_and_2.svg"

REFERENCE_STYLE_PATH = Path(
    "/home/jellepoland/ownCloud/phd/code/dissertation_plot_styling.py"
)
USE_REFERENCE_PLOT_STYLE = True
USE_GENERATED_LEVEL_2_GEOMETRY = True

FIGSIZE = (9.0, 4.6)
VIEW_ELEVATION = 20
VIEW_AZIMUTH = -140
MARGIN_FRACTION = 0.06
TRANSPARENT_BACKGROUND = True

LINEWIDTH_WING = 0.7
LINEWIDTH_LE = 4.0 * LINEWIDTH_WING
LINEWIDTH_STRUT = 3.0 * LINEWIDTH_WING
LINEWIDTH_BRIDLE = 0.75
LINEWIDTH_PULLEY = LINEWIDTH_BRIDLE

NODE_SIZE = 0
PULLEY_NODE_SIZE = 0
TITLE_SIZE = 10
LEGEND_FONT_SIZE = 10

LEVEL_1_TITLE = "Level 1 manual"
LEVEL_2_TITLE = "Level 2 manual expanded"

COLORS = {
    "leading edge": "black",
    "struts": "black",
    "trailing edge": "#0072B2",
    "diagonals": "#0072B2",
    "wing": "#0072B2",
    "bridles": "0.45",
    "pulley lines": "red",
    "pulley nodes": "#0072B2",
    "nodes": "red",
}

PLOT_ORDER = [
    "trailing edge",
    "diagonals",
    "wing",
    "bridles",
    "pulley lines",
    "struts",
    "leading edge",
]

SAVE_FIGURE = True
SHOW = True


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


def apply_reference_plot_style() -> None:
    if not USE_REFERENCE_PLOT_STYLE or not REFERENCE_STYLE_PATH.exists():
        return

    spec = importlib.util.spec_from_file_location(
        "dissertation_plot_styling", REFERENCE_STYLE_PATH
    )
    if spec is None or spec.loader is None:
        return
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.set_plot_style()


def load_yaml(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"Structural geometry file not found: {path}")
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def coords_from_raw_yaml(cfg: dict) -> np.ndarray:
    rows = []
    rows.extend(cfg.get("wing_particles", {}).get("data", []))
    rows.extend(cfg.get("bridle_particles", {}).get("data", []))

    max_id = max([0, *(int(row[0]) for row in rows)])
    coords = np.full((max_id + 1, 3), np.nan, dtype=float)
    coords[0] = np.asarray(cfg.get("bridle_point_node", [0.0, 0.0, 0.0]), dtype=float)

    for row in rows:
        coords[int(row[0])] = [float(row[1]), float(row[2]), float(row[3])]

    return coords


# ---------------------------------------------------------------------------
# Connectivity
# ---------------------------------------------------------------------------


def add_segment(
    connectivity: dict[str, list[tuple[int, int]]],
    group: str,
    id1: int,
    id2: int,
) -> None:
    connectivity.setdefault(group, []).append((int(id1), int(id2)))


def add_sequence_segments(
    connectivity: dict[str, list[tuple[int, int]]],
    group: str,
    node_ids: list[int],
) -> None:
    for id1, id2 in zip(node_ids[:-1], node_ids[1:]):
        add_segment(connectivity, group, id1, id2)


def add_bridle_connection(
    connectivity: dict[str, list[tuple[int, int]]],
    row: list,
) -> None:
    node_ids = [int(value) for value in row[1:]]
    if len(node_ids) == 3:
        add_segment(connectivity, "pulley lines", node_ids[0], node_ids[1])
        add_segment(connectivity, "pulley lines", node_ids[1], node_ids[2])
        connectivity.setdefault("__pulley_nodes__", []).append(
            (node_ids[1], node_ids[1])
        )
    elif len(node_ids) >= 2:
        add_segment(connectivity, "bridles", node_ids[0], node_ids[1])


def group_from_wing_connection_name(name: str) -> str:
    name = name.lower()
    if name.startswith("le"):
        return "leading edge"
    if name.startswith("strut"):
        return "struts"
    if name.startswith("te"):
        return "trailing edge"
    if name.startswith("dia"):
        return "diagonals"
    return "wing"


def level_1_connectivity_by_group(cfg: dict) -> dict[str, list[tuple[int, int]]]:
    connectivity: dict[str, list[tuple[int, int]]] = {}

    for name, ci, cj in cfg.get("wing_connections", {}).get("data", []):
        add_segment(connectivity, group_from_wing_connection_name(str(name)), ci, cj)

    for row in cfg.get("bridle_connections", {}).get("data", []):
        add_bridle_connection(connectivity, row)

    return connectivity


def level_2_connectivity_by_group(cfg: dict) -> dict[str, list[tuple[int, int]]]:
    connectivity: dict[str, list[tuple[int, int]]] = {}

    for row in cfg.get("strut_tubes", {}).get("data", []):
        name = str(row[0]).lower()
        group = "struts" if name.startswith("strut") else "leading edge"
        node_ids = row[6] if len(row) > 6 and isinstance(row[6], list) else row[1:3]
        for id1, id2 in zip(node_ids[:-1], node_ids[1:]):
            add_segment(connectivity, group, id1, id2)

    for row in cfg.get("leading_edge_tubes", {}).get("data", []):
        add_segment(connectivity, "leading edge", row[1], row[2])

    wing_rows = cfg.get("wing_particles", {}).get("data", [])
    wing_pairs = [
        (int(le_row[0]), int(te_row[0]))
        for le_row, te_row in zip(wing_rows[0::2], wing_rows[1::2])
    ]
    for (_, te1), (_, te2) in zip(wing_pairs[:-1], wing_pairs[1:]):
        add_segment(connectivity, "trailing edge", te1, te2)
    for (le1, te1), (le2, te2) in zip(wing_pairs[:-1], wing_pairs[1:]):
        add_segment(connectivity, "diagonals", le1, te2)
        add_segment(connectivity, "diagonals", te1, le2)

    for row in cfg.get("bridle_connections", {}).get("data", []):
        add_bridle_connection(connectivity, row)

    return connectivity


def add_generated_canopy_mesh(
    connectivity: dict[str, list[tuple[int, int]]],
    canopy_sections: list[list[int]],
    strut_sections: list[list[int]],
) -> None:
    for canopy_section in canopy_sections:
        add_sequence_segments(connectivity, "wing", canopy_section)

    all_sections = canopy_sections + strut_sections
    all_sections.sort(key=lambda section: section[0])

    for section1, section2 in zip(all_sections[:-1], all_sections[1:]):
        for idx in range(1, min(len(section1), len(section2))):
            group = (
                "trailing edge"
                if idx == len(section1) - 1 and idx == len(section2) - 1
                else "wing"
            )
            add_segment(connectivity, group, section1[idx], section2[idx])
            add_segment(connectivity, "diagonals", section1[idx], section2[idx - 1])
            add_segment(connectivity, "diagonals", section1[idx - 1], section2[idx])


def generated_level_2_connectivity_by_group(
    cfg: dict,
    canopy_sections: list[list[int]],
    strut_sections: list[list[int]],
    bridle_connectivity_arr: np.ndarray,
    n_wing_connections: int,
    pulley_node_ids: list[int],
    pulley_line_indices: list[int],
) -> dict[str, list[tuple[int, int]]]:
    connectivity: dict[str, list[tuple[int, int]]] = {}

    for row, section in zip(cfg.get("strut_tubes", {}).get("data", []), strut_sections):
        name = str(row[0]).lower()
        group = "struts" if name.startswith("strut") else "leading edge"
        add_sequence_segments(connectivity, group, section)

    for row in cfg.get("leading_edge_tubes", {}).get("data", []):
        add_segment(connectivity, "leading edge", row[1], row[2])

    add_generated_canopy_mesh(connectivity, canopy_sections, strut_sections)

    pulley_line_indices = set(int(idx) for idx in pulley_line_indices)
    for local_idx, (id1, id2) in enumerate(bridle_connectivity_arr):
        group = (
            "pulley lines"
            if n_wing_connections + local_idx in pulley_line_indices
            else "bridles"
        )
        add_segment(connectivity, group, id1, id2)

    connectivity["__pulley_nodes__"] = [
        (int(node_id), int(node_id)) for node_id in pulley_node_ids
    ]
    return connectivity


def level_2_plot_data(cfg: dict) -> tuple[np.ndarray, dict[str, list[tuple[int, int]]]]:
    if not USE_GENERATED_LEVEL_2_GEOMETRY:
        return coords_from_raw_yaml(cfg), level_2_connectivity_by_group(cfg)

    if read_struc_geometry_yaml_level_2 is None:
        raise ImportError(
            "Could not import kitesim.read_struc_geometry_yaml_level_2. "
            "Set USE_GENERATED_LEVEL_2_GEOMETRY = False to plot only raw YAML nodes."
        )

    generated_cfg = copy.deepcopy(cfg)
    (
        struc_nodes,
        _m_arr,
        _struc_node_le_indices,
        _struc_node_te_indices,
        _power_tape_index,
        _steering_tape_indices,
        pulley_node_ids,
        canopy_sections,
        strut_sections,
        _simplified_bridle_points,
        kite_connectivity_arr,
        bridle_connectivity_arr,
        _bridle_diameter_arr,
        _l0_arr,
        _k_arr,
        _c_arr,
        _linktype_arr,
        pulley_line_indices,
        _pulley_line_to_other_node_pair_dict,
    ) = read_struc_geometry_yaml_level_2.main(generated_cfg)

    n_wing_connections = len(kite_connectivity_arr) - len(bridle_connectivity_arr)
    connectivity = generated_level_2_connectivity_by_group(
        generated_cfg,
        canopy_sections,
        strut_sections,
        bridle_connectivity_arr,
        n_wing_connections,
        pulley_node_ids,
        pulley_line_indices,
    )

    return np.asarray(struc_nodes, dtype=float), connectivity


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------


def valid_points(coords: np.ndarray) -> np.ndarray:
    return coords[np.isfinite(coords).all(axis=1)]


def set_equal_3d_limits(ax, coords: np.ndarray) -> None:
    pts = valid_points(coords)
    if len(pts) == 0:
        return

    mins = pts.min(axis=0)
    maxs = pts.max(axis=0)
    center = 0.5 * (mins + maxs)
    half_range = 0.5 * float(np.max(maxs - mins))
    half_range *= 1.0 + MARGIN_FRACTION
    if half_range == 0:
        half_range = 1.0

    ax.set_xlim(center[0] - half_range, center[0] + half_range)
    ax.set_ylim(center[1] - half_range, center[1] + half_range)
    ax.set_zlim(center[2] - half_range, center[2] + half_range)
    ax.set_box_aspect((1.0, 1.0, 1.0))


def segment_is_valid(coords: np.ndarray, id1: int, id2: int) -> bool:
    if id1 < 0 or id2 < 0 or id1 >= len(coords) or id2 >= len(coords):
        return False
    return bool(np.isfinite(coords[[id1, id2]]).all())


def linewidth_for_group(group: str) -> float:
    if group == "bridles":
        return LINEWIDTH_BRIDLE
    if group == "pulley lines":
        return LINEWIDTH_PULLEY
    if group == "leading edge":
        return LINEWIDTH_LE
    if group == "struts":
        return LINEWIDTH_STRUT
    return LINEWIDTH_WING


def plot_connectivity(
    ax,
    coords: np.ndarray,
    connectivity: dict[str, list[tuple[int, int]]],
) -> None:
    plot_order = [
        *PLOT_ORDER,
        *(
            group
            for group in connectivity
            if group not in PLOT_ORDER and not group.startswith("__")
        ),
    ]

    for group in plot_order:
        for id1, id2 in connectivity.get(group, []):
            if not segment_is_valid(coords, id1, id2):
                continue
            p1 = coords[id1]
            p2 = coords[id2]
            ax.plot(
                [p1[0], p2[0]],
                [p1[1], p2[1]],
                [p1[2], p2[2]],
                color=COLORS.get(group, COLORS["wing"]),
                lw=linewidth_for_group(group),
                alpha=1.0,
                zorder=2,
            )


def pulley_node_indices(connectivity: dict[str, list[tuple[int, int]]]) -> list[int]:
    return sorted({idx for idx, _ in connectivity.get("__pulley_nodes__", [])})


def plot_nodes(
    ax,
    coords: np.ndarray,
    pulley_nodes: list[int],
) -> None:
    pts = valid_points(coords)
    if len(pts):
        ax.scatter(
            pts[:, 0],
            pts[:, 1],
            pts[:, 2],
            s=NODE_SIZE,
            color=COLORS["nodes"],
            alpha=1.0,
            depthshade=False,
            zorder=3,
        )

    pulley = [
        idx
        for idx in pulley_nodes
        if 0 <= idx < len(coords) and np.isfinite(coords[idx]).all()
    ]
    if not pulley:
        return

    pulley_pts = coords[pulley]
    ax.scatter(
        pulley_pts[:, 0],
        pulley_pts[:, 1],
        pulley_pts[:, 2],
        s=PULLEY_NODE_SIZE,
        color=COLORS["pulley nodes"],
        depthshade=False,
        zorder=12,
    )


def clean_3d_axis(ax) -> None:
    ax.computed_zorder = False
    ax.grid(False)
    ax.set_axis_off()
    ax.set_facecolor("none")
    ax.view_init(elev=VIEW_ELEVATION, azim=VIEW_AZIMUTH)


def draw_model(
    ax,
    title: str,
    coords: np.ndarray,
    connectivity: dict[str, list[tuple[int, int]]],
) -> None:
    plot_connectivity(ax, coords, connectivity)
    plot_nodes(ax, coords, pulley_node_indices(connectivity))
    set_equal_3d_limits(ax, coords)
    clean_3d_axis(ax)
    ax.set_title(title, fontsize=TITLE_SIZE)


def legend_handles() -> list[Line2D]:
    return [
        Line2D([0], [0], color=COLORS["leading edge"], lw=LINEWIDTH_LE, label="LE"),
        Line2D([0], [0], color=COLORS["struts"], lw=LINEWIDTH_STRUT, label="Struts"),
        Line2D(
            [0], [0], color=COLORS["wing"], lw=LINEWIDTH_WING, label="Wing elements"
        ),
        Line2D(
            [0],
            [0],
            color=COLORS["bridles"],
            lw=LINEWIDTH_BRIDLE,
            label="Bridle elements",
        ),
        Line2D(
            [0],
            [0],
            color=COLORS["pulley lines"],
            lw=LINEWIDTH_PULLEY,
            label="Pulley lines",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=COLORS["nodes"],
            markeredgecolor=COLORS["nodes"],
            markersize=4,
            label="Particles",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=COLORS["pulley nodes"],
            markeredgecolor=COLORS["pulley nodes"],
            markersize=5,
            label="Pulley nodes",
        ),
    ]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main(show: bool = SHOW, save_figure: bool = SAVE_FIGURE) -> None:
    apply_reference_plot_style()

    level_1_cfg = load_yaml(LEVEL_1_YAML_PATH)
    level_2_cfg = load_yaml(LEVEL_2_YAML_PATH)

    level_1_coords = coords_from_raw_yaml(level_1_cfg)
    level_2_coords, level_2_connectivity = level_2_plot_data(level_2_cfg)

    fig = plt.figure(figsize=FIGSIZE)
    if TRANSPARENT_BACKGROUND:
        fig.patch.set_alpha(0.0)

    axes = [
        fig.add_subplot(1, 2, 1, projection="3d"),
        fig.add_subplot(1, 2, 2, projection="3d"),
    ]

    draw_model(
        axes[0],
        LEVEL_1_TITLE,
        level_1_coords,
        level_1_connectivity_by_group(level_1_cfg),
    )
    draw_model(
        axes[1],
        LEVEL_2_TITLE,
        level_2_coords,
        level_2_connectivity,
    )

    fig.legend(
        handles=legend_handles(),
        loc="lower center",
        ncol=5,
        frameon=False,
        fontsize=LEGEND_FONT_SIZE,
    )
    plt.tight_layout(rect=[0, 0.08, 1, 1])

    if save_figure:
        RESULTS_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(
            OUTPUT_PDF,
            bbox_inches="tight",
            transparent=TRANSPARENT_BACKGROUND,
        )
        fig.savefig(
            OUTPUT_SVG,
            bbox_inches="tight",
            transparent=TRANSPARENT_BACKGROUND,
        )
        print(f"[compare_raw_struc_geometry_yaml_level_1_and_2] wrote {OUTPUT_PDF}")
        print(f"[compare_raw_struc_geometry_yaml_level_1_and_2] wrote {OUTPUT_SVG}")

    if show:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    args = {arg.lower() for arg in sys.argv[1:]}
    show = SHOW
    save_figure = SAVE_FIGURE
    if "--show" in args:
        show = True
    if "--no-show" in args:
        show = False
    if "--save" in args:
        save_figure = True
    if "--no-save" in args:
        save_figure = False

    main(show=show, save_figure=save_figure)
