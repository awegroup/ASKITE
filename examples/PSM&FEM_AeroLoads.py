# This example visualizes the initial geometry and aerodynamic loads of the PSM and FE models 
# It produces 3D renderings of the kite structure with aerodynamic load arrows, colored by element category. 
# The aerodynamic loads are scaled to have a common resultant magnitude between the two models for better visual comparison. 
# The resulting figures are saved in PDF, SVG, and PNG formats in the specified output directory.


from pathlib import Path
import copy

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from kitesim.utils import load_sim_output
from kitesim.utils import load_yaml
from kitesim import read_struc_geometry_yaml_level_1, read_struc_geometry_yaml_level_2


PROJECT_DIR = Path(__file__).resolve().parents[1]
RESULTS_ROOT = PROJECT_DIR / "results" / "TUDELFT_V3_KITE"
OUTPUT_DIR = PROJECT_DIR / "results" / "AWEC_2026_initial_load_visuals"

PSM_RUN = RESULTS_ROOT / "2026_06_01_1049h" / "sim_output.h5"
FE_RUN = RESULTS_ROOT / "2026_06_01_1055h" / "sim_output.h5"

VIEW_ELEV = 15
VIEW_AZIM = -120
FORCE_COLOR = "#D62728"
NODE_COLOR = "#1F1F1F"
FORCE_DISPLAY_EXPONENT = 0.65
REFERENCE_NODE_COUNT = 38
REFERENCE_NODE_MARKER_SIZE = 18.0
REFERENCE_FORCE_LINEWIDTH = 0.9

CATEGORY_STYLE = {
    "Inflatable tube spring": {"color": "#0072B2", "linewidth": 1.35, "alpha": 0.84},
    "Power tape": {"color": "#009E73", "linewidth": 1.2, "alpha": 0.94},
    "Pulley spring": {"color": "#E69F00", "linewidth": 0.85, "alpha": 0.82},
    "Non-compressive spring": {"color": "#56B4E9", "linewidth": 0.55, "alpha": 0.46},
    "Inflatable beams": {"color": "#0072B2", "linewidth": 1.7, "alpha": 0.96},
    "Aerodynamic forces": {"color": FORCE_COLOR, "linewidth": 0.01, "alpha": 1.0},
    "Free node point mass": {"color": NODE_COLOR, "linewidth": 0.0, "alpha": 1.0},
    "Pulley node point mass": {"color": "#009E73", "linewidth": 0.0, "alpha": 1.0},
    "KCU / fixed node point mass": {"color": "#B8860B", "linewidth": 0.0, "alpha": 1.0},
    "Free node": {"color": NODE_COLOR, "linewidth": 0.0, "alpha": 1.0},
    "Pulley node": {"color": "#009E73", "linewidth": 0.0, "alpha": 1.0},
    "KCU / fixed node": {"color": "#B8860B", "linewidth": 0.0, "alpha": 1.0},
}

NODE_AND_LOAD_CATEGORIES = {
    "Free node point mass",
    "Pulley node point mass",
    "KCU / fixed node point mass",
    "Free node",
    "Pulley node",
    "KCU / fixed node",
    "Aerodynamic forces",
}

PSM_CATEGORY_ORDER = [
    "Inflatable tube spring",
    "Non-compressive spring",
    "Power tape",
    "Pulley spring",
    "Free node point mass",
    "Pulley node point mass",
    "KCU / fixed node point mass",
    "Aerodynamic forces",
]
FE_CATEGORY_ORDER = [
    "Inflatable beams",
    "Non-compressive spring",
    "Power tape",
    "Pulley spring",
    "Free node",
    "Pulley node",
    "KCU / fixed node",
    "Aerodynamic forces",
]


def get_particle_masses(run_dir, n_nodes):
    struc_geometry = load_yaml(run_dir / "struc_geometry.yaml")
    if n_nodes == 38:
        return read_struc_geometry_yaml_level_1.main(copy.deepcopy(struc_geometry))[1]
    if n_nodes == 268:
        return read_struc_geometry_yaml_level_2.main(copy.deepcopy(struc_geometry))[1]
    raise ValueError(f"Unexpected node count for AWEC visual export: {n_nodes}")


def model_output(geometry, n_nodes):
    if n_nodes == 38:
        return read_struc_geometry_yaml_level_1.main(copy.deepcopy(geometry))
    if n_nodes == 268:
        return read_struc_geometry_yaml_level_2.main(copy.deepcopy(geometry))
    raise ValueError(f"Unexpected node count for AWEC visual export: {n_nodes}")


def expanded_bridle_categories(geometry, bridle_key, pulley_label, spring_label):
    line_lookup = {
        row[0]: dict(zip(geometry[bridle_key]["headers"][1:], row[1:]))
        for row in geometry[bridle_key]["data"]
    }
    categories = []
    for conn_data in geometry["bridle_connections"]["data"]:
        linktype = line_lookup[conn_data[0]]["linktype"]
        if conn_data[0] == "Power Tape":
            category = "Power tape"
        else:
            category = pulley_label if linktype == "pulley" else spring_label
        categories.append(category)
        if len(conn_data[1:]) == 3:
            categories.append(category)
    return categories


def psm_element_categories(geometry):
    categories = []
    for conn_name, *_ in geometry["wing_connections"]["data"]:
        lower_name = conn_name.lower()
        if "le" in lower_name or "strut" in lower_name:
            categories.append("Inflatable tube spring")
        else:
            categories.append("Non-compressive spring")
    categories.extend(
        expanded_bridle_categories(
            geometry,
            "bridle_elements",
            pulley_label="Pulley spring",
            spring_label="Non-compressive spring",
        )
    )
    return categories


def fe_element_categories(geometry, connectivity):
    model_output = read_struc_geometry_yaml_level_2.main(copy.deepcopy(geometry))
    power_tape_index = int(model_output[4])
    linktypes = model_output[16]

    categories = []
    for idx, linktype in enumerate(linktypes):
        if idx == power_tape_index:
            category = "Power tape"
        elif linktype == "pulley":
            category = "Pulley spring"
        elif linktype == "inflatable_beam":
            category = "Inflatable beams"
        else:
            category = "Non-compressive spring"
        categories.append(category)
    return categories


def element_categories(run_dir, connectivity):
    geometry = load_yaml(run_dir / "struc_geometry.yaml")
    if len(connectivity) == 89:
        return psm_element_categories(geometry)
    if len(connectivity) == 881:
        return fe_element_categories(geometry, connectivity)
    raise ValueError(f"Unexpected element count for AWEC visual export: {len(connectivity)}")


def element_widths(geometry, n_nodes, categories):
    if n_nodes != 268:
        return np.ones(len(categories), dtype=float)

    output = read_struc_geometry_yaml_level_2.main(copy.deepcopy(geometry))
    diameters = np.asarray(output[14], dtype=float)
    widths = np.ones(len(categories), dtype=float)
    beam_mask = np.asarray(categories, dtype=object) == "Inflatable beams"
    beam_diameters = diameters[beam_mask]
    if len(beam_diameters) > 0:
        d_min = beam_diameters.min()
        d_max = beam_diameters.max()
        if d_max > d_min:
            widths[beam_mask] = 1.4 + 2.2 * (beam_diameters - d_min) / (d_max - d_min)
        else:
            widths[beam_mask] = 2.2
    return widths


def node_metadata(geometry, n_nodes):
    output = model_output(geometry, n_nodes)
    masses = np.asarray(output[1], dtype=float)
    fixed = {int(i) for i in geometry.get("fixed_point_indices", [])}
    pulley = {int(i) for i in output[6]}
    return masses, fixed, pulley


def load_initial_geometry_with_first_load(h5_path):
    run_dir = h5_path.parent
    meta, tracking = load_sim_output(h5_path)
    positions = np.asarray(tracking["positions"], dtype=float)
    loads = np.asarray(tracking["f_ext"], dtype=float)
    config = load_yaml(run_dir / "config.yaml")
    geometry = load_yaml(run_dir / "struc_geometry.yaml")
    masses, fixed_nodes, pulley_nodes = node_metadata(geometry, positions.shape[1])
    connectivity = np.asarray(meta["kite_connectivity"], dtype=int)
    categories = element_categories(run_dir, connectivity)

    load_step = 0
    load_norms = np.linalg.norm(loads, axis=2).sum(axis=1)
    nonzero_steps = np.flatnonzero(load_norms > 1e-9)
    if len(nonzero_steps) > 0:
        load_step = int(nonzero_steps[0])

    plotted_loads = loads[load_step].copy()
    if config.get("is_with_gravity", False):
        gravity = np.asarray(config["grav_constant"], dtype=float) * masses[:, None]
        plotted_loads = plotted_loads - gravity

    return {
        "nodes": positions[0],
        "n_nodes": positions.shape[1],
        "loads": plotted_loads,
        "connectivity": connectivity,
        "categories": categories,
        "element_width_scale": element_widths(geometry, positions.shape[1], categories),
        "masses": masses,
        "fixed_nodes": fixed_nodes,
        "pulley_nodes": pulley_nodes,
        "load_step": load_step,
    }


def equalize_total_force_magnitude(*datasets):
    resultants = [np.sum(data["loads"], axis=0) for data in datasets]
    magnitudes = [np.linalg.norm(resultant) for resultant in resultants]
    target = float(np.mean([mag for mag in magnitudes if mag > 1e-12]))

    for data, magnitude in zip(datasets, magnitudes):
        data["force_scale"] = 1.0
        if magnitude > 1e-12:
            data["force_scale"] = target / magnitude
            data["loads"] = data["loads"] * data["force_scale"]
        data["total_force_magnitude"] = target
    for data in datasets:
        data["force_arrow_reference"] = np.linalg.norm(data["loads"], axis=1).max()


def set_axes_equal(ax, points, pad=0.08):
    mins = points.min(axis=0)
    maxs = points.max(axis=0)
    center = (mins + maxs) / 2.0
    radius = np.max(maxs - mins) * (0.5 + pad)

    ax.set_xlim(center[0] - radius, center[0] + radius)
    ax.set_ylim(center[1] - radius, center[1] + radius)
    ax.set_zlim(center[2] - radius, center[2] + radius)
    ax.set_box_aspect((1, 1, 1))


def draw_loads(
    ax,
    nodes,
    loads,
    arrow_fraction=0.16,
    max_arrows=None,
    max_norm=None,
    linewidth=REFERENCE_FORCE_LINEWIDTH,
):
    norms = np.linalg.norm(loads, axis=1)
    valid = np.flatnonzero(norms > 1e-9)
    if max_arrows is not None and len(valid) > max_arrows:
        valid = valid[np.argsort(norms[valid])[-max_arrows:]]

    span = np.max(nodes.max(axis=0) - nodes.min(axis=0))
    max_norm = max_norm if max_norm is not None else (norms[valid].max() if len(valid) else 1.0)
    if max_norm <= 0:
        max_norm = 1.0
    directions = loads[valid] / norms[valid, None]
    relative_lengths = (norms[valid] / max_norm) ** FORCE_DISPLAY_EXPONENT
    scaled = directions * relative_lengths[:, None] * span * arrow_fraction

    ax.quiver(
        nodes[valid, 0],
        nodes[valid, 1],
        nodes[valid, 2],
        scaled[:, 0],
        scaled[:, 1],
        scaled[:, 2],
        length=1.0,
        arrow_length_ratio=0.28,
        linewidth=linewidth,
        color=FORCE_COLOR,
        normalize=False,
    )


def node_resolution_scale(data):
    return REFERENCE_NODE_COUNT / float(data["n_nodes"])


def node_marker_size(data):
    return REFERENCE_NODE_MARKER_SIZE * node_resolution_scale(data)


def force_linewidth(data):
    return REFERENCE_FORCE_LINEWIDTH * np.sqrt(node_resolution_scale(data))


def legend_handles(category_order):
    handles = []
    for category in category_order:
        style = CATEGORY_STYLE[category]
        if category in NODE_AND_LOAD_CATEGORIES - {"Aerodynamic forces"}:
            handles.append(
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    color="none",
                    markerfacecolor=style["color"],
                    markeredgecolor="white",
                    markeredgewidth=0.7,
                    markersize=5.5,
                    label=category,
                )
            )
        else:
            handles.append(
                Line2D(
                    [0],
                    [0],
                    color=style["color"],
                    linewidth=max(style["linewidth"] * 2.0, 1.5),
                    label=category,
                )
            )
    return handles


def draw_nodes(ax, data, model_kind):
    nodes = data["nodes"]
    size = node_marker_size(data)
    if model_kind == "psm":
        fixed = sorted(data["fixed_nodes"])
        pulley = sorted(data["pulley_nodes"] - data["fixed_nodes"])
        free = [
            idx
            for idx in range(len(nodes))
            if idx not in data["fixed_nodes"] and idx not in data["pulley_nodes"]
        ]
        ax.scatter(
            nodes[free, 0],
            nodes[free, 1],
            nodes[free, 2],
            s=size,
            color=NODE_COLOR,
            edgecolors="none",
            linewidths=0.7,
            depthshade=False,
            alpha=0.98,
            zorder=10,
        )
        if pulley:
            ax.scatter(
                nodes[pulley, 0],
                nodes[pulley, 1],
                nodes[pulley, 2],
                s=size,
                color=CATEGORY_STYLE["Pulley node point mass"]["color"],
                edgecolors="none",
                linewidths=0.7,
                depthshade=False,
                alpha=1.0,
                zorder=11,
            )

        if fixed:
            ax.scatter(
                nodes[fixed, 0],
                nodes[fixed, 1],
                nodes[fixed, 2],
                s=size,
                color=CATEGORY_STYLE["KCU / fixed node point mass"]["color"],
                edgecolors="none",
                linewidths=0.9,
                depthshade=False,
                alpha=1.0,
                zorder=12,
            )
    else:
        pulley = sorted(data["pulley_nodes"])
        fixed = sorted(data["fixed_nodes"])
        free = [
            idx
            for idx in range(len(nodes))
            if idx not in data["fixed_nodes"] and idx not in data["pulley_nodes"]
        ]
        ax.scatter(
            nodes[free, 0],
            nodes[free, 1],
            nodes[free, 2],
            s=size,
            color=CATEGORY_STYLE["Free node"]["color"],
            edgecolors="none",
            linewidths=0.7,
            depthshade=False,
            alpha=0.98,
            zorder=6,
        )
        if pulley:
            ax.scatter(
                nodes[pulley, 0],
                nodes[pulley, 1],
                nodes[pulley, 2],
                s=size,
                color=CATEGORY_STYLE["Pulley node"]["color"],
                edgecolors="none",
                linewidths=0.7,
                depthshade=False,
                alpha=1.0,
                zorder=7,
            )
        if fixed:
            ax.scatter(
                nodes[fixed, 0],
                nodes[fixed, 1],
                nodes[fixed, 2],
                s=size,
                color=CATEGORY_STYLE["KCU / fixed node"]["color"],
                edgecolors="none",
                linewidths=0.7,
                depthshade=False,
                alpha=1.0,
                zorder=8,
            )


def draw_model(ax, data, title, category_order, model_kind, load_arrow_limit=None):
    nodes = data["nodes"]
    conn = data["connectivity"]
    categories = data["categories"]
    element_width_scale = data["element_width_scale"]

    for (i, j, *_), category, width_scale in zip(conn, categories, element_width_scale):
        p1 = nodes[int(i)]
        p2 = nodes[int(j)]
        style = CATEGORY_STYLE[category]
        ax.plot(
            [p1[0], p2[0]],
            [p1[1], p2[1]],
            [p1[2], p2[2]],
            color=style["color"],
            linewidth=style["linewidth"] * width_scale,
            alpha=style["alpha"],
        )

    draw_nodes(ax, data, model_kind)
    draw_loads(
        ax,
        nodes,
        data["loads"],
        max_arrows=load_arrow_limit,
        max_norm=data.get("force_arrow_reference"),
        linewidth=force_linewidth(data),
    )
    set_axes_equal(ax, nodes)
    ax.view_init(elev=VIEW_ELEV, azim=VIEW_AZIM)
    ax.set_axis_off()
    ax.set_title(title, fontsize=13, pad=0)
    element_categories = [
        category for category in category_order if category not in NODE_AND_LOAD_CATEGORIES
    ]
    node_load_categories = [
        category for category in category_order if category in NODE_AND_LOAD_CATEGORIES
    ]

    element_legend = ax.legend(
        handles=legend_handles(element_categories),
        loc="lower left",
        bbox_to_anchor=(-0.02, -0.05),
        frameon=False,
        fontsize=7.0,
        handlelength=1.8,
        title="Elements",
        title_fontsize=7.2,
    )
    ax.add_artist(element_legend)
    ax.legend(
        handles=legend_handles(node_load_categories),
        loc="lower left",
        bbox_to_anchor=(0.44, -0.05),
        frameon=False,
        fontsize=7.0,
        handlelength=1.8,
        title="Nodes and loads",
        title_fontsize=7.2,
    )


def save_figure(fig, stem):
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "svg", "png"):
        output_path = OUTPUT_DIR / f"{stem}.{ext}"
        try:
            fig.savefig(
                output_path,
                bbox_inches="tight",
                pad_inches=0.02,
                dpi=300,
                transparent=True,
            )
        except PermissionError:
            fallback_path = OUTPUT_DIR / f"{stem}_updated.{ext}"
            fig.savefig(
                fallback_path,
                bbox_inches="tight",
                pad_inches=0.02,
                dpi=300,
                transparent=True,
            )


def main():
    psm = load_initial_geometry_with_first_load(PSM_RUN)
    fe = load_initial_geometry_with_first_load(FE_RUN)
    equalize_total_force_magnitude(psm, fe)

    fig = plt.figure(figsize=(4.8, 4.2))
    ax = fig.add_subplot(111, projection="3d")
    draw_model(ax, psm, "PSM: Particle System Model", PSM_CATEGORY_ORDER, "psm")
    save_figure(fig, "psm_initial_geometry_with_aero_loads")
    plt.close(fig)

    fig = plt.figure(figsize=(4.8, 4.2))
    ax = fig.add_subplot(111, projection="3d")
    draw_model(ax, fe, "FE: inflatable beam structure", FE_CATEGORY_ORDER, "fe")
    save_figure(fig, "fe_initial_geometry_with_aero_loads")
    plt.close(fig)

    fig = plt.figure(figsize=(9.6, 4.2))
    ax1 = fig.add_subplot(121, projection="3d")
    ax2 = fig.add_subplot(122, projection="3d")
    draw_model(ax1, psm, "PSM: Particle System Model", PSM_CATEGORY_ORDER, "psm")
    draw_model(ax2, fe, "FEM: Finite Element Model", FE_CATEGORY_ORDER, "fe")
    fig.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=0.92, wspace=0.02)
    save_figure(fig, "psm_fe_initial_geometry_with_aero_loads")
    plt.close(fig)

    print(f"Saved AWEC visuals to {OUTPUT_DIR}")
    print(f"PSM loads from timestep {psm['load_step']}; FE loads from timestep {fe['load_step']}.")
    print(
        "Aerodynamic forces scaled to common resultant magnitude: "
        f"{psm['total_force_magnitude']:.2f} N "
        f"(PSM scale {psm['force_scale']:.3f}, FE scale {fe['force_scale']:.3f})."
    )


if __name__ == "__main__":
    main()
