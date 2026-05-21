"""Plot powered photogrammetry shape against the ASKITE shape-validation case.

Run from the ASKITE repo root after running the matching run_ script:

python examples/ch9/shape_validation/plot_ch9_3_2_wing_shape_xyz.py \
  --format pdf,png
"""

import argparse
import os
from pathlib import Path
import sys

os.environ.setdefault("MPLCONFIGDIR", "/tmp/askite_matplotlib")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

PROJECT_DIR = Path(__file__).resolve().parents[3]
SRC_DIR = PROJECT_DIR / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from kitesim.analysis_metrics import compute_geometry_metrics
from kitesim.utils import load_sim_output

CODE_DIR = Path("/home/jellepoland/ownCloud/phd/code")
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

try:
    from dissertation_plot_styling import plot_style
except ImportError:
    from dissertation_plot_styling import set_plot_style as plot_style

DEFAULT_CASE_ID = "shape_validation_va_1675_udp_04151"
DEFAULT_MEASUREMENT_CSV = (
    PROJECT_DIR
    / "data"
    / "ch9"
    / "shape_validation"
    / "Torque_paper_data"
    / "powered_flight_for_depowering_plot.csv"
)

VIEW_SPECS = [
    ("Bottom view", (1, 0), "y [m]", "x [m]"),
    ("Side view", (0, 2), "x [m]", "z [m]"),
    ("Front view", (1, 2), "y [m]", "z [m]"),
]
DEFAULT_LIMITS = {
    "x": (-1.6, 1.6),
    "y": (-4.3, 4.3),
    "z": (-3.0, 0.2),
}


def build_parser():
    parser = argparse.ArgumentParser(
        description="Plot Ch. 9 powered photogrammetry vs ASKITE shape validation."
    )
    parser.add_argument(
        "--case-dir",
        type=Path,
        default=PROJECT_DIR
        / "results"
        / "ch9"
        / "shape_validation"
        / "processed_data"
        / DEFAULT_CASE_ID,
    )
    parser.add_argument("--measurement-csv", type=Path, default=DEFAULT_MEASUREMENT_CSV)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_DIR / "results" / "ch9" / "shape_validation",
    )
    parser.add_argument("--format", default="pdf")
    parser.add_argument("--measurement-label", default="Measurement")
    parser.add_argument("--simulation-label", default="Simulation")
    parser.add_argument("--measurement-color", default="C5")
    parser.add_argument("--simulation-color", default="grey")
    parser.add_argument("--simulation-x-shift", type=float, default=0.0)
    parser.add_argument("--measurement-x-shift", type=float, default=0.0)
    parser.add_argument("--show-grid", action="store_true")
    return parser


def parse_formats(text):
    if text is None:
        return ["pdf"]
    return [item.strip().lower() for item in str(text).split(",") if item.strip()]


def save_figure(fig, output_dir, stem, formats):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for fmt in formats:
        path = output_dir / f"{stem}.{fmt}"
        fig.savefig(path, bbox_inches="tight", dpi=300)
        paths.append(path)
    plt.close(fig)
    return paths


def write_json(path, data):
    import json

    def default(value):
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, (np.integer, np.floating)):
            return value.item()
        if isinstance(value, Path):
            return str(value)
        return str(value)

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2, default=default)


def fit_line_3d(points):
    points = np.asarray(points, dtype=float)
    if points.shape[0] < 2:
        return None
    center = points.mean(axis=0)
    _, _, vt = np.linalg.svd(points - center, full_matrices=False)
    direction = vt[0]
    t = (points - center) @ direction
    ts = np.linspace(float(t.min()), float(t.max()), 100)
    return center + np.outer(ts, direction)


def _pairwise_dist(points):
    diffs = points[:, None, :] - points[None, :, :]
    return np.linalg.norm(diffs, axis=2)


def _path_len(points):
    if points.shape[0] < 2:
        return 0.0
    return float(np.linalg.norm(np.diff(points, axis=0), axis=1).sum())


def _two_opt_locked(order, distances):
    order = list(order)
    if len(order) < 4:
        return order
    improved = True
    while improved:
        improved = False
        for i in range(0, len(order) - 3):
            for j in range(i + 2, len(order) - 1):
                a, b = order[i], order[i + 1]
                c, d = order[j], order[j + 1]
                old = distances[a, b] + distances[c, d]
                new = distances[a, c] + distances[b, d]
                if new + 1e-12 < old:
                    order[i + 1 : j + 1] = reversed(order[i + 1 : j + 1])
                    improved = True
    return order


def _le_path(points, endpoints_k=10, jump_factor=2.5):
    points = np.asarray(points, dtype=float)
    n_points = points.shape[0]
    if n_points <= 1:
        return points.copy()

    distances = _pairwise_dist(points)
    nearest = np.partition(distances + np.eye(n_points) * 1e9, 1, axis=1)[:, 1]
    limited = distances.copy()
    limited[limited > float(np.median(nearest) * jump_factor)] = 1e6

    candidates = np.argsort(points[:, 2])[: max(2, min(endpoints_k, n_points))]
    best_path = None
    best_cost = np.inf
    for i, start in enumerate(candidates):
        for stop in candidates[i + 1 :]:
            used = np.zeros(n_points, dtype=bool)
            used[int(start)] = True
            order = [int(start)]
            while len(order) < n_points - 1:
                remaining = np.where(~used)[0]
                remaining = remaining[remaining != int(stop)]
                nxt = int(remaining[np.argmin(limited[order[-1], remaining])])
                order.append(nxt)
                used[nxt] = True
            order.append(int(stop))
            order = _two_opt_locked(order, limited)
            path = points[np.asarray(order, dtype=int)]
            cost = _path_len(path)
            if cost < best_cost:
                best_path = path
                best_cost = cost
    if best_path is None:
        return points.copy()
    if best_path[0, 2] > best_path[-1, 2]:
        best_path = best_path[::-1]
    return best_path


def _pca_dir(points):
    points = np.asarray(points, dtype=float)
    if points.shape[0] < 2:
        return None
    center = points.mean(axis=0)
    _, _, vt = np.linalg.svd(points - center, full_matrices=False)
    direction = vt[0]
    norm = np.linalg.norm(direction)
    return direction / norm if norm > 0 else None


def _best_fit_plane_normal(points):
    center = points.mean(axis=0)
    _, _, vt = np.linalg.svd(points - center, full_matrices=False)
    normal = vt[-1]
    return normal / max(np.linalg.norm(normal), 1e-12)


def _measurement_frame(groups, le_path, force_flip_x=True):
    strut3 = groups.get("strut3")
    strut4 = groups.get("strut4")
    if strut3 is None or strut4 is None or strut3.shape[0] < 2 or strut4.shape[0] < 2:
        raise ValueError("Powered measurement CSV must contain strut3 and strut4.")

    center_struts = np.vstack([strut3, strut4])
    z_axis = _best_fit_plane_normal(center_struts)
    if np.dot(z_axis, np.array([0.0, 0.0, 1.0])) < 0.0:
        z_axis = -z_axis

    d3 = _pca_dir(strut3)
    d4 = _pca_dir(strut4)
    if d3 is None or d4 is None:
        raise ValueError("Could not infer powered measurement chord axis.")
    if np.dot(d3, d4) < 0.0:
        d4 = -d4
    x_axis = d3 + d4
    x_axis = x_axis - np.dot(x_axis, z_axis) * z_axis
    x_axis = x_axis / max(np.linalg.norm(x_axis), 1e-12)

    if le_path is not None and le_path.shape[0] >= 2:
        le_vec = le_path[-1] - le_path[0]
        le_proj = le_vec - np.dot(le_vec, z_axis) * z_axis
        if np.linalg.norm(le_proj) > 1e-12 and np.dot(x_axis, le_proj) < 0.0:
            x_axis = -x_axis
    if force_flip_x:
        x_axis = -x_axis

    y_axis = np.cross(z_axis, x_axis)
    y_axis = y_axis / max(np.linalg.norm(y_axis), 1e-12)
    z_axis = np.cross(x_axis, y_axis)
    z_axis = z_axis / max(np.linalg.norm(z_axis), 1e-12)
    origin = 0.5 * (strut3.mean(axis=0) + strut4.mean(axis=0))
    axes = np.column_stack([x_axis, y_axis, z_axis])
    return origin, axes


def load_powered_measurement(csv_path, x_shift=0.0):
    df = pd.read_csv(csv_path)
    required = {"group", "idx_in_group", "x", "y", "z"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{csv_path} is missing columns: {sorted(missing)}")

    groups = {}
    for group, sub in df.sort_values(["group", "idx_in_group"]).groupby("group"):
        groups[str(group)] = sub[["x", "y", "z"]].to_numpy(float)
    le_path_world = _le_path(groups["LE"]) if "LE" in groups else None
    origin, axes = _measurement_frame(groups, le_path_world, force_flip_x=True)

    offset = np.array([float(x_shift), 0.0, 0.0])
    local_groups = {
        name: (points - origin) @ axes + offset for name, points in groups.items()
    }
    le_path = (
        None if le_path_world is None else (le_path_world - origin) @ axes + offset
    )
    return {
        "groups": local_groups,
        "le_path": le_path,
        "source_csv": str(csv_path),
        "origin": origin,
        "axes": axes,
    }


def _final_valid_iteration(meta, tracking):
    positions = np.asarray(tracking["positions"])
    n_rows = positions.shape[0]
    n_iter = int(meta.get("n_iter", n_rows))
    return max(0, min(n_iter - 1, n_rows - 1))


def load_simulation_shape(case_dir, x_shift=0.0):
    case_dir = Path(case_dir)
    meta, tracking = load_sim_output(case_dir / "sim_output.h5")
    final_idx = _final_valid_iteration(meta, tracking)
    nodes = np.asarray(tracking["positions"], dtype=float)[final_idx]
    le_indices = np.asarray(meta.get("struc_node_le_indices", []), dtype=int)
    te_indices = np.asarray(meta.get("struc_node_te_indices", []), dtype=int)
    if le_indices.size == 0 or te_indices.size == 0:
        raise ValueError(f"{case_dir} does not contain LE/TE node indices.")

    le_nodes = nodes[le_indices]
    te_nodes = nodes[te_indices]
    midpoints = 0.5 * (le_nodes + te_nodes)
    order = np.argsort(midpoints[:, 1])
    le_indices = le_indices[order]
    te_indices = te_indices[order]
    le_nodes = le_nodes[order]
    te_nodes = te_nodes[order]
    midpoints = midpoints[order]

    y_axis = midpoints[-1] - midpoints[0]
    y_axis = y_axis / max(np.linalg.norm(y_axis), 1e-12)
    x_axis = np.nanmean(te_nodes - le_nodes, axis=0)
    x_axis = x_axis - np.dot(x_axis, y_axis) * y_axis
    x_axis = x_axis / max(np.linalg.norm(x_axis), 1e-12)
    z_axis = np.cross(x_axis, y_axis)
    z_axis = z_axis / max(np.linalg.norm(z_axis), 1e-12)
    y_axis = np.cross(z_axis, x_axis)
    y_axis = y_axis / max(np.linalg.norm(y_axis), 1e-12)

    center_order = np.argsort(np.abs(midpoints[:, 1] - np.nanmedian(midpoints[:, 1])))
    origin = np.nanmean(midpoints[center_order[:2]], axis=0)
    axes = np.column_stack([x_axis, y_axis, z_axis])
    local_nodes = (nodes - origin) @ axes
    local_nodes[:, 0] *= -1.0
    local_nodes = local_nodes + np.array([float(x_shift), 0.0, 0.0])

    groups = {"LE": local_nodes[le_indices]}
    for idx, (le_idx, te_idx) in enumerate(zip(le_indices, te_indices)):
        groups[f"strut{idx}"] = local_nodes[[le_idx, te_idx]]

    metrics = compute_geometry_metrics(local_nodes, le_indices, te_indices)
    metrics.update(
        {
            "case_dir": str(case_dir),
            "source_h5": str(case_dir / "sim_output.h5"),
            "final_iteration_index": int(final_idx),
            "le_indices": le_indices.tolist(),
            "te_indices": te_indices.tolist(),
            "local_x_reversed_for_plot": True,
        }
    )
    return {
        "groups": groups,
        "le_path": groups["LE"],
        "metrics": metrics,
        "origin": origin,
        "axes": axes,
        "meta": meta,
    }


def _project(points, proj):
    points = np.asarray(points, dtype=float)
    return points[:, [proj[0], proj[1]]]


def _draw_shape(
    ax,
    shape,
    view_proj,
    color,
    label,
    linewidth,
    alpha,
    zorder,
    marker_size=3.0,
):
    handle = None
    le_points = shape["groups"].get("LE")
    if le_points is not None and len(le_points) > 0:
        projected = _project(le_points, view_proj)
        ax.scatter(
            projected[:, 0],
            projected[:, 1],
            s=marker_size**2,
            color=color,
            alpha=alpha,
            linewidths=0,
            zorder=zorder + 0.1,
        )

    le_path = shape.get("le_path")
    if le_path is not None and len(le_path) > 1:
        projected = _project(le_path, view_proj)
        (handle,) = ax.plot(
            projected[:, 0],
            projected[:, 1],
            color=color,
            lw=linewidth,
            alpha=alpha,
            label=label,
            marker="o",
            markersize=marker_size,
            markerfacecolor=color,
            markeredgewidth=0,
            zorder=zorder,
        )

    for name, points in shape["groups"].items():
        if name == "LE" or len(points) == 0:
            continue
        point_projection = _project(points, view_proj)
        ax.scatter(
            point_projection[:, 0],
            point_projection[:, 1],
            s=(marker_size * 0.85) ** 2,
            color=color,
            alpha=alpha,
            linewidths=0,
            zorder=zorder + 0.1,
        )
        line = fit_line_3d(points)
        if line is None:
            continue
        projected = _project(line, view_proj)
        ax.plot(
            projected[:, 0],
            projected[:, 1],
            color=color,
            lw=linewidth * 0.75,
            alpha=alpha,
            marker="o" if len(points) == 2 else None,
            markersize=marker_size * 0.85,
            markerfacecolor=color,
            markeredgewidth=0,
            zorder=zorder,
        )
    return handle


def _all_points(*shapes):
    points = []
    for shape in shapes:
        for group_points in shape["groups"].values():
            points.append(np.asarray(group_points, dtype=float))
        if shape.get("le_path") is not None:
            points.append(np.asarray(shape["le_path"], dtype=float))
    return np.vstack(points)


def _expanded_limits(points, limits):
    expanded = dict(limits)
    axis_names = ["x", "y", "z"]
    for axis, name in enumerate(axis_names):
        lo, hi = expanded[name]
        axis_values = points[:, axis]
        data_lo = float(np.nanmin(axis_values))
        data_hi = float(np.nanmax(axis_values))
        lo = min(lo, data_lo)
        hi = max(hi, data_hi)
        if lo == hi:
            pad = 0.1
        else:
            pad = 0.04 * (hi - lo)
        expanded[name] = (lo - pad, hi + pad)
    return expanded


def _view_limits(view_proj, limits):
    axis_names = ["x", "y", "z"]
    return limits[axis_names[view_proj[0]]], limits[axis_names[view_proj[1]]]


def main():
    args = build_parser().parse_args()
    plot_style()
    formats = parse_formats(args.format)
    measurement = load_powered_measurement(
        args.measurement_csv, x_shift=args.measurement_x_shift
    )
    simulation = load_simulation_shape(args.case_dir, x_shift=args.simulation_x_shift)
    limits = _expanded_limits(_all_points(measurement, simulation), DEFAULT_LIMITS)

    fig = plt.figure(figsize=(10.0, 3.0))
    span_x = limits["x"][1] - limits["x"][0]
    span_y = limits["y"][1] - limits["y"][0]
    width_ratios = [span_y, span_x, span_y]
    gs = fig.add_gridspec(
        1,
        3,
        width_ratios=width_ratios,
        left=0.06,
        right=0.98,
        bottom=0.18,
        top=0.86,
        wspace=0.18,
    )
    axes = [fig.add_subplot(gs[0, idx]) for idx in range(3)]

    for ax, (title, proj, xlabel, ylabel) in zip(axes, VIEW_SPECS):
        _draw_shape(
            ax,
            simulation,
            proj,
            args.simulation_color,
            args.simulation_label,
            linewidth=3.2,
            alpha=0.5,
            zorder=1,
            marker_size=6.0,
        )
        _draw_shape(
            ax,
            measurement,
            proj,
            args.measurement_color,
            args.measurement_label,
            linewidth=2.0,
            alpha=0.95,
            zorder=3,
        )
        xlim, ylim = _view_limits(proj, limits)
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if args.show_grid:
            ax.grid(True, color="0.88", linewidth=0.5)
        else:
            ax.grid(False)

    handles = [
        Line2D(
            [0],
            [0],
            color=args.measurement_color,
            lw=2.0,
            label=args.measurement_label,
        ),
        Line2D(
            [0],
            [0],
            color=args.simulation_color,
            lw=3.2,
            label=args.simulation_label,
        ),
    ]
    axes[-1].legend(handles=handles, loc="lower center", frameon=True)

    output_paths = save_figure(
        fig,
        args.output_dir,
        "fig_9_3_2_3_powered_shape_vs_askite",
        formats,
    )
    write_json(
        args.output_dir / "fig_9_3_2_3_powered_shape_vs_askite_metrics.json",
        {
            "case_dir": str(args.case_dir),
            "measurement_csv": str(args.measurement_csv),
            "output_files": output_paths,
            "axis_limits": limits,
            "simulation_metrics": simulation["metrics"],
        },
    )
    print("created files:")
    for path in output_paths:
        print(f"- {path}")


if __name__ == "__main__":
    main()
