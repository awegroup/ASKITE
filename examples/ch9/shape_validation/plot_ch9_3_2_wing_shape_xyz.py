"""How to run (from ASKITE repo root):

Single case:
python examples/ch9/shape_validation/plot_ch9_3_2_wing_shape_xyz.py \
    --case-dir results/ch9/force_validation/askite_vwt_001 \
    --output-dir results/ch9/shape_validation/plots \
    --format pdf,png

Batch from summary CSV:
python examples/ch9/shape_validation/plot_ch9_3_2_wing_shape_xyz.py \
    --case-list results/ch9/force_validation/ch9_3_2_askite_case_summary.csv \
    --output-dir results/ch9/shape_validation/plots \
    --overlay-initial
"""

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

CH9_DIR = Path(__file__).resolve().parents[1]
if str(CH9_DIR) not in sys.path:
    sys.path.insert(0, str(CH9_DIR))

from ch9_analysis_utils import (
    PROJECT_DIR,
    equalize_2d_axes,
    final_valid_iteration,
    load_case,
    parse_formats,
    save_figure,
    write_json,
)
from kitesim.analysis_metrics import compute_geometry_metrics


def build_parser():
    parser = argparse.ArgumentParser(
        description="Plot final ASKITE wing shape in x-y, x-z, and y-z projections."
    )
    parser.add_argument("--case-dir", type=Path, default=None)
    parser.add_argument(
        "--case-list",
        type=Path,
        default=PROJECT_DIR
        / "results"
        / "ch9"
        / "force_validation"
        / "ch9_3_2_askite_case_summary.csv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_DIR / "results" / "ch9" / "shape_validation" / "plots",
    )
    parser.add_argument("--frame", choices=["solver", "body", "wind"], default="solver")
    parser.add_argument("--no-plot-connectivity", action="store_true")
    parser.add_argument("--plot-node-labels", action="store_true")
    parser.add_argument("--overlay-initial", action="store_true")
    parser.add_argument("--overlay-photogrammetry", type=Path, default=None)
    parser.add_argument("--format", default="pdf")
    return parser


def _cases_from_args(args):
    if args.case_dir is not None:
        return [(args.case_dir.name, args.case_dir)]
    if args.case_list is None:
        raise SystemExit("Provide --case-dir or --case-list.")
    df = pd.read_csv(args.case_list)
    dir_col = "results_dir" if "results_dir" in df.columns else "case_dir"
    id_col = "case_id" if "case_id" in df.columns else dir_col
    return [(str(row[id_col]), Path(row[dir_col])) for _, row in df.iterrows()]


def _body_frame(nodes, le_indices, te_indices):
    nodes = np.asarray(nodes, dtype=float)
    wing_idx = np.unique(np.concatenate([le_indices, te_indices])).astype(int)
    origin = np.nanmean(nodes[wing_idx], axis=0)
    left = nodes[wing_idx[np.argmin(nodes[wing_idx, 1])]]
    right = nodes[wing_idx[np.argmax(nodes[wing_idx, 1])]]
    span_axis = right - left
    span_axis /= max(np.linalg.norm(span_axis), 1e-12)
    chord_axis = np.nanmean(nodes[te_indices] - nodes[le_indices], axis=0)
    chord_axis = chord_axis - np.dot(chord_axis, span_axis) * span_axis
    chord_axis /= max(np.linalg.norm(chord_axis), 1e-12)
    normal_axis = np.cross(chord_axis, span_axis)
    normal_axis /= max(np.linalg.norm(normal_axis), 1e-12)
    rotation = np.vstack([chord_axis, span_axis, normal_axis])
    return (nodes - origin) @ rotation.T


def _project(nodes, frame, le_indices, te_indices):
    if frame == "body":
        return _body_frame(nodes, le_indices, te_indices)
    return np.asarray(nodes, dtype=float)


def _plot_edges(ax, nodes2d, connectivity, color="0.55", linewidth=0.5, linestyle="-"):
    for ci, cj in np.asarray(connectivity, dtype=int):
        if ci < len(nodes2d) and cj < len(nodes2d):
            ax.plot(
                [nodes2d[ci, 0], nodes2d[cj, 0]],
                [nodes2d[ci, 1], nodes2d[cj, 1]],
                color=color,
                linewidth=linewidth,
                linestyle=linestyle,
                zorder=1,
            )


def _make_shape_figure(
    case_id,
    case_dir,
    output_dir,
    frame,
    formats,
    plot_connectivity=True,
    plot_node_labels=False,
    overlay_initial=False,
):
    meta, tracking, final_idx = load_case(case_dir)
    nodes = np.asarray(tracking["positions"])[final_idx]
    initial_nodes = np.asarray(tracking["positions"])[0]
    connectivity = np.asarray(meta.get("kite_connectivity", []), dtype=int)
    le_indices = np.asarray(meta.get("struc_node_le_indices", []), dtype=int)
    te_indices = np.asarray(meta.get("struc_node_te_indices", []), dtype=int)
    if le_indices.size == 0 or te_indices.size == 0:
        le_indices = np.arange(1, min(len(nodes), 21), 2)
        te_indices = np.arange(2, min(len(nodes), 21), 2)

    nodes_frame = _project(nodes, frame, le_indices, te_indices)
    initial_frame = _project(initial_nodes, frame, le_indices, te_indices)

    views = [
        ((0, 1), "x", "y", "Top"),
        ((0, 2), "x", "z", "Side"),
        ((1, 2), "y", "z", "Front"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(9.0, 3.2))
    wing_indices = np.unique(np.concatenate([le_indices, te_indices])).astype(int)
    for ax, (dims, xlabel, ylabel, title) in zip(axes, views):
        projected = nodes_frame[:, dims]
        if overlay_initial:
            initial_projected = initial_frame[:, dims]
            if plot_connectivity and connectivity.size:
                _plot_edges(
                    ax,
                    initial_projected,
                    connectivity,
                    color="0.78",
                    linewidth=0.5,
                    linestyle="--",
                )
            ax.plot(
                initial_projected[wing_indices, 0],
                initial_projected[wing_indices, 1],
                ".",
                color="0.65",
                markersize=2,
                label="Initial",
            )
        if plot_connectivity and connectivity.size:
            _plot_edges(ax, projected, connectivity, color="0.45", linewidth=0.5)
        ax.plot(
            projected[wing_indices, 0],
            projected[wing_indices, 1],
            ".",
            color="0.05",
            markersize=3,
            label="Final",
        )
        if plot_node_labels:
            for node_id in wing_indices:
                ax.text(
                    projected[node_id, 0],
                    projected[node_id, 1],
                    str(node_id),
                    fontsize=5,
                    color="0.2",
                )
        ax.set_title(title)
        ax.set_xlabel(f"{xlabel} [m]")
        ax.set_ylabel(f"{ylabel} [m]")
        equalize_2d_axes(ax, projected[wing_indices, 0], projected[wing_indices, 1])
        ax.grid(True, color="0.9", linewidth=0.5)
    axes[0].legend(frameon=False, fontsize=7)
    fig.tight_layout()

    stem = f"fig_9_3_2_3_shape_case_{case_id}_xyz"
    save_figure(fig, output_dir, stem, formats)
    metrics = compute_geometry_metrics(nodes, le_indices, te_indices)
    metrics.update(
        {
            "case_id": case_id,
            "case_dir": str(case_dir),
            "frame": frame,
            "final_iteration_index": int(final_idx),
            "le_indices": le_indices.tolist(),
            "te_indices": te_indices.tolist(),
        }
    )
    write_json(output_dir / f"{stem}_metrics.json", metrics)


def main():
    args = build_parser().parse_args()
    formats = parse_formats(args.format)
    for case_id, case_dir in _cases_from_args(args):
        _make_shape_figure(
            case_id=case_id,
            case_dir=case_dir,
            output_dir=args.output_dir,
            frame=args.frame,
            formats=formats,
            plot_connectivity=not args.no_plot_connectivity,
            plot_node_labels=args.plot_node_labels,
            overlay_initial=args.overlay_initial,
        )


if __name__ == "__main__":
    main()
