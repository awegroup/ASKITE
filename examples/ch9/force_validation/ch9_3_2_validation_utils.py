from pathlib import Path
import sys

import numpy as np
import pandas as pd

CH9_DIR = Path(__file__).resolve().parents[1]
if str(CH9_DIR) not in sys.path:
    sys.path.insert(0, str(CH9_DIR))

from ch9_analysis_utils import (
    equalize_2d_axes,
    final_valid_iteration,
    infer_column,
    load_case as load_askite_case,
    save_figure,
)
from kitesim.analysis_metrics import (
    compute_geometry_metrics,
    compute_projected_area_from_le_te,
    compute_projected_span,
)


def load_ekf_harvest(ekf_dir):
    ekf_dir = Path(ekf_dir)
    return {
        "cases": pd.read_csv(ekf_dir / "ch9_3_2_vwt_cases_for_askite.csv"),
        "binned_by_va": pd.read_csv(ekf_dir / "ch9_3_2_binned_by_va.csv"),
        "binned_by_tether_length": pd.read_csv(
            ekf_dir / "ch9_3_2_binned_by_tether_length.csv"
        ),
    }


def load_askite_summary(path):
    return pd.read_csv(path)


def project_nodes(nodes, frame="solver"):
    if frame not in ("solver", "body", "wind"):
        raise ValueError("frame must be solver, body, or wind")
    return np.asarray(nodes, dtype=float)


def plot_binned_trend(ax, bins, x_candidates, y_candidates, label, color):
    x_col = infer_column(bins, x_candidates)
    y_col = infer_column(bins, y_candidates)
    if x_col is None or y_col is None:
        return False
    x = pd.to_numeric(bins[x_col], errors="coerce")
    y = pd.to_numeric(bins[y_col], errors="coerce")
    ax.plot(x, y, "-", color=color, linewidth=1.2, label=label)
    low_col = infer_column(bins, [f"{y_col}_p10", f"{y_col}_q10", f"{y_col}_lo"])
    high_col = infer_column(bins, [f"{y_col}_p90", f"{y_col}_q90", f"{y_col}_hi"])
    if low_col and high_col:
        ax.fill_between(
            x,
            pd.to_numeric(bins[low_col], errors="coerce"),
            pd.to_numeric(bins[high_col], errors="coerce"),
            color=color,
            alpha=0.15,
            linewidth=0,
        )
    return True


def plot_sim_markers(ax, summary, x_col, y_col, label, color, marker="o"):
    if x_col not in summary or y_col not in summary:
        return False
    ax.plot(
        pd.to_numeric(summary[x_col], errors="coerce"),
        pd.to_numeric(summary[y_col], errors="coerce"),
        linestyle="none",
        marker=marker,
        color=color,
        markersize=4,
        label=label,
    )
    return True
