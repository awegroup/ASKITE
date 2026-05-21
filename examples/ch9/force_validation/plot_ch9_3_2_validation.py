"""How to run (from ASKITE repo root):

python examples/ch9/force_validation/plot_ch9_3_2_validation.py \
    --ekf-dir /home/jellepoland/ownCloud/phd/code/EKF-AWE/data/ch9_3_2_straight_vwt \
    --askite-summary results/ch9/force_validation/processed_data/ch9_3_2_askite_case_summary.csv \
    --output-dir results/ch9/force_validation \
    --format pdf,png
"""

import argparse
import os
from pathlib import Path
import sys

os.environ.setdefault("MPLCONFIGDIR", "/tmp/askite_matplotlib")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
import numpy as np
import pandas as pd

CH9_DIR = Path(__file__).resolve().parents[1]
if str(CH9_DIR) not in sys.path:
    sys.path.insert(0, str(CH9_DIR))

from ch9_analysis_utils import PROJECT_DIR, infer_column, parse_formats, save_figure

CODE_DIR = Path("/home/jellepoland/ownCloud/phd/code")
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

try:
    from dissertation_plot_styling import plot_style
except ImportError:
    from dissertation_plot_styling import set_plot_style as plot_style

DEFAULT_EKF_DIR = (
    Path("/home/jellepoland/ownCloud/phd/code/EKF-AWE")
    / "data"
    / "ch9_3_2_straight_vwt"
)


def build_parser():
    parser = argparse.ArgumentParser(
        description="Plot Section 9.3.2 EKF/ASKITE validation figures."
    )
    parser.add_argument("--ekf-dir", type=Path, default=DEFAULT_EKF_DIR)
    parser.add_argument(
        "--askite-summary",
        type=Path,
        default=PROJECT_DIR
        / "results"
        / "ch9"
        / "force_validation"
        / "processed_data"
        / "ch9_3_2_askite_case_summary.csv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_DIR / "results" / "ch9" / "force_validation",
    )
    parser.add_argument("--format", default="pdf")
    parser.add_argument("--campaigns", default="2019,2025")
    parser.add_argument(
        "--force-source",
        choices=["wing_only", "aero_total", "reaction_total"],
        default="aero_total",
    )
    parser.add_argument(
        "--ekf-occurrence-source",
        choices=["samples", "windows"],
        default="samples",
        help="Raw EKF-AWE table used for occurrence heatmaps.",
    )
    parser.add_argument(
        "--ekf-occurrence-file",
        type=Path,
        default=None,
        help="Optional explicit EKF-AWE occurrence CSV. Overrides --ekf-occurrence-source.",
    )
    parser.add_argument(
        "--heatmap-bins",
        type=int,
        default=45,
        help="Number of bins per axis for EKF occurrence heatmaps.",
    )
    return parser


def _read_required(path):
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def _campaign_groups(df, campaigns):
    if "campaign" in df.columns:
        return [
            (camp, df[df["campaign"].astype(str) == str(camp)]) for camp in campaigns
        ]
    if "year" in df.columns:
        return [(camp, df[df["year"].astype(str) == str(camp)]) for camp in campaigns]
    return [("all", df)]


def _plot_binned(ax, df, x_candidates, y_candidates, label, color):
    x_col = infer_column(df, x_candidates)
    y_col = infer_column(df, y_candidates)
    if x_col is None or y_col is None:
        return False
    x = pd.to_numeric(df[x_col], errors="coerce")
    y = pd.to_numeric(df[y_col], errors="coerce")
    ax.plot(x, y, "-", color=color, linewidth=1.2, label=label)
    low_col = infer_column(
        df,
        [
            y_col.replace("_mean", "_p10"),
            y_col.replace("_mean", "_q10"),
            f"{y_col}_p10",
            f"{y_col}_lo",
        ],
    )
    high_col = infer_column(
        df,
        [
            y_col.replace("_mean", "_p90"),
            y_col.replace("_mean", "_q90"),
            f"{y_col}_p90",
            f"{y_col}_hi",
        ],
    )
    if low_col and high_col:
        ax.fill_between(
            x,
            pd.to_numeric(df[low_col], errors="coerce"),
            pd.to_numeric(df[high_col], errors="coerce"),
            color=color,
            alpha=0.15,
            linewidth=0,
        )
    return True


def _plot_sim(ax, summary, x_col, y_col, label, color, marker):
    if x_col not in summary or y_col not in summary:
        return
    ax.plot(
        pd.to_numeric(summary[x_col], errors="coerce"),
        pd.to_numeric(summary[y_col], errors="coerce"),
        linestyle="none",
        marker=marker,
        markersize=4,
        color=color,
        label=label,
    )


def _filter_converged(summary):
    if "converged" not in summary.columns:
        return summary
    return summary[summary["converged"].astype(str).str.lower().isin(["true", "1"])]


def _force_col(force_source):
    return {
        "wing_only": "sim_wing_force_N",
        "aero_total": "sim_aero_total_force_N",
        "reaction_total": "sim_tether_or_reaction_force_N",
    }[force_source]


def _read_ekf_occurrences(args):
    if args.ekf_occurrence_file is not None:
        path = args.ekf_occurrence_file
    elif args.ekf_occurrence_source == "samples":
        path = args.ekf_dir / "ch9_3_2_straight_samples.csv"
    else:
        path = args.ekf_dir / "ch9_3_2_straight_windows.csv"

    df = _read_required(path)
    if "straight_filter_pass" in df.columns:
        passed = df["straight_filter_pass"].astype(str).str.lower().isin(["true", "1"])
        df = df[passed].copy()
    return df


def _numeric(df, column):
    if column is None or column not in df.columns:
        return pd.Series(dtype=float)
    return pd.to_numeric(df[column], errors="coerce")


def _finite_values(*series):
    values = []
    for item in series:
        arr = np.asarray(item, dtype=float).reshape(-1)
        arr = arr[np.isfinite(arr)]
        if arr.size:
            values.append(arr)
    if not values:
        return np.array([], dtype=float)
    return np.concatenate(values)


def _robust_range(values, lower=1.0, upper=99.0):
    values = np.asarray(values, dtype=float).reshape(-1)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return (0.0, 1.0)
    lo, hi = np.nanpercentile(values, [lower, upper])
    if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
        lo = float(np.nanmin(values))
        hi = float(np.nanmax(values))
    if lo == hi:
        pad = abs(lo) * 0.05 if lo != 0 else 1.0
        return (lo - pad, hi + pad)
    pad = 0.04 * (hi - lo)
    return (float(lo - pad), float(hi + pad))


def _robust_range_with_required(background_values, required_values):
    lo, hi = _robust_range(background_values)
    required = np.asarray(required_values, dtype=float).reshape(-1)
    required = required[np.isfinite(required)]
    if required.size == 0:
        return lo, hi
    lo = min(lo, float(np.nanmin(required)))
    hi = max(hi, float(np.nanmax(required)))
    if lo == hi:
        pad = abs(lo) * 0.05 if lo != 0 else 1.0
    else:
        pad = 0.04 * (hi - lo)
    return (float(lo - pad), float(hi + pad))


def _campaign_df(df, campaign):
    if "campaign" in df.columns:
        return df[df["campaign"].astype(str) == str(campaign)]
    if "year" in df.columns:
        return df[df["year"].astype(str) == str(campaign)]
    return df


def _hist2d_for_panel(df, x_col, y_col, x_range, y_range, bins):
    x = _numeric(df, x_col)
    y = _numeric(df, y_col)
    mask = np.isfinite(x) & np.isfinite(y)
    if not mask.any():
        return None
    hist, xedges, yedges = np.histogram2d(
        x[mask],
        y[mask],
        bins=[bins, bins],
        range=[x_range, y_range],
    )
    return hist, xedges, yedges


def _heatmap_norm(max_count):
    if max_count > 1:
        return LogNorm(vmin=1, vmax=max_count)
    return Normalize(vmin=0, vmax=1)


def _plot_occurrence_heatmap(
    ax,
    hist_result,
    norm,
    cmap="Blues",
):
    if hist_result is None:
        return None
    hist, xedges, yedges = hist_result
    masked = np.ma.masked_where(hist.T <= 0, hist.T)
    return ax.pcolormesh(xedges, yedges, masked, cmap=cmap, norm=norm, shading="auto")


def _plot_askite_points(ax, df, x_col, y_col, color, marker, label):
    if x_col not in df.columns or y_col not in df.columns:
        return
    x = _numeric(df, x_col)
    y = _numeric(df, y_col)
    mask = np.isfinite(x) & np.isfinite(y)
    if not mask.any():
        return
    ax.plot(
        x[mask],
        y[mask],
        linestyle="none",
        marker=marker,
        markersize=4.5,
        markerfacecolor=color,
        markeredgewidth=0.7,
        color=color,
        label=label,
    )


def _plot_force_coefficient_heatmap_grid(
    ekf_occurrences,
    summary,
    campaigns,
    args,
    x_kind,
    x_candidates,
    sim_x_col,
    x_label,
    output_stem,
    formats,
):
    panels = [
        (
            "Kite tether force",
            [
                "tether_force_kite_N",
                "mean_tether_force_kite_N",
                "tether_force_kite_N_mean",
                "preferred_tether_force_N",
                "mean_preferred_tether_force_N",
                "tether_force_preferred_N_mean",
                "ground_tether_force_N",
                "mean_ground_tether_force_N",
            ],
            [_force_col(args.force_source)],
            "Kite-side tether force [N]",
        ),
        (
            "$C_{L,\\mathrm{kite}}$",
            [
                "CL_kite_ekf",
                "mean_CL_kite_ekf",
                "C_L_kite",
                "CL_ekf",
                "mean_CL_ekf",
                "C_L_mean",
                "CL_mean",
                "CL_ekf_mean",
            ],
            ["sim_CL_kite", "sim_CL_total", "sim_CL_wing"],
            "$C_{L,\\mathrm{kite}}$ [-]",
        ),
        (
            "$C_{D,\\mathrm{kite}}$",
            [
                "CD_kite_ekf",
                "mean_CD_kite_ekf",
                "C_D_kite",
                "CD_ekf",
                "mean_CD_ekf",
                "C_D_mean",
                "CD_mean",
                "CD_ekf_mean",
            ],
            ["sim_CD_kite", "sim_CD_total", "sim_CD_wing"],
            "$C_{D,\\mathrm{kite}}$ [-]",
        ),
        (
            "$C_{L,\\mathrm{kite}}/C_{D,\\mathrm{kite}}$",
            [
                "L_over_D_kite_ekf",
                "mean_L_over_D_kite_ekf",
                "L_over_D_ekf",
                "mean_L_over_D_ekf",
                "mean_CL_over_CD_from_window_means",
                "L_over_D_mean",
                "CL_over_CD_mean",
                "glide_ratio_mean",
            ],
            ["sim_L_over_D_kite", "sim_L_over_D_total", "sim_L_over_D_wing"],
            "$C_{L,\\mathrm{kite}}/C_{D,\\mathrm{kite}}$ [-]",
        ),
    ]

    x_col = infer_column(ekf_occurrences, x_candidates, required=True)
    x_range = _robust_range_with_required(
        _numeric(ekf_occurrences, x_col),
        _numeric(summary, sim_x_col),
    )

    resolved_panels = []
    hist_results = {}
    max_count = 0
    for panel_idx, (title, y_candidates, sim_y_candidates, y_label) in enumerate(panels):
        y_col = infer_column(ekf_occurrences, y_candidates, required=True)
        sim_y_col = infer_column(summary, sim_y_candidates)
        y_range = _robust_range_with_required(
            _numeric(ekf_occurrences, y_col),
            _numeric(summary, sim_y_col),
        )
        resolved_panels.append((title, y_col, sim_y_col, y_label, y_range))
        for campaign in campaigns:
            hist_result = _hist2d_for_panel(
                _campaign_df(ekf_occurrences, campaign),
                x_col,
                y_col,
                x_range,
                y_range,
                args.heatmap_bins,
            )
            hist_results[(str(campaign), panel_idx)] = hist_result
            if hist_result is not None:
                max_count = max(max_count, int(np.nanmax(hist_result[0])))

    norm = _heatmap_norm(max_count)
    fig, axes = plt.subplots(
        len(campaigns),
        4,
        figsize=(13.2, 5.9),
        sharex=True,
        squeeze=False,
        constrained_layout=True,
    )
    last_mesh = None

    for row_idx, campaign in enumerate(campaigns):
        summary_campaign = _campaign_df(summary, campaign)
        for col_idx, (title, _, sim_y_col, y_label, y_range) in enumerate(
            resolved_panels
        ):
            ax = axes[row_idx, col_idx]
            mesh = _plot_occurrence_heatmap(
                ax,
                hist_results[(str(campaign), col_idx)],
                norm,
            )
            if mesh is not None:
                last_mesh = mesh
            _plot_askite_points(
                ax,
                summary_campaign,
                sim_x_col,
                sim_y_col,
                "black",
                "o",
                "Simulation" if col_idx == 0 else "_nolegend_",
            )
            ax.set_xlim(x_range)
            ax.set_ylim(y_range)
            ax.grid(True, color="0.88", linewidth=0.45)
            ax.set_title(title if row_idx == 0 else "")
            ax.set_ylabel(y_label)
            if col_idx == 0:
                handles, _ = ax.get_legend_handles_labels()
                if handles:
                    ax.legend(frameon=False, fontsize=7, loc="best")
        axes[row_idx, 0].set_ylabel(f"{campaign}\n{resolved_panels[0][3]}")

    for ax in axes[-1, :]:
        ax.set_xlabel(x_label)
    if last_mesh is not None:
        fig.colorbar(
            last_mesh,
            ax=axes.ravel().tolist(),
            label="EKF occurrence count",
            fraction=0.025,
            pad=0.015,
        )
    fig.suptitle(f"Force and coefficients vs {x_kind}")
    save_figure(fig, args.output_dir, output_stem, formats)


def main():
    args = build_parser().parse_args()
    plot_style()
    formats = parse_formats(args.format)
    campaigns = [item.strip() for item in args.campaigns.split(",") if item.strip()]
    ekf_occurrences = _read_ekf_occurrences(args)
    summary = _filter_converged(pd.read_csv(args.askite_summary))

    _plot_force_coefficient_heatmap_grid(
        ekf_occurrences=ekf_occurrences,
        summary=summary,
        campaigns=campaigns,
        args=args,
        x_kind="$V_a$",
        x_candidates=[
            "V_a_ms",
            "mean_V_a_ms",
            "Va_ms",
            "va_ms",
            "V_a_bin_center_ms",
            "V_a_bin_center",
            "va_bin_center_ms",
            "bin_center",
        ],
        sim_x_col="requested_V_a_ms",
        x_label="$V_a$ [m s$^{-1}$]",
        output_stem="fig_9_3_2_2_force_coefficients_vs_va",
        formats=formats,
    )

    _plot_force_coefficient_heatmap_grid(
        ekf_occurrences=ekf_occurrences,
        summary=summary,
        campaigns=campaigns,
        args=args,
        x_kind="tether length",
        x_candidates=[
            "tether_length_m",
            "mean_tether_length_m",
            "tether_length_bin_center_m",
            "tether_length_bin_center",
            "bin_center",
        ],
        sim_x_col="requested_tether_length_m",
        x_label="Tether length [m]",
        output_stem="fig_9_3_2_3_force_coefficients_vs_tether_length",
        formats=formats,
    )

    if {
        "requested_V_a_ms",
        "measured_tether_force_N",
        _force_col(args.force_source),
    }.issubset(summary.columns):
        fig, ax = plt.subplots(figsize=(4.2, 3.1))
        ax.plot(
            pd.to_numeric(summary["requested_V_a_ms"], errors="coerce"),
            pd.to_numeric(summary[_force_col(args.force_source)], errors="coerce")
            - pd.to_numeric(summary["measured_tether_force_N"], errors="coerce"),
            "o",
            color="0.25",
        )
        ax.set_xlabel("$V_a$ [m s$^{-1}$]")
        ax.set_ylabel("ASKITE - EKF force [N]")
        ax.grid(True, color="0.88", linewidth=0.6)
        fig.tight_layout()
        save_figure(fig, args.output_dir, "sim_minus_ekf_force_vs_va", formats)

    sim_cl_col = infer_column(summary, ["sim_CL_kite", "sim_CL_total", "sim_CL_wing"])
    sim_cd_col = infer_column(summary, ["sim_CD_kite", "sim_CD_total", "sim_CD_wing"])
    if (
        "requested_V_a_ms" in summary.columns
        and sim_cl_col is not None
        and sim_cd_col is not None
        and {"measured_CL_ekf", "measured_CD_ekf"}.issubset(summary.columns)
    ):
        fig, ax = plt.subplots(figsize=(4.2, 3.1))
        x = pd.to_numeric(summary["requested_V_a_ms"], errors="coerce")
        ax.plot(
            x,
            pd.to_numeric(summary[sim_cl_col], errors="coerce")
            - pd.to_numeric(summary["measured_CL_ekf"], errors="coerce"),
            "o",
            label="$C_L$",
        )
        ax.plot(
            x,
            pd.to_numeric(summary[sim_cd_col], errors="coerce")
            - pd.to_numeric(summary["measured_CD_ekf"], errors="coerce"),
            "s",
            label="$C_D$",
        )
        ax.set_xlabel("$V_a$ [m s$^{-1}$]")
        ax.set_ylabel("ASKITE - EKF kite coefficient")
        ax.grid(True, color="0.88", linewidth=0.6)
        ax.legend(frameon=False)
        fig.tight_layout()
        save_figure(fig, args.output_dir, "sim_minus_ekf_CL_CD_vs_va", formats)

    sim_glide_col = infer_column(
        summary, ["sim_L_over_D_kite", "sim_L_over_D_total", "sim_L_over_D_wing"]
    )
    if (
        "requested_tether_length_m" in summary.columns
        and sim_glide_col is not None
        and "measured_L_over_D_ekf" in summary.columns
    ):
        fig, ax = plt.subplots(figsize=(4.2, 3.1))
        ax.plot(
            pd.to_numeric(summary["requested_tether_length_m"], errors="coerce"),
            pd.to_numeric(summary[sim_glide_col], errors="coerce")
            - pd.to_numeric(summary["measured_L_over_D_ekf"], errors="coerce"),
            "o",
            color="0.25",
        )
        ax.set_xlabel("Tether length [m]")
        ax.set_ylabel("ASKITE - EKF kite $C_L/C_D$")
        ax.grid(True, color="0.88", linewidth=0.6)
        fig.tight_layout()
        save_figure(
            fig, args.output_dir, "sim_minus_ekf_glide_vs_tether_length", formats
        )


if __name__ == "__main__":
    main()
