"""How to run (from ASKITE repo root):

python examples/ch9/force_validation/plot_ch9_3_2_validation.py \
    --ekf-dir /home/jellepoland/ownCloud/phd/code/EKF-AWE/data/ch9_3_2_straight_vwt \
    --askite-summary results/ch9/force_validation/ch9_3_2_askite_case_summary.csv \
    --output-dir results/ch9/force_validation/plots \
    --format pdf,png
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

from ch9_analysis_utils import PROJECT_DIR, infer_column, parse_formats, save_figure

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
        / "ch9_3_2_askite_case_summary.csv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_DIR / "results" / "ch9" / "force_validation" / "plots",
    )
    parser.add_argument("--format", default="pdf")
    parser.add_argument("--campaigns", default="2019,2025")
    parser.add_argument(
        "--force-source",
        choices=["wing_only", "aero_total", "reaction_total"],
        default="wing_only",
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


def main():
    args = build_parser().parse_args()
    formats = parse_formats(args.format)
    campaigns = [item.strip() for item in args.campaigns.split(",") if item.strip()]
    binned_va = _read_required(args.ekf_dir / "ch9_3_2_binned_by_va.csv")
    binned_tether = _read_required(args.ekf_dir / "ch9_3_2_binned_by_tether_length.csv")
    summary = _filter_converged(pd.read_csv(args.askite_summary))

    colors = {"2019": "C0", "2025": "C1", "all": "0.2"}
    markers = {"2019": "o", "2025": "s", "all": "o"}
    x_va_candidates = [
        "V_a_ms",
        "Va_ms",
        "va_ms",
        "V_a_bin_center_ms",
        "V_a_bin_center",
        "va_bin_center_ms",
        "bin_center",
    ]

    fig, axes = plt.subplots(2, 2, figsize=(7.5, 5.7), sharex=True)
    panels = [
        (
            axes[0, 0],
            [
                "tether_force_preferred_N_mean",
                "tether_force_kite_N_mean",
                "tether_force_N_mean",
                "force_N_mean",
                "tether_force_mean_N",
            ],
            _force_col(args.force_source),
            "Tether / aero force [N]",
        ),
        (
            axes[0, 1],
            ["C_L_mean", "CL_mean", "CL_ekf_mean"],
            "sim_CL_wing",
            "$C_L$",
        ),
        (
            axes[1, 0],
            ["C_D_mean", "CD_mean", "CD_ekf_mean"],
            "sim_CD_wing",
            "$C_D$",
        ),
        (
            axes[1, 1],
            ["L_over_D_mean", "CL_over_CD_mean", "glide_ratio_mean"],
            "sim_L_over_D_wing",
            "$C_L/C_D$",
        ),
    ]
    for campaign, df_campaign in _campaign_groups(binned_va, campaigns):
        color = colors.get(str(campaign), None)
        for ax, y_candidates, _, ylabel in panels:
            _plot_binned(
                ax,
                df_campaign,
                x_va_candidates,
                y_candidates,
                f"EKF {campaign}",
                color,
            )
            ax.set_ylabel(ylabel)
    for campaign, df_campaign in _campaign_groups(summary, campaigns):
        color = colors.get(str(campaign), None)
        marker = markers.get(str(campaign), "o")
        for ax, _, sim_col, _ in panels:
            _plot_sim(
                ax,
                df_campaign,
                "requested_V_a_ms",
                sim_col,
                f"ASKITE {campaign}",
                color,
                marker,
            )
    for ax in axes.flat:
        ax.grid(True, color="0.88", linewidth=0.6)
        handles, _ = ax.get_legend_handles_labels()
        if handles:
            ax.legend(frameon=False, fontsize=7)
    for ax in axes[1, :]:
        ax.set_xlabel("$V_a$ [m s$^{-1}$]")
    fig.tight_layout()
    save_figure(fig, args.output_dir, "fig_9_3_2_2_force_coefficients_vs_va", formats)

    fig, axes = plt.subplots(1, 3, figsize=(9.0, 3.1))
    for campaign, df_campaign in _campaign_groups(binned_tether, campaigns):
        color = colors.get(str(campaign), None)
        _plot_binned(
            axes[0],
            df_campaign,
            [
                "tether_length_m",
                "tether_length_bin_center_m",
                "tether_length_bin_center",
                "bin_center",
            ],
            ["L_over_D_mean", "CL_over_CD_mean", "glide_ratio_mean"],
            f"EKF {campaign}",
            color,
        )
    for campaign, df_campaign in _campaign_groups(summary, campaigns):
        color = colors.get(str(campaign), None)
        marker = markers.get(str(campaign), "o")
        _plot_sim(
            axes[0],
            df_campaign,
            "requested_tether_length_m",
            "sim_L_over_D_wing",
            f"ASKITE {campaign}",
            color,
            marker,
        )
        _plot_sim(
            axes[1],
            df_campaign,
            "requested_V_a_ms",
            "sim_L_over_D_wing",
            f"ASKITE {campaign}",
            color,
            marker,
        )
    axes[0].set_xlabel("Tether length [m]")
    axes[0].set_ylabel("$C_L/C_D$")
    axes[1].set_xlabel("$V_a$ [m s$^{-1}$]")
    axes[1].set_ylabel("$C_L/C_D$")
    if {
        "sim_L_over_D_wing",
        "measured_L_over_D_ekf",
        "requested_tether_length_m",
    }.issubset(summary.columns):
        axes[2].plot(
            pd.to_numeric(summary["requested_tether_length_m"], errors="coerce"),
            pd.to_numeric(summary["sim_L_over_D_wing"], errors="coerce")
            - pd.to_numeric(summary["measured_L_over_D_ekf"], errors="coerce"),
            "o",
            color="0.25",
        )
    axes[2].set_xlabel("Tether length [m]")
    axes[2].set_ylabel("ASKITE - EKF $C_L/C_D$")
    for ax in axes:
        ax.grid(True, color="0.88", linewidth=0.6)
        handles, _ = ax.get_legend_handles_labels()
        if handles:
            ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    save_figure(
        fig, args.output_dir, "fig_9_3_2_3_glide_ratio_vs_tether_length", formats
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

    if {
        "requested_V_a_ms",
        "sim_CL_wing",
        "measured_CL_ekf",
        "sim_CD_wing",
        "measured_CD_ekf",
    }.issubset(summary.columns):
        fig, ax = plt.subplots(figsize=(4.2, 3.1))
        x = pd.to_numeric(summary["requested_V_a_ms"], errors="coerce")
        ax.plot(
            x,
            pd.to_numeric(summary["sim_CL_wing"], errors="coerce")
            - pd.to_numeric(summary["measured_CL_ekf"], errors="coerce"),
            "o",
            label="$C_L$",
        )
        ax.plot(
            x,
            pd.to_numeric(summary["sim_CD_wing"], errors="coerce")
            - pd.to_numeric(summary["measured_CD_ekf"], errors="coerce"),
            "s",
            label="$C_D$",
        )
        ax.set_xlabel("$V_a$ [m s$^{-1}$]")
        ax.set_ylabel("ASKITE - EKF coefficient")
        ax.grid(True, color="0.88", linewidth=0.6)
        ax.legend(frameon=False)
        fig.tight_layout()
        save_figure(fig, args.output_dir, "sim_minus_ekf_CL_CD_vs_va", formats)

    if {
        "requested_tether_length_m",
        "sim_L_over_D_wing",
        "measured_L_over_D_ekf",
    }.issubset(summary.columns):
        fig, ax = plt.subplots(figsize=(4.2, 3.1))
        ax.plot(
            pd.to_numeric(summary["requested_tether_length_m"], errors="coerce"),
            pd.to_numeric(summary["sim_L_over_D_wing"], errors="coerce")
            - pd.to_numeric(summary["measured_L_over_D_ekf"], errors="coerce"),
            "o",
            color="0.25",
        )
        ax.set_xlabel("Tether length [m]")
        ax.set_ylabel("ASKITE - EKF $C_L/C_D$")
        ax.grid(True, color="0.88", linewidth=0.6)
        fig.tight_layout()
        save_figure(
            fig, args.output_dir, "sim_minus_ekf_glide_vs_tether_length", formats
        )


if __name__ == "__main__":
    main()
