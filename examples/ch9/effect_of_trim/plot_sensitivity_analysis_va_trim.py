import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

CH9_DIR = Path(__file__).resolve().parents[1]
if str(CH9_DIR) not in sys.path:
    sys.path.insert(0, str(CH9_DIR))

from ch9_analysis_utils import PROJECT_DIR, parse_formats, save_figure, write_markdown_table


def build_parser():
    parser = argparse.ArgumentParser(
        description="Plot Section 9.3.1 V_a and depower sensitivity summaries."
    )
    parser.add_argument(
        "--summary-csv",
        type=Path,
        default=PROJECT_DIR
        / "results"
        / "ch9"
        / "effect_of_trim"
        / "processed_data"
        / "sensitivity_va_trim_summary.csv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_DIR / "results" / "ch9" / "effect_of_trim",
    )
    parser.add_argument("--plot-format", default="pdf")
    parser.add_argument("--include-failed", action="store_true")
    parser.add_argument("--case-root", type=Path, default=None)
    return parser


def _numeric(df, column):
    return pd.to_numeric(df[column], errors="coerce")


def _plot_family(df, x_col, x_label, stem, output_dir, formats):
    if df.empty:
        return []
    df = df.sort_values(x_col)
    x = _numeric(df, x_col)
    fig, axes = plt.subplots(2, 2, figsize=(7.5, 5.6), sharex=True)

    ax = axes[0, 0]
    ax.plot(x, _numeric(df, "C_L_wing"), "o-", label="$C_L$")
    ax.plot(x, _numeric(df, "C_D_wing"), "s-", label="$C_D$")
    ax.plot(x, _numeric(df, "glide_ratio_wing"), "^-", label="$C_L/C_D$")
    ax.set_ylabel("Coefficient")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[0, 1]
    ax.plot(x, _numeric(df, "resultant_force_N"), "o-", label="Aero total")
    if "tether_or_reaction_force_N" in df:
        ax.plot(x, _numeric(df, "tether_or_reaction_force_N"), "s-", label="Reaction")
    ax.set_ylabel("Force [N]")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 0]
    ax.plot(x, _numeric(df, "projected_span_m"), "o-", label="Span")
    ax.plot(x, _numeric(df, "projected_area_m2"), "s-", label="Area")
    ax.set_ylabel("Geometry [m, m$^2$]")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 1]
    ax.plot(x, _numeric(df, "midspan_pitch_deg"), "o-", label="Midspan pitch")
    ax.plot(x, _numeric(df, "mean_twist_deg"), "s-", label="Mean twist")
    ax.plot(x, _numeric(df, "max_abs_twist_deg"), "^-", label="Max |twist|")
    ax.set_ylabel("Angle [deg]")
    ax.legend(frameon=False, fontsize=8)

    for ax in axes.flat:
        ax.grid(True, color="0.88", linewidth=0.6)
    for ax in axes[1, :]:
        ax.set_xlabel(x_label)
    fig.tight_layout()
    return save_figure(fig, output_dir, stem, formats)


def _plot_metric_grid(df, x_col, x_label, stem, output_dir, formats):
    if df.empty:
        return []
    df = df.sort_values(x_col)
    x = _numeric(df, x_col)
    force_col = (
        "tether_or_reaction_force_N"
        if "tether_or_reaction_force_N" in df.columns
        else "resultant_force_N"
    )
    panels = [
        ("C_L_wing", "$C_L$", "[-]"),
        ("C_D_wing", "$C_D$", "[-]"),
        ("glide_ratio_wing", "$C_L/C_D$", "[-]"),
        (force_col, "Reaction force", "[N]"),
        ("projected_span_m", "Projected span", "[m]"),
        ("projected_area_m2", "Projected area", "[m$^2$]"),
        ("midspan_pitch_deg", "Midspan pitch", "[deg]"),
        ("mean_twist_deg", "Mean twist", "[deg]"),
        ("max_abs_twist_deg", "Max |twist|", "[deg]"),
    ]
    fig, axes = plt.subplots(3, 3, figsize=(8.2, 6.8), sharex=True)
    for ax, (column, title, units) in zip(axes.flat, panels):
        y = _numeric(df, column)
        ax.plot(x, y, "o-", color="0.15", linewidth=1.0, markersize=3.5)
        ax.set_title(title, fontsize=9)
        ax.set_ylabel(units)
        ax.grid(True, color="0.88", linewidth=0.6)
    for ax in axes[-1, :]:
        ax.set_xlabel(x_label)
    fig.tight_layout()
    return save_figure(fig, output_dir, stem, formats)


def _table_rows(df, x_col):
    cols = [
        x_col,
        "converged",
        "C_L_wing",
        "C_D_wing",
        "glide_ratio_wing",
        "resultant_force_N",
        "projected_span_m",
        "mean_twist_deg",
    ]
    cols = [col for col in cols if col in df.columns]
    return df.sort_values(x_col)[cols].to_dict("records")


def main():
    args = build_parser().parse_args()
    print("creating plots")
    df = pd.read_csv(args.summary_csv)
    if not args.include_failed and "converged" in df.columns:
        df = df[df["converged"].astype(str).str.lower().isin(["true", "1"])]
    formats = parse_formats(args.plot_format)

    va_df = df[df["sweep_type"] == "va"] if "sweep_type" in df else pd.DataFrame()
    depower_df = (
        df[df["sweep_type"] == "depower"] if "sweep_type" in df else pd.DataFrame()
    )
    created_files = []
    created_files.extend(
        _plot_family(
            va_df,
            "requested_va_ms",
            "$V_a$ [m s$^{-1}$]",
            "fig_9_3_1_5_effect_of_va",
            args.output_dir,
            formats,
        )
    )
    created_files.extend(
        _plot_family(
            depower_df,
            "requested_u_dp",
            "$u_{dp}$ [-]",
            "fig_9_3_1_6_effect_of_trim_depower",
            args.output_dir,
            formats,
        )
    )
    created_files.extend(
        _plot_metric_grid(
            va_df,
            "requested_va_ms",
            "$V_a$ [m s$^{-1}$]",
            "fig_9_3_1_5_effect_of_va_3x3",
            args.output_dir,
            formats,
        )
    )
    created_files.extend(
        _plot_metric_grid(
            depower_df,
            "requested_u_dp",
            "$u_{dp}$ [-]",
            "fig_9_3_1_6_effect_of_trim_depower_3x3",
            args.output_dir,
            formats,
        )
    )
    if not va_df.empty:
        table_path = args.output_dir / "table_9_3_1_5_va_sweep.md"
        write_markdown_table(
            table_path,
            _table_rows(va_df, "requested_va_ms"),
        )
        created_files.append(table_path)
    if not depower_df.empty:
        table_path = args.output_dir / "table_9_3_1_6_trim_sweep.md"
        write_markdown_table(
            table_path,
            _table_rows(depower_df, "requested_u_dp"),
        )
        created_files.append(table_path)

    created_files = sorted(created_files)
    if created_files:
        print("created files:")
        for path in created_files:
            print(f"- {path}")
    else:
        print(f"no files were created in {args.output_dir}")


if __name__ == "__main__":
    main()
