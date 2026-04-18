"""
Plot statistics from a steering/wind sweep CSV produced by run_sweep_wind_steering.py.

Creates separate figures grouped by theme:
  - Aerodynamics: AoA, side-slip, Cl, Cd
  - Kite attitude: roll, pitch, yaw, course rate
  - Kite speed: velocity magnitude

Non-converged points are marked with 'x' symbols.

Usage:
    python examples/plot_sweep_wind_steering.py
    python examples/plot_sweep_wind_steering.py --csv results/TUDELFT_V3_KITE/sweep_wind_steering.csv
    python examples/plot_sweep_wind_steering.py --csv results/TUDELFT_V3_KITE/sweep_wind_steering.csv --save
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd

PROJECT_DIR = Path(__file__).resolve().parents[1]
DEFAULT_CSV = PROJECT_DIR / "results" / "TUDELFT_V3_KITE" / "sweep_wind_steering.csv"

# Columns that uniquely identify one simulation run.
CASE_KEYS = [
    "wind_speed_wind_ref_ms",
    "depower_tape_final_extension_m",
    "steering_tape_final_extension_m",
]

# Quantities to plot against steering extension.
PLOT_VARS = [
    ("opt_kite_speed", "Kite speed [m/s]"),
    ("opt_course_rate_body", "Course rate [rad/s]"),
]


def load_and_deduplicate(csv_path: Path) -> pd.DataFrame:
    """
    Load CSV and keep only the last row for each unique (wind, depower, steering)
    combination — handles repeated appends from multiple sweep runs.
    """
    df = pd.read_csv(csv_path)
    df = df.drop_duplicates(subset=CASE_KEYS, keep="last").reset_index(drop=True)
    df = df.sort_values(CASE_KEYS).reset_index(drop=True)
    return df


def make_label(wind: float, depower: float) -> str:
    return f"V_w={wind:.1f} m/s  d_p={depower*1000:.0f} mm"


def plot_sweep(df: pd.DataFrame, save_dir: Path | None = None) -> None:
    groups = df.groupby(["wind_speed_wind_ref_ms", "depower_tape_final_extension_m"])
    n_groups = len(groups)

    colors = cm.tab10(np.linspace(0, 0.9, n_groups))

    # Create single figure with subplots for speed and course_rate
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    for ax, (col, ylabel) in zip(axes, PLOT_VARS):
        for color, ((wind, depower), gdf) in zip(colors, groups):
            # Skip if column doesn't exist
            if col not in gdf.columns:
                continue
            x = gdf["steering_tape_final_extension_m"] * 1000  # mm
            y = gdf[col]
            label = make_label(wind, depower)
            ax.plot(
                x,
                y,
                marker="o",
                markersize=5,
                linewidth=2,
                color=color,
                label=label,
            )
        ax.set_ylabel(ylabel, fontsize=11, fontweight="bold")
        ax.grid(True, alpha=0.3)
        ax.axhline(0, color="k", linewidth=0.6, linestyle="--", alpha=0.4)

    axes[-1].set_xlabel("Steering tape extension [mm]", fontsize=11, fontweight="bold")

    # Single legend below
    handles, labels = axes[0].get_legend_handles_labels()

    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=min(n_groups, 3),
        fontsize=9,
        bbox_to_anchor=(0.5, -0.05),
        frameon=True,
    )

    fig.suptitle(
        "Steering sweep — Speed & Course Rate (QSM)",
        fontsize=13,
        fontweight="bold",
        y=0.995,
    )
    fig.tight_layout()

    if save_dir is not None:
        save_dir.mkdir(parents=True, exist_ok=True)
        out = save_dir / "sweep_speed_course_rate.png"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"Saved: {out}")

    plt.show()


def plot_scatter_speed_vs_course_rate(
    df: pd.DataFrame, save_dir: Path | None = None
) -> None:
    """Scatter plot: X=speed*steering, Y=-course_rate."""
    fig, ax = plt.subplots(figsize=(10, 7))

    # Compute x-axis: speed * steering deflection
    x = df["opt_kite_speed"] * df["steering_tape_final_extension_m"]
    y = -df["opt_course_rate_body"]  # Negate course rate

    # Create scatter plot
    ax.scatter(
        x,
        y,
        s=100,
        c="steelblue",
        alpha=0.7,
        edgecolors="black",
        linewidth=0.5,
    )

    ax.set_xlabel("Speed × Steering deflection [m²/s]", fontsize=12, fontweight="bold")
    ax.set_ylabel("Course rate [rad/s]", fontsize=12, fontweight="bold")
    ax.set_title("Speed × Steering vs Course Rate", fontsize=13, fontweight="bold")
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color="k", linewidth=0.6, linestyle="--", alpha=0.4)

    fig.tight_layout()

    if save_dir is not None:
        save_dir.mkdir(parents=True, exist_ok=True)
        out = save_dir / "sweep_scatter_speed_course_rate.png"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"Saved: {out}")

    plt.show()


def plot_convergence_summary(df: pd.DataFrame, save_dir: Path | None = None) -> None:
    """Bar chart showing how many cases converged per (wind, depower) group."""
    groups = df.groupby(["wind_speed_wind_ref_ms", "depower_tape_final_extension_m"])
    labels, conv_counts, total_counts = [], [], []

    for (wind, depower), gdf in groups:
        labels.append(make_label(wind, depower))
        total_counts.append(len(gdf))
        conv_counts.append(int(gdf["converged"].astype(bool).sum()))

    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(max(6, len(labels) * 1.6), 3.5))
    ax.bar(x, total_counts, color="lightsteelblue", label="total")
    ax.bar(x, conv_counts, color="steelblue", label="converged")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=15, ha="right", fontsize=8)
    ax.set_ylabel("Number of cases")
    ax.set_title("Convergence summary")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()

    if save_dir is not None:
        out = save_dir / "sweep_convergence_summary.png"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"Saved: {out}")

    plt.show()


def main():
    parser = argparse.ArgumentParser(description="Plot steering sweep CSV statistics.")
    parser.add_argument(
        "--csv",
        type=Path,
        default=DEFAULT_CSV,
        help="Path to the sweep CSV file.",
    )
    parser.add_argument(
        "--save",
        action="store_true",
        help="Save figures next to the CSV file.",
    )
    args = parser.parse_args()

    csv_path = args.csv
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV not found: {csv_path}")

    df = load_and_deduplicate(csv_path)
    print(f"Loaded {len(df)} unique cases from {csv_path}")
    print(df[CASE_KEYS].to_string(index=False))

    save_dir = csv_path.parent if args.save else None

    plot_sweep(df, save_dir=save_dir)
    plot_scatter_speed_vs_course_rate(df, save_dir=save_dir)


if __name__ == "__main__":
    main()
