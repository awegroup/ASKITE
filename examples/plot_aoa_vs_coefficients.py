"""
Plot lift and drag coefficients as a function of angle of attack.

Reads a sweep CSV file and plots all data points.
X-axis: angle of attack [deg]
Y-axes: lift coefficient (Cl) and drag coefficient (Cd)

Usage:
    python examples/plot_aoa_vs_coefficients.py
    python examples/plot_aoa_vs_coefficients.py --csv <path_to_csv>
    python examples/plot_aoa_vs_coefficients.py --csv <path_to_csv> --save
    python examples/plot_aoa_vs_coefficients.py --csv results/TUDELFT_V3_KITE/sweep_course_steering_depower.csv --save

Examples:
    # Plot default sweep_wind_steering.csv
    python examples/plot_aoa_vs_coefficients.py

    # Plot a specific sweep CSV
    python examples/plot_aoa_vs_coefficients.py --csv results/TUDELFT_V3_KITE/sweep_course_steering_depower.csv

    # Save the plot
    python examples/plot_aoa_vs_coefficients.py --csv results/TUDELFT_V3_KITE/sweep_course_steering_depower.csv --save
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT_DIR = Path(__file__).resolve().parents[1]
DEFAULT_CSV = PROJECT_DIR / "results" / "TUDELFT_V3_KITE" / "sweep_c3_s9_d4.csv"


def load_csv(csv_path: Path) -> pd.DataFrame:
    """Load CSV with all data."""
    df = pd.read_csv(csv_path)
    return df


def plot_aoa_vs_coefficients(df: pd.DataFrame, save_dir: Path | None = None) -> None:
    """
    Create subplots: separate axes for Cl and Cd vs AoA.
    Colors by course angle to show overlapping points from different parameter combinations.

    Args:
        df: DataFrame with sweep results
        save_dir: Directory to save figure (if None, don't save)
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Get unique course angles for color coding
    if "angle_course_deg" in df.columns:
        courses = df["angle_course_deg"].unique()
        print(f"Course angles found in data: {sorted(courses)}")
        print(f"Data points per course:")
        for c in sorted(courses):
            count = len(df[df["angle_course_deg"] == c])
            print(f"  Course {int(c)}°: {count} points")

        colors = plt.cm.Set2(np.linspace(0, 1, len(courses)))
        course_map = {course: colors[i] for i, course in enumerate(courses)}
    else:
        course_map = None

    # Plot Cl vs AoA (left subplot)
    if course_map is not None:
        for course in sorted(courses):
            mask = df["angle_course_deg"] == course
            ax1.scatter(
                df[mask]["aoa_deg"],
                df[mask]["cl"],
                s=120,
                alpha=0.6,
                edgecolors="black",
                linewidth=0.5,
                label=f"Course {int(course)}°",
                color=course_map[course],
            )
    else:
        ax1.scatter(
            df["aoa_deg"],
            df["cl"],
            s=120,
            alpha=0.6,
            edgecolors="black",
            linewidth=0.5,
            label="Cl",
            color="steelblue",
        )

    ax1.set_xlabel("Angle of attack [deg]", fontsize=11, fontweight="bold")
    ax1.set_ylabel("Lift Coefficient (Cl)", fontsize=11, fontweight="bold")
    ax1.set_title(f"Cl vs AoA ({len(df)} data points)", fontsize=12, fontweight="bold")
    ax1.grid(True, alpha=0.3)
    ax1.axhline(0, color="k", linewidth=0.6, linestyle="--", alpha=0.4)
    ax1.legend(fontsize=9, loc="best")

    # Plot Cd vs AoA (right subplot)
    if course_map is not None:
        for course in sorted(courses):
            mask = df["angle_course_deg"] == course
            ax2.scatter(
                df[mask]["aoa_deg"],
                df[mask]["cd"],
                s=120,
                alpha=0.6,
                edgecolors="black",
                linewidth=0.5,
                marker="s",
                label=f"Course {int(course)}°",
                color=course_map[course],
            )
    else:
        ax2.scatter(
            df["aoa_deg"],
            df["cd"],
            s=120,
            alpha=0.6,
            edgecolors="black",
            linewidth=0.5,
            label="Cd",
            color="coral",
            marker="s",
        )

    ax2.set_xlabel("Angle of attack [deg]", fontsize=11, fontweight="bold")
    ax2.set_ylabel("Drag Coefficient (Cd)", fontsize=11, fontweight="bold")
    ax2.set_title(f"Cd vs AoA ({len(df)} data points)", fontsize=12, fontweight="bold")
    ax2.grid(True, alpha=0.3)
    ax2.axhline(0, color="k", linewidth=0.6, linestyle="--", alpha=0.4)
    ax2.legend(fontsize=9, loc="best")

    fig.suptitle(
        "Aerodynamic Coefficients vs Angle of Attack", fontsize=13, fontweight="bold"
    )
    fig.tight_layout()

    if save_dir is not None:
        save_dir.mkdir(parents=True, exist_ok=True)
        out = save_dir / "aoa_vs_coefficients.png"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"Saved: {out}")

    plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Plot lift and drag coefficients vs angle of attack from sweep CSV."
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=DEFAULT_CSV,
        help="Path to the sweep CSV file. Default: results/TUDELFT_V3_KITE/sweep_wind_steering.csv",
    )
    parser.add_argument(
        "--save",
        action="store_true",
        help="Save figure in same directory as CSV file.",
    )
    args = parser.parse_args()

    csv_path = Path(args.csv)

    # Handle relative and absolute paths
    if not csv_path.is_absolute():
        PROJECT_DIR = Path(__file__).resolve().parents[1]
        csv_path = PROJECT_DIR / csv_path

    if not csv_path.exists():
        print(f"Error: CSV file not found: {csv_path}")
        print(f"       (resolved to: {csv_path.resolve()})")
        return

    print(f"\nLoading CSV: {csv_path}")
    df = load_csv(csv_path)
    print(f"Loaded {len(df)} data points")

    # Check if required columns exist
    if "aoa_deg" not in df.columns:
        print("Error: 'aoa_deg' column not found in CSV")
        return

    if "cl" not in df.columns or "cd" not in df.columns:
        print("Error: 'cl' or 'cd' columns not found in CSV")
        print(f"Available columns: {list(df.columns)}")
        return

    # Count valid values
    n_aoa = df["aoa_deg"].notna().sum()
    n_cl = df["cl"].notna().sum()
    n_cd = df["cd"].notna().sum()

    print(f"Data summary:")
    print(
        f"  AoA:       {n_aoa:4d} non-NaN values (range: {df['aoa_deg'].min():.1f}° to {df['aoa_deg'].max():.1f}°)"
    )
    print(
        f"  Cl:        {n_cl:4d} non-NaN values (range: {df['cl'].min():.3f} to {df['cl'].max():.3f})"
    )
    print(
        f"  Cd:        {n_cd:4d} non-NaN values (range: {df['cd'].min():.3f} to {df['cd'].max():.3f})"
    )

    if n_cl == 0 or n_cd == 0:
        print("\nWarning: No valid Cl or Cd data found.")
        print("These columns may not be saved from the solver yet.")
        print("Check if aerostructural_coupled_solver_qsm is populating 'cl' and 'cd'.")
        return

    print(f"\nPlotting {n_cl} Cl points and {n_cd} Cd points...\n")

    save_dir = csv_path.parent if args.save else None

    plot_aoa_vs_coefficients(df, save_dir=save_dir)


if __name__ == "__main__":
    main()
