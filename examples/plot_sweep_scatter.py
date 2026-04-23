"""
Scatter plot of sweep_wind_steering data.

X-axis: kite_speed × steering input
Y-axis: course rate
Color: variable (steering input by default)

Usage:
    python examples/plot_sweep_scatter.py
    python examples/plot_sweep_scatter.py --csv results/TUDELFT_V3_KITE/sweep_wind_steering.csv
    python examples/plot_sweep_scatter.py --csv results/TUDELFT_V3_KITE/sweep_wind_steering.csv --color wind_speed_wind_ref_ms
    python examples/plot_sweep_scatter.py --csv results/TUDELFT_V3_KITE/sweep_wind_steering.csv --color opt_kite_speed --save
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT_DIR = Path(__file__).resolve().parents[1]
DEFAULT_CSV = PROJECT_DIR / "results" / "TUDELFT_V3_KITE" / "sweep_c1_s11_d1.csv"

# Case keys for deduplication
CASE_KEYS = [
    "wind_speed_wind_ref_ms",
    "depower_tape_final_extension_m",
    "steering_tape_final_extension_m",
]

# Available color variables
COLOR_VARS = {
    "steering_tape_final_extension_m": "Steering input [m]",
    "wind_speed_wind_ref_ms": "Wind speed [m/s]",
    "opt_kite_speed": "Kite speed [m/s]",
    "depower_tape_final_extension_m": "Depower input [m]",
    "aoa_deg": "Angle of attack [deg]",
}


def load_csv(csv_path: Path) -> pd.DataFrame:
    """Load CSV with all data (no deduplication)."""
    df = pd.read_csv(csv_path)
    return df


def plot_scatter(
    df: pd.DataFrame,
    color_var: str = "steering_tape_final_extension_m",
    save_dir: Path | None = None,
) -> None:
    """
    Create scatter plot: X = speed * steering, Y = course rate, colored by color_var.

    Args:
        df: DataFrame with sweep results
        color_var: Column name to use for coloring (default: steering)
        save_dir: Directory to save figure (if None, don't save)
    """
    if color_var not in COLOR_VARS and color_var not in df.columns:
        raise ValueError(f"Color variable '{color_var}' not found in CSV")

    # Compute axes
    x = df["va"] * df["steering_tape_final_extension_m"]
    y = df["opt_course_rate_body"]
    c = df[color_var]

    # Create figure
    fig, ax = plt.subplots(figsize=(11, 8))

    # Scatter plot
    scatter = ax.scatter(
        x,
        y,
        c=c,
        s=120,
        cmap="viridis",
        alpha=0.7,
        edgecolors="black",
        linewidth=0.5,
    )

    # Labels and title
    ax.set_xlabel("Kite speed × Steering input [m²/s]", fontsize=12, fontweight="bold")
    ax.set_ylabel("Course rate [rad/s]", fontsize=12, fontweight="bold")

    color_label = COLOR_VARS.get(color_var, color_var)
    ax.set_title(
        f"Steering sweep: Speed × Steering vs Course Rate\n(colored by {color_label})",
        fontsize=13,
        fontweight="bold",
    )

    ax.grid(True, alpha=0.3)
    ax.axhline(0, color="k", linewidth=0.6, linestyle="--", alpha=0.4)

    # Colorbar
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label(color_label, fontsize=11, fontweight="bold")

    fig.tight_layout()

    if save_dir is not None:
        save_dir.mkdir(parents=True, exist_ok=True)
        out = save_dir / f"sweep_scatter_{color_var}.png"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"Saved: {out}")

    plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Scatter plot of steering sweep results (speed × steering vs course rate)."
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=DEFAULT_CSV,
        help="Path to the sweep CSV file.",
    )
    parser.add_argument(
        "--color",
        type=str,
        default="steering_tape_final_extension_m",
        help=f"Column to use for coloring. Options: {', '.join(COLOR_VARS.keys())}",
    )
    parser.add_argument(
        "--save",
        action="store_true",
        help="Save figure next to the CSV file.",
    )
    args = parser.parse_args()

    csv_path = args.csv
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV not found: {csv_path}")

    df = load_csv(csv_path)
    print(f"Loaded {len(df)} data points from {csv_path}")
    print(f"Coloring by: {COLOR_VARS.get(args.color, args.color)}")

    save_dir = csv_path.parent if args.save else None

    plot_scatter(df, color_var=args.color, save_dir=save_dir)


if __name__ == "__main__":
    main()
