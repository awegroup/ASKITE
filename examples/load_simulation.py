"""
### Info

Author: Jelle Poland \
Citing: https://doi.org/10.3390/en16145264 \
License: ... \
Github: ...
"""

from pathlib import Path
from kitesim.logging_config import *
from kitesim.utils import load_sim_output
from kitesim import plotting  # Add this import
import numpy as np


def _last_finite_value(values):
    """Return the last finite scalar in a tracking array, or NaN if absent."""
    if values is None:
        return np.nan
    values = np.asarray(values, dtype=float).reshape(-1)
    finite = values[np.isfinite(values)]
    return float(finite[-1]) if finite.size else np.nan


def _last_finite_vector(values):
    """Return the last all-finite 3-vector in a tracking array, or None."""
    if values is None:
        return None
    values = np.asarray(values, dtype=float)
    if values.ndim == 1 and values.size == 3 and np.all(np.isfinite(values)):
        return values
    values = values.reshape((-1, values.shape[-1]))
    finite_rows = values[np.all(np.isfinite(values), axis=1)]
    return finite_rows[-1] if finite_rows.size else None


def _metadata_value(meta_data_dict, *keys):
    """Return the first finite scalar metadata value found for the given keys."""
    for key in keys:
        if key not in meta_data_dict:
            continue
        try:
            value = float(np.asarray(meta_data_dict[key]).reshape(-1)[0])
        except (TypeError, ValueError, IndexError):
            continue
        if np.isfinite(value):
            return value
    return np.nan


def _print_final_aero_summary(meta_data_dict, tracking_data):
    """Print final tether force and aerodynamic coefficients from saved results."""
    aero_force_total = _last_finite_vector(tracking_data.get("aero_force_total"))
    if aero_force_total is None:
        aero_force_total = _last_finite_vector(
            tracking_data.get("aero_force_wing_total")
        )

    fixed_node_residual = None
    if "f_residual" in tracking_data:
        fixed_node_residual = np.asarray(
            tracking_data["f_residual"][-1, 0], dtype=float
        )
    elif "f_int" in tracking_data and "f_ext" in tracking_data:
        fixed_node_residual = np.asarray(
            tracking_data["f_int"][-1, 0] - tracking_data["f_ext"][-1, 0],
            dtype=float,
        )

    if aero_force_total is not None:
        tether_force = float(np.linalg.norm(aero_force_total))
        tether_force_source = "||aero_force_total||"
    elif fixed_node_residual is not None and np.all(np.isfinite(fixed_node_residual)):
        tether_force = float(np.linalg.norm(fixed_node_residual))
        tether_force_source = "||fixed node residual||"
    else:
        tether_force = _metadata_value(meta_data_dict, "tether_force")
        tether_force_source = "metadata tether_force"

    cl = _metadata_value(meta_data_dict, "final_C_L_wing", "cl")
    if not np.isfinite(cl):
        cl = _last_finite_value(tracking_data.get("C_L_wing"))

    cd = _metadata_value(meta_data_dict, "final_C_D_wing", "cd")
    if not np.isfinite(cd):
        cd = _last_finite_value(tracking_data.get("C_D_wing"))

    cl_cd = _metadata_value(meta_data_dict, "final_glide_ratio_wing")
    if (
        not np.isfinite(cl_cd)
        and np.isfinite(cl)
        and np.isfinite(cd)
        and abs(cd) > 1e-12
    ):
        cl_cd = cl / cd

    aoa_deg = _metadata_value(
        meta_data_dict,
        "aoa_deg",
        "final_center_or_midspan_alpha_deg",
    )

    print("\nFinal aerodynamic summary")
    print(f"  Tether force [N]: {tether_force:.3f} ({tether_force_source})")
    print(f"  Angle of attack [deg]: {aoa_deg:.3f}")
    print(f"  CL [-]: {cl:.6f}")
    print(f"  CD [-]: {cd:.6f}")
    print(f"  CL/CD [-]: {cl_cd:.6f}")


# Import modules
def main():
    """Main function"""
    PROJECT_DIR = Path(__file__).resolve().parents[1]

    # load files
    results_dir = (
        Path(PROJECT_DIR) / "results" / f"TUDELFT_V3_KITE" / f"2026_05_20_1707h"
    )
    # results_dir = (
    #     Path(PROJECT_DIR) / "results" / f"TUDELFT_V3_KITE" / f"2026_02_10_1128h"
    # )
    # results_dir = Path(PROJECT_DIR) / "results" / f"3plate_kite" / f"2025_10_16_1227h"

    h5_path = Path(results_dir) / "sim_output.h5"
    meta_data_dict, tracking_data = load_sim_output(h5_path)
    _print_final_aero_summary(meta_data_dict, tracking_data)

    logging.info(f"meta_data: {meta_data_dict}")
    print(f"tracking_data keys: {tracking_data.keys()}")
    print(f'final_node_positions: {len(tracking_data["positions"])}')
    print(f'shape: {tracking_data["positions"].shape}')
    print(f'final_node_positions: {tracking_data["positions"][-1]}')

    # TODO:
    # - here you could add functions to plot the tracking of f_int, f_ext and f_residual over the iterations
    # - functions that make an animation of the kite going through the iterations
    # - etc.
    f_residual = tracking_data["f_int"] - tracking_data["f_ext"]

    print(f'f_ext final: {tracking_data["f_ext"][-1]}')

    # --- Interactive plot ---
    plotting.interactive_plot(
        tracking_data=tracking_data,
        kite_connectivity_arr=np.array(meta_data_dict["kite_connectivity"]),
        rest_lengths=np.array(meta_data_dict["rest_lengths"]),
        f_ext=tracking_data["f_ext"],
        title="PSM Interactive",
        # elev=0,
        # azim=0,
        t_per_step=0.1,
    )

    _print_final_aero_summary(meta_data_dict, tracking_data)


if __name__ == "__main__":
    main()
