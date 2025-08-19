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


# Import modules
def main():
    """Main function"""
    PROJECT_DIR = Path(__file__).resolve().parents[1]

    # load files
    results_dir = Path(PROJECT_DIR) / "results" / f"V3_25" / f"2025_06_03_1027h"
    results_dir = (
        Path(PROJECT_DIR) / "results" / f"TUDELFT_V3_KITE" / f"2025_08_15_1122h"
    )
    h5_path = Path(results_dir) / "sim_output.h5"
    meta_data_dict, tracking_data = load_sim_output(h5_path)

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

    # --- Interactive plot ---
    plotting.interactive_plot(
        tracking_data=tracking_data,
        kite_connectivity=np.array(meta_data_dict["kite_connectivity"]),
        rest_lengths=np.array(meta_data_dict["rest_lengths"]),
        f_ext=tracking_data["f_ext"],
        title="PSM Interactive",
        # elev=0,
        # azim=0,
        t_per_step=0.1,
    )


if __name__ == "__main__":
    main()
