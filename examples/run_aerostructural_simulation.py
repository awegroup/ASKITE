"""
### Info

Author: Jelle Poland \
Citing: https://doi.org/10.3390/en16145264 \
License: ... \
Github: ...
"""

from pathlib import Path
from kitesim import solver
from kitesim.logging_config import *
from kitesim.utils import load_and_save_config_files, load_sim_output, save_results


# Import modules
def main():
    """Main function"""
    PROJECT_DIR = Path(__file__).resolve().parents[1]

    # load and save config files
    config, config_kite, results_dir = load_and_save_config_files(PROJECT_DIR)
    logging.info(f"config files saved in {results_dir}\n")

    # run AeroStructural simulation
    tracking_data, meta = solver.main(config, config_kite)
    h5_path = Path(results_dir) / "sim_output.h5"
    save_results(tracking_data, meta, h5_path)

    # Load results
    meta_data_dict, tracking_data = load_sim_output(h5_path)

    logging.info(f"meta_data: {meta_data_dict}")

    # TODO:
    # - here you could add functions to plot the tracking of f_int, f_ext and f_residual over the iterations
    # - functions that make an animation of the kite going through the iterations
    # - etc.
    f_residual = tracking_data["f_int"] - tracking_data["f_ext"]


if __name__ == "__main__":
    main()
