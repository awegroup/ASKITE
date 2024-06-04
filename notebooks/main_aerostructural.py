"""
### Info

Author: Jelle Poland \
Citing: https://doi.org/10.3390/en16145264 \
License: ... \
Github: ...
"""

### Initialisation

# Making things autoreload - needed for Jupyter Kernel/Interactive env.
# %load_ext autoreload
# %autoreload 2
# %matplotlib widget

import os
import sys
import matplotlib.pyplot as plt

# from IPython.display import display, Latex

# TODO: can we remove this?!?
# Define the right path
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(f"{project_root}")  # Needed for running in terminal
sys.path.insert(0, f"{project_root}")  # Needed for running in terminal
os.chdir(f"{project_root}")  # Needed for running in interactive python environment

from kitesim import parent_adapter


# Import modules
def main():
    """Main function"""

    kite_name = "V3_25"

    # loading immutable variables
    config_path = "data/config.yaml"
    # TODO: Remove the output_path from config?
    folder_output_path = f"results/{kite_name}"
    folder_path_kite = f"processed_data/{kite_name}"
    folder_path_kite_data = f"processed_data/{kite_name}/processed_design_files"
    case_path_folder = "../kitesim/src/kitesim/cases"
    folder_name_results = "results/"

    # Importing modules
    (
        InputVSM,
        InputBridleAero,
        setup_config,
        get_mutable_variables,
        solver_main,
        post_processing_main,
    ) = parent_adapter.module_importer()

    config = setup_config(
        config_path,
        folder_output_path,
        folder_path_kite,
        folder_path_kite_data,
        case_path_folder,
    )
    input_VSM = InputVSM.create(config)
    input_bridle_aero = InputBridleAero.create(config)

    # Get mutable variables
    points_ini, vel_app, params_dict, psystem = get_mutable_variables(config)

    # AeroStructural Simulation
    points, df_position, post_processing_data = solver_main.run_aerostructural_solver(
        points_ini,
        vel_app,
        psystem,
        params_dict,
        config,
        input_VSM,
        input_bridle_aero,
    )

    # TODO: Should this be placed inside the solver_main loop?
    # Saving non-interpretable results
    folder_name = post_processing_main.save_non_interpretable_results(
        config,
        input_VSM,
        input_bridle_aero,
        points_ini,
        vel_app,
        params_dict,
        psystem,
        points,
        df_position,
        post_processing_data,
        folder_name_results,
    )
    # Saving interpretable results (generated from the saved non-interpretable results)
    # post_processing_main.save_interpretable_results(folder_name)

    if config.is_with_printing:
        post_processing_main.print_results(
            points,
            post_processing_data["print_data"],
            config,
        )
    if config.is_with_plotting:
        post_processing_main.plot(
            post_processing_data["plot_data"],
            points,
            vel_app,
            config,
        )
    plt.show()
    if config.is_with_animation:
        print(f"")
        print("--> Generating ANIMATION \{*_*}/")
        post_processing_main.animate(
            post_processing_data["animation_data"],
            vel_app,
            config,
            input_VSM,
        )

    # TODO: this should not longer be needed?
    # Post-processing output
    # np.save(
    #     f"data/output/{config.kite_name}/points/{config.sim_name}_power_{1e3*config.depower_tape_final_extension:.0f}_steer_{1e3*np.abs(config.steering_tape_final_extension):.0f}.npy",
    #     points,
    # )


if __name__ == "__main__":
    main()
