"""
### Info

Author: Jelle Poland \
Citing: https://doi.org/10.3390/en16145264 \
License: ... \
Github: ...
"""

from pathlib import Path
from kitesim.logging_config import *
from kitesim.utils import load_and_save_config_files, load_sim_output, save_results

from kitesim import (
    aerodynamic,
    initialisation,
    struc2aero,
    aero2struc,
    structural,
    tracking,
    plotting,
    aerostructural_coupled_solver,
)


# Import modules
def main():
    """Main function"""
    PROJECT_DIR = Path(__file__).resolve().parents[1]
    kite_name = "TUDELFT_V3_KITE"  # the dir name with the relevant .yaml files

    # load config.yaml & geometry.yaml, save both, and return them as dicts
    config, geometry, results_dir = load_and_save_config_files(kite_name, PROJECT_DIR)
    logging.info(f"config files saved in {results_dir}\n")

    ## AERO
    n_struc_ribs = len(struc_geometry["wing_nodes"]["data"]) / 2
    n_panels_aero = (n_struc_ribs - 1) * config["aerodynamic"][
        "n_aero_panels_per_struc_section"
    ]
    body_aero, vsm_solver, vel_app, initial_polar_data = aerodynamic.initialize(
        kite_name,
        PROJECT_DIR,
        config,
        n_panels_aero,
    )

    ## GENERAL
    (
        struc_nodes,
        wing_ci,
        wing_cj,
        bridle_ci,
        bridle_cj,
        struc_node_le_indices,  # was le_node_indices
        struc_node_te_indices,  # was te_node_indices
        pulley_point_indices,
        tubular_frame_line_idx_list,
        te_line_idx_list,
        n_struc_ribs,
        wing_connectivity,
        bridle_connectivity,
        kite_connectivity,
        m_array,
        bridle_rest_lengths_initial,
        wing_rest_lengths_initial,
        rest_lengths,
        pulley_point_indices,
        pulley_line_indices,
        pulley_line_to_other_node_pair_dict,
        power_tape_index,
        steering_tape_indices,
    ) = initialisation.main(geometry)

    ## STRUC

    pulley_line_to_other_node_pair_dict
    rest_lengths
    struc_nodes
    kite_connectivity
    m_array
    power_tape_index --> defined in struc_geometry
    steering_tape_index --> defined in struc_geometry
    fixed_nodes_indices --> defined in struc_geometry


    psystem, params, pss_kite_connectivity = structural.instantiate_psystem(
        config,
        geometry,
        struc_nodes,
        wing_connectivity,
        kite_connectivity,
        rest_lengths,
        m_array,
        tubular_frame_line_idx_list,
        te_line_idx_list,
        pulley_point_indices,
        pulley_line_indices,
        pulley_line_to_other_node_pair_dict,
        power_tape_index,
    )

    ## AERO2STRUC
    aero2struc_mapping = aero2struc.initialize_mapping(
        body_aero.panels,
        struc_nodes,
        struc_node_le_indices,
        struc_node_te_indices,
    )

    ## ACTUATION
    initial_length_power_tape = params["l0"][power_tape_index]
    power_tape_extension_step = config["power_tape_extension_step"]
    power_tape_final_extension = config["power_tape_final_extension"]
    if power_tape_extension_step != 0:
        n_power_tape_steps = int(power_tape_final_extension / power_tape_extension_step)
    else:
        n_power_tape_steps = 0
    logging.info(
        f"Initial depower tape length: {psystem.extract_rest_length[power_tape_index]:.3f}m"
    )
    logging.info(
        f"Desired depower tape length: {initial_length_power_tape + power_tape_final_extension:.3f}m"
    )
    initial_length_steering_left = params["l0"][steering_tape_indices[0]]
    initial_length_steering_right = params["l0"][steering_tape_indices[1]]
    steering_tape_extension_step = config["steering_tape_extension_step"]
    steering_tape_final_extension = config["steering_tape_final_extension"]
    if steering_tape_extension_step != 0:
        n_steering_tape_steps = int(
            steering_tape_final_extension / steering_tape_extension_step
        )
    else:
        n_steering_tape_steps = 0

    psystem.update_rest_length(steering_tape_indices[0], -steering_tape_final_extension)
    psystem.update_rest_length(steering_tape_indices[1], steering_tape_final_extension)

    ## Run AeroStructural simulation
    tracking_data, meta = aerostructural_coupled_solver.main(
        m_array,
        struc_nodes,
        psystem,
        struc_node_le_indices,
        struc_node_te_indices,
        body_aero,
        vsm_solver,
        vel_app,
        initial_polar_data,
        aero2struc_mapping,
        pss_kite_connectivity,
        params,
        power_tape_index,
        initial_length_power_tape,
        n_power_tape_steps,
        power_tape_final_extension,
        power_tape_extension_step,
        config,
    )
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
