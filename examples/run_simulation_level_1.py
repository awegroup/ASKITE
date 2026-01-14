"""
### Info

Author: Jelle Poland \
Citing: https://doi.org/10.3390/en16145264 \
License: ... \
Github: ...
"""

import numpy as np
from pathlib import Path
from kitesim.logging_config import *
from datetime import datetime
from kitesim.utils import (
    load_and_save_config_files,
    load_yaml,
    load_sim_output,
    save_results,
    printing_rest_lengths,
)
from kitesim import (
    aero2struc,
    aerodynamic_vsm,
    structural_kite_fem,
    structural_pss,
    aerostructural_coupled_solver,
    read_struc_geometry_yaml,
)


# Import modules
def main():
    """Main function"""
    PROJECT_DIR = Path(__file__).resolve().parents[1]
    kite_name = "TUDELFT_V3_KITE"  # the dir name with the relevant .yaml files
    # kite_name = "3plate_kite"  # the dir name with the relevant .yaml files
    # load config.yaml & geometry.yaml, save both, and return them as dicts
    config_path = Path(PROJECT_DIR) / "data" / f"{kite_name}" / "config.yaml"
    struc_geometry_path = (
        Path(PROJECT_DIR)
        / "data"
        / f"{kite_name}"
        / "struc_geometry_level_1_manual.yaml"
    )
    aero_geometry_path = (
        Path(PROJECT_DIR) / "data" / f"{kite_name}" / "aero_geometry.yaml"
    )

    results_dir = (
        Path(PROJECT_DIR)
        / "results"
        / f"{kite_name}"
        / f'{datetime.now().strftime("%Y_%m_%d_%H%M")}h'
    )
    config, struc_geometry, aero_geometry, results_dir = load_and_save_config_files(
        config_path, struc_geometry_path, aero_geometry_path, results_dir
    )
    logging.info(f"config files saved in {results_dir}\n")

    ###################
    ### AERODYNAMIC ###
    ###################
    n_wing_struc_nodes = len(struc_geometry["wing_particles"]["data"])
    n_struc_ribs = n_wing_struc_nodes / 2
    n_panels_aero = (n_struc_ribs - 1) * config["aerodynamic"][
        "n_aero_panels_per_struc_section"
    ]
    body_aero, vsm_solver, vel_app, initial_polar_data = aerodynamic_vsm.initialize(
        kite_name,
        PROJECT_DIR,
        config,
        n_panels_aero,
    )

    ##################
    ### STRUCTURAL ###
    ##################
    (
        # node level
        struc_nodes,
        m_arr,
        struc_node_le_indices,
        struc_node_te_indices,
        power_tape_index,
        steering_tape_indices,
        pulley_node_indices,
        # element level
        kite_connectivity_arr,
        bridle_connectivity_arr,
        bridle_diameter_arr,
        l0_arr,
        k_arr,
        c_arr,
        linktype_arr,
        pulley_line_indices,
        pulley_line_to_other_node_pair_dict,
    ) = read_struc_geometry_yaml.main(struc_geometry)

    # logging initial conditions
    logging.info(f"\n\nINITIAL CONDITIONS, NODES \n")
    for idx, (node_i, m_i) in enumerate(zip(struc_nodes, m_arr)):
        logging.info(f"node_idx: {idx}: node: {node_i}, mass: {m_i}")

    logging.info(f"\n\nINITIAL CONDITIONS, ELEMENTS \n")
    for idx, conn in enumerate(kite_connectivity_arr):
        logging.info(
            f"conn_idx: {idx}: conn: {conn}, l0: {l0_arr[idx]}, k: {k_arr[idx]}, c: {c_arr[idx]}, linktype: {linktype_arr[idx]}"
        )

    if config["structural_solver"] == "pss":
        ## pss -- https://github.com/awegroup/Particle_System_Simulator
        ##TODO: Fix the comment below, it SHOULD read l0
        # Note: ParticleSystem doesn’t read l0_arr. SpringDamper sets l0
        # from the initial particle positions.
        # So l0_arr is a bookkeeping array for you, not used at instantiation.
        (psystem, pss_initial_conditions, pss_params, struc_nodes_initial) = (
            structural_pss.instantiate(
                # yaml files
                config,
                # node level
                struc_nodes,
                m_arr,
                # element_level
                kite_connectivity_arr,
                l0_arr,
                k_arr,
                c_arr,
                linktype_arr,
                pulley_line_to_other_node_pair_dict,
            )
        )
        if config["is_with_initial_structure_plot"]:
            structural_pss.plot_3d_kite_structure(
                struc_nodes,
                kite_connectivity_arr,
                power_tape_index,
                k_arr=k_arr,
                c_arr=c_arr,
                linktype_arr=linktype_arr,
                pulley_nodes=pulley_node_indices,
            )
        # setting kite_fem related output to None
        kite_fem_structure = None

    elif config["structural_solver"] == "kite_fem":
        ### kite_fem -- https://github.com/awegroup/kite_fem
        (
            kite_fem_structure,
            kite_fem_initial_conditions,
            kite_fem_pulley_matrix,
            kite_fem_spring_matrix,
            struc_nodes_initial,
        ) = structural_kite_fem.instantiate(
            config,
            struc_geometry,
            struc_nodes,
            kite_connectivity_arr,
            l0_arr,
            k_arr,
            c_arr,
            m_arr,
            linktype_arr,
            pulley_line_to_other_node_pair_dict,
        )
        # setting psm related output to None
        psystem = None

    else:
        raise ValueError("Invalid structural solver specified, either pss or pyfe3d")

    ##################
    ### AERO2STRUC ###
    ##################
    aero2struc_mapping = aero2struc.initialize_mapping(
        body_aero.panels,
        struc_nodes,
        struc_node_le_indices,
        struc_node_te_indices,
    )

    #################
    ### ACTUATION ###
    #################
    initial_length_power_tape = l0_arr[power_tape_index]
    power_tape_extension_step = config["power_tape_extension_step"]
    power_tape_final_extension = config["power_tape_final_extension"]
    if power_tape_extension_step != 0:
        n_power_tape_steps = int(power_tape_final_extension / power_tape_extension_step)
    else:
        n_power_tape_steps = 0
    logging.info(f"Initial depower tape length: {l0_arr[power_tape_index]:.3f}m")
    logging.info(
        f"Desired depower tape length: {initial_length_power_tape + power_tape_final_extension:.3f}m"
    )

    if config["steering_tape_extension_step"] != 0:
        raise NotImplementedError(
            "Steering tape actuation not implemented yet, set steering_tape_extension_step to 0."
        )

        ##TODO: get steering tape actuation back in
        # initial_length_steering_left = l0_arr[steering_tape_indices[0]]
        # initial_length_steering_right = l0_arr[steering_tape_indices[1]]
        # steering_tape_extension_step = config["steering_tape_extension_step"]
        # steering_tape_final_extension = config["steering_tape_final_extension"]
        # if steering_tape_extension_step != 0:
        #     n_steering_tape_steps = int(
        #         steering_tape_final_extension / steering_tape_extension_step
        #     )
        # else:
        #     n_steering_tape_steps = 0

        # if config["structural_solver"] == "pss":
        #     psystem.update_rest_length(
        #         steering_tape_indices[0], -steering_tape_final_extension
        #     )
        #     psystem.update_rest_length(
        #         steering_tape_indices[1], steering_tape_final_extension
        #     )
        # elif config["structural_solver"] == "kite_fem":
        #     raise ValueError("kite_fem solver is not implemented yet")

    ########################################
    ### AEROSTUCTURAL COUPLED SIMULATION ###
    ########################################
    tracking_data, meta = aerostructural_coupled_solver.main(
        m_arr=m_arr,
        struc_nodes=struc_nodes,
        struc_nodes_initial=struc_nodes_initial,
        config=config,
        ### ACTUATION
        initial_length_power_tape=initial_length_power_tape,
        n_power_tape_steps=n_power_tape_steps,
        power_tape_final_extension=power_tape_final_extension,
        power_tape_extension_step=power_tape_extension_step,
        ### CONNECTIVITY
        kite_connectivity_arr=kite_connectivity_arr,
        bridle_connectivity_arr=bridle_connectivity_arr,
        pulley_line_indices=pulley_line_indices,
        pulley_line_to_other_node_pair_dict=pulley_line_to_other_node_pair_dict,
        ### STRUC --> AERO
        struc_node_le_indices=struc_node_le_indices,
        struc_node_te_indices=struc_node_te_indices,
        ### AERO
        body_aero=body_aero,
        vsm_solver=vsm_solver,
        vel_app=vel_app,
        initial_polar_data=initial_polar_data,
        bridle_diameter_arr=bridle_diameter_arr,
        ### AERO --> STRUC
        aero2struc_mapping=aero2struc_mapping,
        power_tape_index=power_tape_index,
        ### STRUC
        psystem=psystem,
        kite_fem_structure=kite_fem_structure,
    )

    # Save results
    h5_path = Path(results_dir) / "sim_output.h5"
    save_results(tracking_data, meta, h5_path)

    # Load results
    meta_data_dict, tracking_data = load_sim_output(h5_path)

    # logging.info(f"meta_data: {meta_data_dict}")
    # - here you could add functions to plot the tracking of f_int, f_ext and f_residual over the iterations
    # - functions that make an animation of the kite going through the iterations
    # - etc.
    f_residual = tracking_data["f_int"] - tracking_data["f_ext"]

    printing_rest_lengths(tracking_data, struc_geometry)


if __name__ == "__main__":
    main()
