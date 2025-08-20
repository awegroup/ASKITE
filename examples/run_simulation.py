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
from kitesim.utils import load_and_save_config_files, load_sim_output, save_results
from kitesim import (
    aerodynamic,
    # struc2aero,
    aero2struc,
    structural_pss,
    structural_pyfe3d,
    # tracking,
    # plotting,
    aerostructural_coupled_solver,
    read_struc_geometry_yaml,
)


# Import modules
def main():
    """Main function"""
    PROJECT_DIR = Path(__file__).resolve().parents[1]
    kite_name = "TUDELFT_V3_KITE"  # the dir name with the relevant .yaml files

    # load config.yaml & geometry.yaml, save both, and return them as dicts
    config, struc_geometry, aero_geometry, results_dir = load_and_save_config_files(
        kite_name, PROJECT_DIR
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
    body_aero, vsm_solver, vel_app, initial_polar_data = aerodynamic.initialize(
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
        conn_arr,
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
    for idx, conn in enumerate(conn_arr):
        logging.info(
            f"conn_idx: {idx}: conn: {conn}, l0: {l0_arr[idx]}, k: {k_arr[idx]}, c: {c_arr[idx]}, linktype: {linktype_arr[idx]}"
        )

    if config["structural_solver"] == "pss":
        ## pss -- https://github.com/awegroup/Particle_System_Simulator
        ##TODO: Fix the comment below, it SHOULD read l0
        # Note: ParticleSystem doesn’t read l0_arr. SpringDamper sets l0
        # from the initial particle positions.
        # So l0_arr is a bookkeeping array for you, not used at instantiation.
        (
            psystem,
            pss_connectivity,
            pss_initial_conditions,
            pss_params,
        ) = structural_pss.instantiate(
            # yaml files
            config,
            struc_geometry,
            # node level
            struc_nodes,
            m_arr,
            # element_level
            conn_arr,
            l0_arr,
            k_arr,
            c_arr,
            linktype_arr,
            pulley_line_to_other_node_pair_dict,
        )
        if config["is_with_initial_structure_plot"]:
            structural_pss.plot_3d_kite_structure(
                struc_nodes,
                pss_connectivity,
                power_tape_index,
                fixed_nodes=struc_geometry["fixed_point_indices"],
                pulley_nodes=pulley_node_indices,
            )
    elif config["structural_solver"] == "pyfe3d":
        ### py3fe -- https://github.com/saullocastro/pyfe3d
        structural_pyfe3d.instantiate(
            config,
            struc_geometry,
            struc_nodes,
            conn_arr,
            l0_arr,
            k_arr,
            c_arr,
            m_arr,
            linktype_arr,
            pulley_line_to_other_node_pair_dict,
        )

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
    initial_length_steering_left = l0_arr[steering_tape_indices[0]]
    initial_length_steering_right = l0_arr[steering_tape_indices[1]]
    steering_tape_extension_step = config["steering_tape_extension_step"]
    steering_tape_final_extension = config["steering_tape_final_extension"]
    if steering_tape_extension_step != 0:
        n_steering_tape_steps = int(
            steering_tape_final_extension / steering_tape_extension_step
        )
    else:
        n_steering_tape_steps = 0

    if config["structural_solver"] == "pss":
        psystem.update_rest_length(
            steering_tape_indices[0], -steering_tape_final_extension
        )
        psystem.update_rest_length(
            steering_tape_indices[1], steering_tape_final_extension
        )
    elif config["structural_solver"] == "pyfe3d":
        raise ValueError("pyfe3d solver is not implemented yet")

    ########################################
    ### AEROSTUCTURAL COUPLED SIMULATION ###
    ########################################
    tracking_data, meta = aerostructural_coupled_solver.main(
        m_arr,
        struc_nodes,
        psystem,
        struc_node_le_indices,
        struc_node_te_indices,
        body_aero,
        vsm_solver,
        vel_app,
        initial_polar_data,
        aero2struc_mapping,
        pss_connectivity,
        pss_params,
        power_tape_index,
        initial_length_power_tape,
        n_power_tape_steps,
        power_tape_final_extension,
        power_tape_extension_step,
        config,
        pulley_line_indices,
        pulley_line_to_other_node_pair_dict,
    )

    # Save results
    h5_path = Path(results_dir) / "sim_output.h5"
    save_results(tracking_data, meta, h5_path)

    # Load results
    meta_data_dict, tracking_data = load_sim_output(h5_path)

    # logging.info(f"meta_data: {meta_data_dict}")

    # TODO:
    # - here you could add functions to plot the tracking of f_int, f_ext and f_residual over the iterations
    # - functions that make an animation of the kite going through the iterations
    # - etc.
    f_residual = tracking_data["f_int"] - tracking_data["f_ext"]

    def printing_rest_lengths(tracking_data, struc_geometry):
        """
        Print the current and initial rest lengths of all bridle lines defined in bridle_elements by
        averaging the lengths of all their segments in bridle_connections.

        For each connection:
        - if 3 nodes, sum ci-cj and cj-ck
        - if 2 nodes, sum ci-cj
        """
        positions = tracking_data["positions"]
        struc_nodes = positions[-1]  # current positions
        initial_struc_nodes = positions[0]  # initial positions

        bridle_elements_data = struc_geometry["bridle_elements"]["data"]
        bridle_line_names = [row[0] for row in bridle_elements_data]
        bridle_connections_data = struc_geometry["bridle_connections"]["data"]

        # YAML l0 lookup dictionary
        bridle_elements_yaml = {row[0]: row[1] for row in bridle_elements_data}

        print("\nCurrent bridle line lengths:")
        results = []
        for line_name in bridle_line_names:
            total_length = 0.0
            initial_total_length = 0.0
            count = 0

            for conn in bridle_connections_data:
                if conn[0] == line_name:
                    ci = int(conn[1])
                    cj = int(conn[2])
                    if len(conn) > 3 and conn[3] not in (None, "", 0):
                        ck = int(conn[3])
                        # current
                        total_length += np.linalg.norm(
                            struc_nodes[ci] - struc_nodes[cj]
                        )
                        total_length += np.linalg.norm(
                            struc_nodes[cj] - struc_nodes[ck]
                        )
                        # initial
                        initial_total_length += np.linalg.norm(
                            initial_struc_nodes[ci] - initial_struc_nodes[cj]
                        )
                        initial_total_length += np.linalg.norm(
                            initial_struc_nodes[cj] - initial_struc_nodes[ck]
                        )
                        count += 1
                    else:
                        # current
                        total_length += np.linalg.norm(
                            struc_nodes[ci] - struc_nodes[cj]
                        )
                        # initial
                        initial_total_length += np.linalg.norm(
                            initial_struc_nodes[ci] - initial_struc_nodes[cj]
                        )
                        count += 1

            if count > 0:
                avg_length = total_length / count
                initial_avg_length = initial_total_length / count
                delta_pct = (
                    100.0 * (avg_length - initial_avg_length) / initial_avg_length
                    if initial_avg_length != 0
                    else 0.0
                )

                # yaml l0 value
                yaml_l0 = bridle_elements_yaml.get(line_name, None)
                try:
                    yaml_l0_val = float(yaml_l0) if yaml_l0 is not None else None
                except Exception:
                    yaml_l0_val = yaml_l0

                delta_yaml_pct = (
                    100.0 * (avg_length - yaml_l0_val) / yaml_l0_val
                    if yaml_l0_val not in (None, 0)
                    else 0.0
                )

                results.append(
                    f"{line_name}: curr: {avg_length:.3f} m, "
                    f"initial: {initial_avg_length:.3f} m, "
                    f"delta: {delta_pct:+.2f} %, "
                    f"yaml: {yaml_l0_val} m, "
                    f"delta_yaml: {delta_yaml_pct:+.2f} %"
                )

        # print once
        for result in results:
            print(result)

    printing_rest_lengths(tracking_data, struc_geometry)


if __name__ == "__main__":
    main()
