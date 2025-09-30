import time
from tqdm import tqdm
import numpy as np
import logging
from pathlib import Path
import copy
from kitesim import (
    aerodynamic,
    struc2aero,
    aero2struc,
    structural_pss,
    tracking,
    plotting,
    structural_kite_fem,
)


# TODO: remove hardcoded values, when changing away from V3
def forcing_symmetry(struc_nodes):
    """
    Forcing symmetry in the y-direction for the kite structure nodes.
    This is a temporary solution to ensure symmetry in the simulation.
    """
    symmetry_pairs_dict = {
        1: 19,
        2: 20,
        3: 17,
        4: 18,
        5: 15,
        6: 16,
        7: 13,
        8: 14,
        9: 11,
        10: 12,
        # bridles
        21: 24,
        22: 23,
        25: 26,
        27: 30,
        28: 29,
        31: 32,
        33: 35,
        36: 37,
    }

    for key, value in symmetry_pairs_dict.items():
        struc_nodes[value] = np.array(
            [struc_nodes[key][0], -struc_nodes[key][1], struc_nodes[key][2]]
        )
    struc_nodes[34][1] = 0
    return struc_nodes


def main(
    m_arr=None,
    struc_nodes=None,
    psystem=None,
    struc_node_le_indices=None,
    struc_node_te_indices=None,
    body_aero=None,
    vsm_solver=None,
    vel_app=None,
    initial_polar_data=None,
    aero2struc_mapping=None,
    pss_connectivity=None,
    pss_params=None,
    power_tape_index=None,
    initial_length_power_tape=None,
    n_power_tape_steps=None,
    power_tape_final_extension=None,
    power_tape_extension_step=None,
    config=None,
    pulley_line_indices=None,
    pulley_line_to_other_node_pair_dict=None,
    ### kite_fem specifics ###
    kite_fem_structure=None,
    kite_fem_initial_conditions=None,
    kite_fem_pulley_matrix=None,
    kite_fem_spring_matrix=None,
    l0_arr=None,
):
    """
    Runs the aero-structural solver for the given input parameters.

    Args:
        config (dict): Main configuration dictionary.
        PROJECT_DIR (Path): Path to the project directory.
        results_dir (Path): Path to the results directory.

    Returns:
        tracking_data (dict): Dictionary containing time histories of positions, forces, etc.
        meta (dict): Dictionary with meta information about the simulation (timing, convergence, etc).
    """

    ## PRELOOP
    if config["is_with_gravity"]:
        f_ext_gravity = np.array(
            [np.array(config["grav_constant"]) * m_pt for m_pt in m_arr]
        )
    else:
        f_ext_gravity = np.zeros(struc_nodes.shape)

    if config["structural_solver"] == "pss":
        initial_particles = copy.deepcopy(psystem.particles)
    elif config["structural_solver"] == "kite_fem":
        initial_particles = copy.deepcopy(kite_fem_initial_conditions)
        rest_lengths = kite_fem_structure.modify_get_spring_rest_length()
    else:
        raise ValueError("Invalid structural solver specified in config.")

    t_vector = np.linspace(
        1,
        config["aero_structural_solver"]["max_iter"],
        config["aero_structural_solver"]["max_iter"],
    )
    tracking_data = tracking.setup_tracking_arrays(len(struc_nodes), t_vector)
    is_convergence = False
    f_residual_list = []
    f_tether_drag = np.zeros(3)
    is_residual_below_tol = False
    struc_nodes_prev = None  # Initialize previous points for tracking
    start_time = time.time()
    plotting.set_plot_style()

    ## track initial state
    # Update unified tracking dataframe (replaces position update)
    tracking.update_tracking_arrays(
        tracking_data,
        0,
        struc_nodes,
        np.zeros(np.shape(struc_nodes.flatten())),
        np.zeros(np.shape(struc_nodes.flatten())),
    )

    ######################################################################
    # SIMULATION LOOP
    ######################################################################
    ## propagating the simulation for each timestep and saving results
    with tqdm(total=len(t_vector), desc="Simulating", leave=True) as pbar:
        for i, step in enumerate(t_vector):
            if i > 0:
                struc_nodes_prev = struc_nodes.copy()

            ## external force
            begin_time_f_ext = time.time()

            ### STRUC --> AERO
            le_arr, te_arr = struc2aero.main(
                struc_nodes,
                struc_node_le_indices,
                struc_node_te_indices,
                config["aerodynamic"]["n_aero_panels_per_struc_section"],
            )

            ### AERO
            f_aero_wing_vsm_format, body_aero, results_aero = (
                aerodynamic.run_vsm_package(
                    body_aero=body_aero,
                    solver=vsm_solver,
                    le_arr=le_arr,
                    te_arr=te_arr,
                    va_vector=vel_app,
                    aero_input_type="reuse_initial_polar_data",
                    initial_polar_data=initial_polar_data,
                    is_with_plot=config["is_with_aero_plot_per_iteration"],
                )
            )
            logging.debug(
                f"Aero symmetry check, f_aero_y: {np.sum([force[1] for force in f_aero_wing_vsm_format])}"
            )
            ### AERO --> STRUC
            f_aero_wing = aero2struc.main(
                config["aero2struc"]["coupling_method"],
                f_aero_wing_vsm_format,
                struc_nodes,
                np.array(results_aero["panel_cp_locations"]),
                aero2struc_mapping,
                p=2,
                eps=1e-6,
                is_with_coupling_plot=config["is_with_coupling_plot_per_iteration"],
            )

            ## TODO: get bridle line forces back in to play
            f_aero_bridle = np.zeros(struc_nodes.shape)
            f_aero = f_aero_wing + f_aero_bridle

            ## summing up
            f_ext = f_aero + f_ext_gravity
            f_ext = np.round(f_ext, 5)

            if config["is_with_struc_plot_per_iteration"]:
                if config["structural_solver"] == "pss":
                    rest_lengths = psystem.extract_rest_length
                elif config["structural_solver"] == "kite_fem":
                    rest_lengths = kite_fem_structure.modify_get_spring_rest_length()
                plotting.main(
                    struc_nodes,
                    pss_connectivity,
                    rest_lengths,
                    f_ext=f_ext,
                    title=f"i: {i}",
                    body_aero=body_aero,
                    chord_length=2.6,
                    is_with_node_indices=False,
                    pulley_line_indices=pulley_line_indices,
                    pulley_line_to_other_node_pair_dict=pulley_line_to_other_node_pair_dict,
                )

            ### STRUC
            f_ext_flat = f_ext.flatten()
            end_time_f_ext = time.time()
            begin_time_f_int = time.time()
            if config["structural_solver"] == "pss":
                psystem, converged = structural_pss.run_pss(
                    psystem, pss_params, f_ext_flat
                )
                logging.debug(f"PS converged: {converged}")
                # logging.debug(f"position.loc[step].shape: {position.loc[step].shape}")
                logging.debug(f"internal force: {psystem.f_int}")
                logging.debug(f"external force: {f_ext}")
                # Updating the points
                struc_nodes = np.array([particle.x for particle in psystem.particles])
                # Extracting internal force
                f_int = psystem.f_int

            elif config["structural_solver"] == "kite_fem":
                # reset to this iteration (sets coords_current to initial coords)
                kite_fem_structure.reset()

                # [fx, fy, fz, mx, my, mz] for each node
                f_ext_reshaped = f_ext_flat.reshape(-1, 3)
                fe_6d = [[fe[0], fe[1], fe[2], 0, 0, 0] for fe in f_ext_reshaped]
                fe_6d = np.array(fe_6d).flatten()

                # plotting
                # kite_fem_structure.plot_3D()
                # print(f"fe_6d: {fe_6d}")
                # breakpoint()

                kite_fem_structure.solve(
                    fe=fe_6d,
                    max_iterations=config["structural_kite_fem"]["max_iterations"],
                    tolerance=config["structural_kite_fem"]["tolerance"],
                    step_limit=config["structural_kite_fem"]["step_limit"],
                    relax_init=config["structural_kite_fem"]["relax_init"],
                    relax_update=config["structural_kite_fem"]["relax_update"],
                    k_update=config["structural_kite_fem"]["k_update"],
                    I_stiffness=config["structural_kite_fem"]["I_stiffness"],
                )
                struc_nodes = kite_fem_structure.coords_current
                # reshape from flat to (n_nodes, 3)
                struc_nodes = struc_nodes.reshape(-1, 3)
                f_int = kite_fem_structure.fi
                # remove moments
                f_int = f_int.reshape(-1, 6)[:, :3].flatten()
                ##TODO: get convergence status from kite_fem_structure.solve
                converged = True

                ##TODO: implement kite_fem solver
                # need to change the fe to a 6d vector, with fx,fy,fz and mx, my, mz
                # 2. fe_6d = [[fe[0],fe[1],fe[2],0,0,0] for fe in f_ext_flat]
                # 3. fem_structure.solve(
                #   fe=f_6d,
                #   config['structural_kite_fem']['max_iterations'],
                #   config['structural_kite_fem']['tolerance'],
                #   config['structural_kite_fem']['step_limit'],
                #   config['structural_kite_fem']['relax_init'],
                #   config['structural_kite_fem']['relax_update'],
                #   config['structural_kite_fem']['k_update'],
                #   config['structural_kite_fem']['I_stiffness'],
                # )
                # --> when convergence
                # coordinates_flat_array = fem_structure.coords_current
                # 4. fem_structure.coords_init()
                # 5. change current with re-initilasiatoinaoginoeng function

                # raise ValueError("kite_fem solver is not implemented yet")
                # kite_fem_structure.plot_convergence()
                # kite_fem_structure.plot_3D(fe=fe_6d)

            end_time_f_int = time.time()

            # Forcing symmetry
            if config["is_with_forcing_symmetry"]:
                logging.info("Forcing symmetry in y-direction")
                struc_nodes = forcing_symmetry(struc_nodes)

            # Computing residual forces
            f_residual = f_int + f_ext_flat
            f_residual_list.append(np.linalg.norm(np.abs(f_residual)))
            logging.debug(
                f"residual force in y-direction: {np.sum([f_residual[1::3]]):.3f}N"
            )

            # Update unified tracking dataframe (replaces position update)
            tracking.update_tracking_arrays(
                tracking_data,
                i,
                struc_nodes,
                f_ext_flat,
                f_residual,
            )
            ## calculating delta tape lengths
            if config["structural_solver"] == "pss":
                delta_power_tape = (
                    psystem.extract_rest_length[power_tape_index]
                    - initial_length_power_tape
                )
            ##TODO: implement depower tape length change for kite_fem
            ##TODO: fix this
            elif config["structural_solver"] == "kite_fem":
                delta_power_tape = (
                    rest_lengths[power_tape_index] + power_tape_extension_step
                )
                rest_lengths = kite_fem_structure.modify_get_spring_rest_length(
                    spring_ids=[power_tape_index], new_l0s=[new_l0_depower_tape]
                )
                delta_power_tape = 0
                print(f"NOT IMPLEMENTED: delta_power_tape for kite_fem")

            pbar.set_postfix(
                {
                    "res": f"{np.linalg.norm(f_residual):.3f}N",
                    "aero": f"{end_time_f_ext-begin_time_f_ext:.2f}s",
                    "struc": f"{end_time_f_int-begin_time_f_int:.2f}s",
                }
            )
            pbar.update(1)

            ### All the convergence checks, are be done in if-elif because only 1 should hold at once
            # if convergence (residual below set tolerance)
            if (
                i > n_power_tape_steps
                and np.linalg.norm(f_residual)
                <= config["aero_structural_solver"]["tol"]
            ):
                is_residual_below_tol = True
                is_convergence = True
            # if np.linalg.norm(f_residual) <= pss_params["aerostructural_tol"]l and i > 50:
            #     is_residual_below_tol = True
            #     is_convergence = True

            # if residual forces are NaN
            elif np.isnan(np.linalg.norm(f_residual)):
                is_convergence = False
                logging.info("Classic PS diverged - residual force is NaN")
                break
            # if residual forces are not changing anymore
            elif (
                i > 15
                and np.abs(np.mean(f_residual_list[i - 15]) - f_residual_list[i - 10])
                < 1
            ):
                is_convergence = False
                logging.info("Classic PS non-converging - residual no longer changes")
                break
            # if to many iterations are needed
            elif i > config["aero_structural_solver"]["max_iter"]:
                is_convergence = False
                logging.info(
                    f"Classic PS non-converging - more than max ({config['aero_structural_solver']['max_iter']}) iterations needed"
                )
                break
            # special case for running the simulation for only one timestep
            elif config["is_run_only_1_time_step"]:
                break
            # when aero does not converge
            elif np.sum([force[1] for force in f_aero_wing_vsm_format]) == np.nan:
                is_convergence = False
                logging.info("Classic PS non-converging - aero forces are NaN")
                break

            ## if convergenced and not yet at desired depower-tape length
            elif (
                is_residual_below_tol and delta_power_tape <= power_tape_final_extension
            ):
                psystem.update_rest_length(power_tape_index, power_tape_extension_step)
                delta_power_tape = (
                    psystem.extract_rest_length[power_tape_index]
                    - initial_length_power_tape
                )
                logging.info(
                    f"||--- delta l_d: {delta_power_tape:.3f}m | new l_d: {psystem.extract_rest_length[power_tape_index]:.3f}m | Steps required: {n_power_tape_steps}"
                )
                is_residual_below_tol = False

            ## if convergenced and all changes are made
            elif is_residual_below_tol:
                is_convergence = True

            if is_convergence:
                break
    ######################################################################
    ## END OF SIMULATION FOR LOOP
    ######################################################################
    # TODO: move reshaping and so on to structural solver files
    if config["structural_solver"] == "pss":
        current_particles = np.array([particle.x for particle in psystem.particles])
        print(f"psystem: {psystem}")
        rest_lengths = psystem.extract_rest_length
        particles_ini = np.array(initial_particles)
        # reshape to current_particles shape
        particles_ini = initial_particles.reshape(current_particles.shape)
    elif config["structural_solver"] == "kite_fem":
        current_particles = struc_nodes
        rest_lengths = kite_fem_structure.modify_get_spring_rest_length()
        # For kite_fem, initial_particles is already a flat array, reshape it properly
        if isinstance(initial_particles, np.ndarray):
            # If it's already a numpy array, reshape it to match current_particles
            particles_ini = initial_particles.reshape(-1, 6)[:, :3]
            particles_ini = particles_ini.reshape(current_particles.shape)
        else:
            # If it's still a list/other format, convert appropriately
            particles_ini = np.array(initial_particles).reshape(-1, 6)[:, :3]
            particles_ini = particles_ini.reshape(current_particles.shape)
    if config["is_with_final_plot"]:

        plotting.main(
            current_particles,
            pss_connectivity,
            f_ext=f_ext,
            rest_lengths=rest_lengths,
            struc_nodes_initial=particles_ini,
            title="Initial  vs final",
            pulley_line_indices=pulley_line_indices,
            pulley_line_to_other_node_pair_dict=pulley_line_to_other_node_pair_dict,
        )
    # Convert kite_connectivity to a numeric array for HDF5 compatibility
    kite_connectivity_numeric = np.array(pss_connectivity)
    # If it is not 2D and numeric, try to extract only the first two columns (node indices)
    if kite_connectivity_numeric.dtype == object or kite_connectivity_numeric.ndim != 2:
        kite_connectivity_numeric = np.array(
            [[int(row[0]), int(row[1])] for row in pss_connectivity],
            dtype=np.int32,
        )
    meta = {
        "total_time_s": time.time() - start_time,
        "n_iter": i + 1,
        "converged": is_convergence,
        "rest_lengths": rest_lengths,  # ensure numeric array
        "kite_connectivity": kite_connectivity_numeric,
    }

    return tracking_data, meta
