from src.post_processing import functions_print, functions_plot
from src.post_processing import post_processing_utils as post_processing_utils
from src.coupling import coupling_struc2aero, coupling_aero2struc
from src.aerodynamic import VSM, bridle_line_system_aero
from src.solver import solver_main
from src.structural import structural_model
from src.initialisation.yaml_loader import config
from src.initialisation.path_functions import load_module_from_path
import numpy as np
import imageio
from PIL import Image


def print_results(
    # mutables
    points,
    print_data,
    # immutables
    config,
):
    # Unpacking the print_data
    (
        is_convergence,
        num_of_iterations,
        aero_structural_total_time,
        vel_app,
        residual_f_including_fixed_nodes,
        residual_f,
        f_internal,
        f_external,
        aero_force_print_data,
        f_tether_drag,
        f_gravity,
        wing_rest_lengths,
        bridle_rest_lengths,
    ) = print_data
    force_aero, force_aero_wing, force_aero_bridle = aero_force_print_data

    ## Printing out the solution nicely
    functions_print.print_aerostructural(
        is_convergence,
        num_of_iterations,
        config.aero_structural.max_iter,
        aero_structural_total_time,
    )
    functions_print.print_settings(vel_app, config)
    functions_print.print_initial_kite_dimensions(config)
    functions_print.print_forces(
        f_internal,
        f_external,
        force_aero,
        residual_f_including_fixed_nodes,
        residual_f,
        f_tether_drag,
        f_gravity,
        points,
    )
    functions_print.print_aero(
        points, vel_app, force_aero_wing, force_aero_bridle, config
    )
    functions_print.print_structural(
        points, bridle_rest_lengths, wing_rest_lengths, config
    )


def plot(
    # mutables
    plot_data,
    points,
    vel_app,
    # immutables
    config,
):
    # Unpacking the plot_data
    wingpanels, controlpoints, rings, coord_L, F_rel = plot_data[0]
    wing_rest_lengths = plot_data[1]
    bridle_rest_lengths = plot_data[2]

    ##TODO: not working, doesn't show elongation
    ## extracting the elongation values of the bridle lines
    elongation_values = post_processing_utils.calculate_elongation(
        points,
        wing_rest_lengths,
        bridle_rest_lengths,
        config,
    )[2]
    functions_plot.plot_aero(
        points,
        elongation_values,
        vel_app,
        wingpanels,
        controlpoints,
        rings,
        coord_L,
        F_rel,
        config,
        elev=9.14,  # 10
        azim=0,  # 230
    )
    # functions_plot.plot_aero(config.vel_wind, config.vel_kite, vel_app, wingpanels, controlpoints, rings, coord_L, config.kite.n_segments, F_rel, elev=0, azim=-90)

    if config.is_with_plotly_plot:
        # Plotting the solution - offline PLOTLY
        plot_points_ini = [2, config.kite.points_ini, "black", 2]
        plot_points = [2, points, "green", 4]
        plot_lines_old = [
            config.kite.points_ini,
            np.concatenate(
                (config.kite.connectivity.bridle_ci, config.kite.connectivity.wing_ci)
            ),
            np.concatenate(
                (config.kite.connectivity.bridle_cj, config.kite.connectivity.wing_cj)
            ),
            "grey",
            1.5,
            [],
        ]
        plot_lines_new = [
            points,
            np.concatenate(
                (config.kite.connectivity.bridle_ci, config.kite.connectivity.wing_ci)
            ),
            np.concatenate(
                (config.kite.connectivity.bridle_cj, config.kite.connectivity.wing_cj)
            ),
            "purple",
            4,
            elongation_values,
        ]
        plot_surface_old = [
            True,
            config.kite.points_ini,
            config.kite.connectivity.plate_point_indices,
            "purple",
            0.0,
        ]
        plot_surface_new = [
            True,
            points,
            config.kite.connectivity.plate_point_indices,
            "lightgrey",
            0.3,
        ]
        functions_plot.plot_kite(
            [plot_points, plot_points_ini],
            [plot_lines_new, plot_lines_old],
            f" ",
            [plot_surface_new, plot_surface_old],
        )
        # functions_plot.plot_kite([plot_points],[plot_lines_new],f"Light-grey: surfplan, Purple: simulation",[plot_surface_new])


def animate(animation_data, vel_app, operating_condition, config, input_VSM):
    # Unpacking animation_data
    position, num_of_iterations, wing_rest_lengths, bridle_rest_lengths = animation_data

    n = len(config.kite.points_ini)
    num_frames = num_of_iterations
    print(f"Number of frames: {num_frames}")

    # Generate each frame
    for frame in range(num_frames):
        points = np.array(
            [
                [
                    position[f"x{n_i + 1}"].iloc[frame],
                    position[f"y{n_i + 1}"].iloc[frame],
                    position[f"z{n_i + 1}"].iloc[frame],
                ]
                for n_i in range(n)
            ]
        )

        # Update the data arguments for your method
        # This will call your method again with new data
        # computing the aero-again as its somehow not correctly stored ##TODO: fix-this
        # Struc --> aero
        points_left_to_right = coupling_struc2aero.order_struc_nodes_right_to_left(
            points, config.kite.connectivity.plate_point_indices
        )
        # Wing Aerodynamic
        (
            force_aero_wing_VSM,
            moment_aero_wing_VSM,
            F_rel,
            ringvec,
            controlpoints,
            wingpanels,
            rings,
            coord_L,
            coord_refined,
        ) = VSM.calculate_force_aero_wing_VSM(points_left_to_right, vel_app, input_VSM)

        elongation_values = post_processing_utils.calculate_elongation(
            points,
            wing_rest_lengths,
            bridle_rest_lengths,
            config,
        )[2]

        functions_plot.plot_aero(
            points,
            elongation_values,
            vel_app,
            wingpanels,
            controlpoints,
            rings,
            coord_L,
            F_rel,
            config,
            elev=10,
            azim=230,
            it_number=frame,
        )

    # Use pillow to save all frames as an animation in a gif file

    output_path = config.output_path

    # get_centroid = load_module_from_path(
    #     KITE_NAME, f"{folder_path_initialisation}/functions_wing.py"
    #     ).get_centroid

    images = [
        Image.open(f"{output_path}/animation/plot_iteration_{frame}.png")
        for frame in range(num_frames)
    ]

    # images[0].save(f"results/{config.kite_name}/animation/{operating_condition}_animation_1.gif", save_all=True, append_images=images[1:], duration=100, loop=0)

    # Assuming images is a list of PIL Image objects
    imageio.mimsave(
        f"{output_path}/animation/{operating_condition}_animation_1.mp4",
        images,
    )
