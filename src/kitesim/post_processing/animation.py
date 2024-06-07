# import numpy as np
# import imageio
# from PIL import Image

# from kitesim.coupling import coupling_struc2aero
# from kitesim.aerodynamic import VSM
# from kitesim.post_processing import post_processing_utils
# from kitesim.post_processing import plotting


# def make_animation(loaded_data, path_run_results_folder):

#     # Unpacking loaded_data
#     vel_app = loaded_data["vel_app"]
#     config = loaded_data["config"]
#     input_VSM = loaded_data["input_VSM"]
#     position = loaded_data["position"]
#     num_of_iterations = loaded_data["num_of_iterations"]
#     wing_rest_lengths = loaded_data["wing_rest_lengths"]
#     bridle_rest_lengths = loaded_data["bridle_rest_lengths"]

#     n = len(config.kite.points_ini)
#     num_frames = num_of_iterations
#     print(f"Number of frames: {num_frames}")
#     vel_app_norm = np.linalg.norm(vel_app)

#     # Generate each frame
#     for frame in range(num_frames):
#         points = np.array(
#             [
#                 [
#                     position[f"x{n_i + 1}"].iloc[frame],
#                     position[f"y{n_i + 1}"].iloc[frame],
#                     position[f"z{n_i + 1}"].iloc[frame],
#                 ]
#                 for n_i in range(n)
#             ]
#         )

#         # Update the data arguments for your method
#         # This will call your method again with new data
#         # computing the aero-again as its somehow not correctly stored ##TODO: fix-this
#         # Struc --> aero
#         points_left_to_right = coupling_struc2aero.order_struc_nodes_right_to_left(
#             points, config.kite.connectivity.plate_point_indices
#         )
#         # Wing Aerodynamic
#         (
#             force_aero_wing_VSM,
#             moment_aero_wing_VSM,
#             F_rel,
#             ringvec,
#             controlpoints,
#             wingpanels,
#             rings,
#             coord_L,
#             coord_refined,
#         ) = VSM.calculate_force_aero_wing_VSM(points_left_to_right, vel_app, input_VSM)

#         elongation_values = post_processing_utils.calculate_elongation(
#             points,
#             wing_rest_lengths,
#             bridle_rest_lengths,
#             config,
#         )[2]

#         plotting.plot_aero(
#             points,
#             elongation_values,
#             vel_app,
#             wingpanels,
#             controlpoints,
#             rings,
#             coord_L,
#             F_rel,
#             config,
#             path_run_results_folder,
#             elev=10,
#             azim=-90,  # 230,
#             it_number=frame,
#         )

#     # Use pillow to save all frames as an animation in a gif file

#     # get_centroid = load_module_from_path(
#     #     KITE_NAME, f"{folder_path_initialisation}/functions_wing.py"
#     #     ).get_centroid

#     images = [
#         Image.open(f"{path_run_results_folder}/animation/plot_iteration_{frame}.png")
#         for frame in range(num_frames)
#     ]

#     # images[0].save(f"results/{config.kite_name}/animation/{config.sim_name}_animation_1.gif", save_all=True, append_images=images[1:], duration=100, loop=0)

#     # Assuming images is a list of PIL Image objects
#     imageio.mimsave(
#         f"{path_run_results_folder}/animation/{config.sim_name}_animation_va_{vel_app_norm:.1f}.mp4",
#         images,
#     )
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from kitesim.coupling import coupling_struc2aero
from kitesim.aerodynamic import VSM
from kitesim.post_processing import post_processing_utils
from kitesim.post_processing import plotting


def make_animation(loaded_data: dict, path_run_results_folder: str):
    """Create an animation of the kite simulation.

    Args:
        loaded_data (dict): Dictionary containing the loaded data.
        path_run_results_folder (str): The path to the folder where the results will be stored.

    Returns:
        None
    """

    # Unpacking loaded_data
    vel_app = loaded_data["vel_app"]
    config = loaded_data["config"]
    input_VSM = loaded_data["input_VSM"]
    position = loaded_data["position"]
    num_of_iterations = loaded_data["num_of_iterations"]
    wing_rest_lengths = loaded_data["wing_rest_lengths"]
    bridle_rest_lengths = loaded_data["bridle_rest_lengths"]

    # Creating the folder
    if not os.path.exists(f"{path_run_results_folder}/animation"):
        os.makedirs(f"{path_run_results_folder}/animation")

    # Other parameters
    n = len(config.kite.points_ini)
    num_frames = num_of_iterations
    print(f"Number of frames: {num_frames}")
    vel_app_norm = np.linalg.norm(vel_app)
    animation_elev = config.animation_elev
    animation_azim = config.animation_azim

    # Handle the case where there are no frames to animate
    if num_frames == 0:
        print("Error: Number of frames is zero. Exiting function.")
        return

    # Precompute connectivity
    plate_point_indices = config.kite.connectivity.plate_point_indices

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

    def update(frame):
        ax.clear()

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

        # Struc --> aero
        points_left_to_right = coupling_struc2aero.order_struc_nodes_right_to_left(
            points, plate_point_indices
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

        plotting.plot_aero(
            points,
            elongation_values,
            vel_app,
            wingpanels,
            controlpoints,
            rings,
            coord_L,
            F_rel,
            config,
            path_run_results_folder,
            elev=animation_elev,
            azim=animation_azim,
            it_number=frame,
            ax=ax,  # Pass the axes object
        )

    anim = FuncAnimation(fig, update, frames=num_frames, repeat=False)
    # Save the animation with desired quality settings
    anim.save(
        f"{path_run_results_folder}/animation/{config.sim_name}_animation_va_{vel_app_norm:.1f}.mp4",
        writer="ffmpeg",
        fps=config.animation_fps,  # Adjust fps for frame rate
        dpi=config.animation_dpi,  # Adjust dpi for resolution
        bitrate=config.animation_bitrate,  # Adjust bitrate for quality
        extra_args=["-vcodec", "libx264"],  # Use H.264 codec for high quality
    )

    plt.close(fig)
