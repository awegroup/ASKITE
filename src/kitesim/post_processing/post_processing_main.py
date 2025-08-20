import dill
import os
from datetime import datetime
from pathlib import Path
import matplotlib.pyplot as plt

from kitesim.post_processing import plotting
from kitesim.post_processing import printing
from kitesim.post_processing import animation


def create_results_folder(sim_input, path_results_folder: str) -> str:
    """Create a folder to store the results of the simulation.

    Args:
        path_results_folder (str): The path to the folder where the results will be stored.
        config (Config): The configuration object containing the simulation settings.

    Returns:
        str: The path to the created results folder.
    """
    config = sim_input["config"]

    path_results_folder_run = (
        Path(path_results_folder)
        / config.kite_name
        / datetime.now().strftime("%Y_%m_%d_%Hh")
        / config.sim_name
    )

    # Ensure the folder exists
    if not os.path.exists(path_results_folder_run):
        os.makedirs(path_results_folder_run)

    return path_results_folder_run


def saving_all_dict_entries(
    to_be_saved_dict: dict, sub_folder_name: str, path_results_folder_run: str
):
    """Save the output data for the simulation.

    Args:
        to_be_saved_dict (dict): Dictionary containing data to be saved.
        sub_folder_name (str): Extension for the path to save the data.
        path_results_folder_run (str): Path to the results folder.

    Returns:
        str: The path to the saved output.

    Raises:
        None

    """
    # Create the path to the folder in which it will be saved
    path_saved_folder = Path(path_results_folder_run) / sub_folder_name

    # Ensure the folder exists
    if not os.path.exists(path_saved_folder):
        os.makedirs(path_saved_folder)

    # Serialize and save all data with dill
    for name, data in to_be_saved_dict.items():
        with open(f"{path_saved_folder}/{name}.pkl", "wb") as f:
            dill.dump(data, f)

    return


def load_all_files_from_folder(folder_path: str) -> dict:
    """Load all .pkl files in the specified folder and return a dictionary with the loaded data.

    Args:
        folder_path (Path): Path to the folder containing the .pkl files.

    Returns:
        dict: Dictionary containing the loaded data. The keys are the filenames without the extension.
    """
    pkl_files = folder_path.glob("*.pkl")
    data = {}
    for file in pkl_files:
        with open(file, "rb") as f:
            filename = file.stem
            data[filename] = dill.load(f)

    return data


def processing_output(path_results_folder_run: str) -> dict:
    """Process the output data from the simulation.

        First load the data from the saved input and output files.
        Then print, plot, and animate the results based on the configuration settings.

    Args:
        path_results_folder_run (str): Path to the results folder.

    Returns:
        dict: Dictionary containing the results of the simulation.

    """
    # Load the input and output data
    path_saved_input = Path(path_results_folder_run) / "input"
    path_saved_output = Path(path_results_folder_run) / "output"
    loaded_data_input = load_all_files_from_folder(path_saved_input)
    loaded_data_output = load_all_files_from_folder(path_saved_output)
    loaded_data = {**loaded_data_input, **loaded_data_output}

    # Get the configuration settings
    config = loaded_data["config"]

    # Print, plot, and animate the results based on the configuration settings
    if config.is_with_printing:
        # TODO: Could think about also saving this as .txt file
        printing.print_results(loaded_data)

    if config.is_with_plotting:
        plotting.make_plot(loaded_data, path_results_folder_run)
        plt.show()

    if config.is_with_animation:
        print(f"")
        print("--> Generating ANIMATION \{*_*}/")
        animation.make_animation(loaded_data, path_results_folder_run)

    return loaded_data_input
