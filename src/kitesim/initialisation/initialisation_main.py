import numpy as np
from pathlib import Path

from kitesim.cases import cases_yaml_reader
from kitesim.initialisation import (
    input_psm,
    setup_kite_dict,
    initialisation_utils,
    input_vsm,
    input_bridle_aero,
)


def setup_config(
    path_config,
    path_processed_data_folder,
):
    """Setup the configuration for the simulation.

        Steps:
        1. Load the yaml defined settings into dicts
        2. Create extra dicts for some of the child classes
        3. Create a big nested dict_config_data with all of these dicts
        4. Dynamically create a frozen attrs nested config class from this big dictionary
    Args:
        path_config (str): The path to the configuration file.
        path_processed_data_folder (str): The path to the processed data folder.

    Returns:
        The configuration class.
    """

    # Load the yaml defined settings into dicts
    config_data = cases_yaml_reader.read_yaml_file(path_config, None)
    config_data.update(cases_yaml_reader.read_yaml_file(None, "default_case.yaml"))
    config_data.update(
        cases_yaml_reader.read_yaml_file(None, f"{config_data['sim_name']}.yaml")
    )
    path_kite_config = (
        Path(path_processed_data_folder)
        / config_data["kite_name"]
        / f"config_kite_{config_data['kite_name']}.yaml"
    )
    path_kite_data = (
        Path(path_processed_data_folder)
        / str(config_data["kite_name"])
        / "processed_design_files"
    )
    config_data_kite = cases_yaml_reader.read_yaml_file(path_kite_config, None)
    dict_kite_config = setup_kite_dict.process_kite_yamls_into_dict(
        config_data,
        config_data_kite,
        path_kite_data,
    )
    config_data.update(dict_kite_config)

    # Correcting certain dict items, ensuring the float and np.array types
    config_data["vel_wind"] = np.array([float(i) for i in config_data["vel_wind"]])
    config_data["vel_kite"] = np.array([float(i) for i in config_data["vel_kite"]])
    config_data["acc_kite"] = np.array([float(i) for i in config_data["acc_kite"]])
    config_data["plot_elev"] = np.array([float(i) for i in config_data["plot_elev"]])
    config_data["plot_azim"] = np.array([float(i) for i in config_data["plot_azim"]])

    return initialisation_utils.create_attr_class_from_dict("config", config_data)


def get_sim_input(
    path_config,
    path_processed_data_folder,
):
    """Get the simulation input.

    Args:
        path_config (str): The path to the configuration file.
        path_processed_data_folder (str): The path to the processed data folder.

    Returns:
        The simulation input, consisting of:
            points: The initial points of the kite.
            vel_app: The apparent velocity of the kite.
            config: The configuration class.
            input_bridle_aero: The input bridle aero class.
            input_VSM: The input VSM class.
            input_PSM: The input PSM class."""

    # Setting up the config class
    config = setup_config(
        path_config,
        path_processed_data_folder,
    )

    # Defining input_dict
    sim_input = {
        "points": config.kite.points_ini,
        "vel_app": np.array(config.vel_wind - config.vel_kite),
        "config": config,
        "input_bridle_aero": input_bridle_aero.create_input_bridle_aero(config),
        "input_tether_aero": input_bridle_aero.create_input_bridle_aero(config),
        "input_VSM": input_vsm.create_input_VSM(config),
        "input_PSM": input_psm.create_input_PSM(config),
    }

    return sim_input
