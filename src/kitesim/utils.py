import yaml
from pathlib import Path
import numpy as np
import os
from datetime import datetime
from scipy.spatial import ConvexHull


def load_yaml(path: Path) -> dict:
    """
    Read a YAML file and return the parsed data as a Python dict.
    Args:
        path (Path): The path to the YAML file.
    Returns:
        dict: The parsed data from the YAML file.
    """
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def load_and_save_config_files(PROJECT_DIR):
    """
    Load and save the configuration files.
    Args:
        PROJECT_DIR (Path): The project directory.
    Returns:
        config (dict): The loaded configuration.
        config_kite (dict): The loaded kite configuration.
    """
    config = load_yaml(Path(PROJECT_DIR) / "data" / "config.yaml")
    config_kite = load_yaml(
        Path(PROJECT_DIR) / "data" / f"{config['kite_name']}" / "config_kite.yaml"
    )

    # create a results folder on this date and time and save the config files
    path_results_folder = (
        Path(PROJECT_DIR)
        / "results"
        / f'{config["kite_name"]}'
        / f'{datetime.now().strftime("%Y_%m_%d_%H")}h'
    )
    path_results_folder.mkdir(parents=True, exist_ok=True)
    with open(path_results_folder / "config.yaml", "w") as f:
        yaml.dump(config, f, sort_keys=False)
    with open(path_results_folder / "config_kite.yaml", "w") as f:
        yaml.dump(config_kite, f, sort_keys=False)

    return config, config_kite


##TODO: at this point unused
def calculate_projected_area(points):
    # Project points onto the x,y plane
    xy_points = points[:, :2]

    # Find the convex hull
    hull = ConvexHull(xy_points)
    hull_points = xy_points[hull.vertices]

    # Using the shoelace formula
    x = hull_points[:, 0]
    y = hull_points[:, 1]

    return 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))
