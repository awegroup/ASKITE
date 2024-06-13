import os
import yaml


def read_yaml_file(
    file_path: str,
    yaml_file_name: str,
) -> dict:
    """Reads a yaml file and returns its content as a dictionary

    Args:
        file_path: Path to the yaml file
        yaml_file_name: Name of the yaml file

    Returns:
        dict: Content of the yaml file
    """
    if file_path is None:
        yaml_folder = os.path.dirname(__file__)
        file_path = os.path.join(yaml_folder, yaml_file_name)

    if not os.path.exists(file_path):
        raise FileNotFoundError(f"No such file: {file_path}")

    with open(file_path, "r") as file:
        dict_filled_with_yaml = yaml.load(file, Loader=yaml.SafeLoader)
    return dict_filled_with_yaml
