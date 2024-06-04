import os
import yaml


def read_yaml_file(yaml_file_name):
    yaml_folder = os.path.dirname(__file__)
    file_path = os.path.join(yaml_folder, yaml_file_name)

    if not os.path.exists(file_path):
        raise FileNotFoundError(f"No such file: {file_path}")

    with open(file_path, "r") as file:
        return yaml.safe_load(file)
