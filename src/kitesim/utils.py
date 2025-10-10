import yaml
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull
import h5py


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


def load_and_save_config_files(
    config_path, struc_geometry_path, aero_geometry_path, results_dir
):
    """
    Load configuration files and save copies to a timestamped results directory.

    Args:
        PROJECT_DIR (Path): The project directory.

    Returns:
        config (dict): The loaded main configuration.
        config_kite (dict): The loaded kite configuration.
        results_dir (Path): Path to the results directory where configs are saved.
    """
    config = load_yaml(config_path)
    struc_geometry = load_yaml(struc_geometry_path)
    aero_geometry = load_yaml(aero_geometry_path)

    results_dir.mkdir(parents=True, exist_ok=True)
    with open(results_dir / "config.yaml", "w") as f:
        yaml.dump(config, f, sort_keys=False)
    with open(results_dir / "struc_geometry.yaml", "w") as f:
        yaml.dump(struc_geometry, f, sort_keys=False)
    with open(results_dir / "aero_geometry.yaml", "w") as f:
        yaml.dump(aero_geometry, f, sort_keys=False)

    return config, struc_geometry, aero_geometry, results_dir


def save_results(tracking, meta, filename):
    """
    Save tracking arrays and metadata to an HDF5 file.

    Args:
        tracking (dict): Dictionary of arrays to save under the "tracking" group.
        meta (dict): Metadata dictionary to save as attributes.
        filename (str or Path): Output HDF5 file path.

    Returns:
        None
    """
    with h5py.File(filename, "w") as f:
        grp = f.create_group("tracking")
        for name, arr in tracking.items():
            grp.create_dataset(name, data=arr[: meta["n_iter"]], compression="gzip")
        for k, v in meta.items():
            grp.attrs[k] = v


def load_sim_output(h5_path):
    """
    Load simulation results and metadata from an HDF5 file written with h5py.

    Args:
        h5_path (str or Path): Path to the .h5 file (e.g. "sim_output.h5").

    Returns:
        tuple:
            metadata (dict): Run-level metadata (attributes from the file).
            track (dict): Dictionary of numpy arrays for each dataset under "tracking".
    """
    h5_path = Path(h5_path)
    if not h5_path.exists():
        raise FileNotFoundError(f"No such file: {h5_path}")

    with h5py.File(h5_path, "r") as f:
        if "tracking" not in f:
            raise KeyError(f"No 'tracking' group in {h5_path}")
        grp = f["tracking"]

        # load metadata
        metadata = {key: grp.attrs[key] for key in grp.attrs}

        # load all datasets under /tracking into numpy arrays
        track = {}
        for name, item in grp.items():
            if isinstance(item, h5py.Dataset):
                track[name] = item[()]  # read the full array into memory

    return metadata, track


# TODO: at this moment unused
def calculate_projected_area(points):
    """
    Calculate the projected area of a set of 3D points onto the XY plane using the convex hull.

    Args:
        points (np.ndarray): Array of 3D points (n_points, 3).

    Returns:
        float: Projected area on the XY plane.
    """
    # Project points onto the x,y plane
    xy_points = points[:, :2]

    # Find the convex hull
    hull = ConvexHull(xy_points)
    hull_points = xy_points[hull.vertices]

    # Using the shoelace formula
    x = hull_points[:, 0]
    y = hull_points[:, 1]

    return 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))


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
                    total_length += np.linalg.norm(struc_nodes[ci] - struc_nodes[cj])
                    total_length += np.linalg.norm(struc_nodes[cj] - struc_nodes[ck])
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
                    total_length += np.linalg.norm(struc_nodes[ci] - struc_nodes[cj])
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
