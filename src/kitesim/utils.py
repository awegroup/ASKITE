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
    Print current and initial lengths of bridle lines by averaging the lengths
    of their segments in bridle_connections.

    Supports both legacy `bridle_elements` and newer `bridle_lines` YAML schemas.

    For each connection:
    - if 3 nodes, sum ci-cj and cj-ck
    - if 2 nodes, sum ci-cj
    """
    positions = tracking_data["positions"]
    struc_nodes = positions[-1]  # current positions
    initial_struc_nodes = positions[0]  # initial positions

    if "bridle_elements" in struc_geometry:
        bridle_defs_data = struc_geometry["bridle_elements"]["data"]
    elif "bridle_lines" in struc_geometry:
        bridle_defs_data = struc_geometry["bridle_lines"]["data"]
    else:
        raise KeyError("Expected 'bridle_elements' or 'bridle_lines' in struc_geometry")

    bridle_line_names = [row[0] for row in bridle_defs_data]
    bridle_connections_data = struc_geometry["bridle_connections"]["data"]

    # YAML l0 lookup dictionary
    bridle_elements_yaml = {row[0]: row[1] for row in bridle_defs_data}

    # Collect data rows: (line_name, curr_length, yaml_l0, delta_pct, initial_nodal_dist)
    rows = []
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
            curr_l = total_length / count
            init_nodal_dist = initial_total_length / count

            yaml_l0 = bridle_elements_yaml.get(line_name, None)
            try:
                yaml_l0_val = float(yaml_l0) if yaml_l0 is not None else None
            except Exception:
                yaml_l0_val = yaml_l0

            delta_pct = (
                100.0 * (curr_l - yaml_l0_val) / yaml_l0_val
                if yaml_l0_val not in (None, 0)
                else 0.0
            )

            rows.append((line_name, curr_l, yaml_l0_val, delta_pct, init_nodal_dist))

    if not rows:
        print("\nNo bridle lines with matching connections found.")
        return

    # Determine column widths for aligned output
    name_w = max(len(r[0]) for r in rows) if rows else 10

    # Print header
    print(
        f"\n{'Line':<{name_w}}   {'current_l':>10}   {'initial_l0_yaml':>16}   {'delta':>8}   {'initial_nodal_distance':>23}"
    )
    print(f"{'-' * name_w}   {'-' * 10}   {'-' * 16}   {'-' * 8}   {'-' * 23}")

    for name, curr_l, yaml_l0, delta_pct, init_dist in rows:
        yaml_str = f"{yaml_l0:.3f} m" if yaml_l0 is not None else "N/A"
        print(
            f"{name:<{name_w}}   {curr_l:>7.3f} m   {yaml_str:>16}   {delta_pct:>+7.2f}%   ({init_dist:>19.3f} m)"
        )


def _axis_to_unit_vector(axis, name="axis"):
    """
    Convert an axis definition to a unit 3D vector.

    Accepted forms:
      - labels: "x", "y", "z"
      - vectors: array-like with 3 numeric components
    """
    if isinstance(axis, str):
        axis_label = axis.strip().lower()
        axis_map = {
            "x": np.array([1.0, 0.0, 0.0]),
            "y": np.array([0.0, 1.0, 0.0]),
            "z": np.array([0.0, 0.0, 1.0]),
        }
        if axis_label not in axis_map:
            raise ValueError(
                f"Invalid {name} label '{axis}'. Allowed labels are 'x', 'y', 'z'."
            )
        vec = axis_map[axis_label]
    else:
        vec = _to_3_vector(axis, name)

    norm = np.linalg.norm(vec)
    if norm <= 1e-15:
        raise ValueError(f"{name} must be non-zero. Got {vec}.")
    return vec / norm


def _normalize_axes(axes):
    """
    Normalize a 3-axis sequence into three unit vectors.

    Supported input layouts:
      - ("x", "y", "z")
      - ([1,0,0], [0,1,0], [0,0,1])
      - np.array([[1,0,0],[0,1,0],[0,0,1]])
    """
    axes_arr = np.asarray(axes, dtype=object)
    if axes_arr.ndim == 1:
        if axes_arr.shape[0] != 3:
            raise ValueError(
                f"`axes` must contain exactly 3 axis definitions. Got {axes_arr.shape[0]}."
            )
        axis_defs = list(axes_arr)
    elif axes_arr.ndim == 2 and axes_arr.shape == (3, 3):
        axis_defs = [axes_arr[i, :] for i in range(3)]
    else:
        raise ValueError(
            "`axes` must be either length-3 labels/vectors or a (3,3) axis matrix."
        )

    return tuple(
        _axis_to_unit_vector(axis_defs[i], name=f"axes[{i}]") for i in range(3)
    )


def _rotation_matrix_from_axis(axis, angle_rad):
    """
    Return a 3x3 right-hand-rule rotation matrix about any 3D axis.

    Args:
        axis (str or array-like): "x"/"y"/"z" or a 3D axis vector.
        angle_rad (float): Rotation angle in radians.
    """
    ux, uy, uz = _axis_to_unit_vector(axis)
    c = np.cos(angle_rad)
    s = np.sin(angle_rad)
    one_minus_c = 1.0 - c

    return np.array(
        [
            [
                c + ux * ux * one_minus_c,
                ux * uy * one_minus_c - uz * s,
                ux * uz * one_minus_c + uy * s,
            ],
            [
                uy * ux * one_minus_c + uz * s,
                c + uy * uy * one_minus_c,
                uy * uz * one_minus_c - ux * s,
            ],
            [
                uz * ux * one_minus_c - uy * s,
                uz * uy * one_minus_c + ux * s,
                c + uz * uz * one_minus_c,
            ],
        ]
    )


def _to_3_vector(values, name):
    """Convert input to a strict 3-vector of floats."""
    arr = np.asarray(values, dtype=float).reshape(-1)
    if arr.shape != (3,):
        raise ValueError(f"{name} must contain exactly 3 values. Got shape {arr.shape}.")
    return arr


def rotate_geometry(
    struc_nodes,
    angle_deg=None,
    angle_rad=None,
    point=(0.0, 0.0, 0.0),
    axes=((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
):
    """
    Rotate structural nodes with three sequential axis-angle rotations.

    Args:
        struc_nodes (np.ndarray): Array of node positions (n_nodes, 3).
        angle_deg (array-like, optional): Three angles in degrees.
        angle_rad (array-like, optional): Three angles in radians.
        point (array-like, optional): Pivot point for rotation. Defaults to origin.
        axes (array-like, optional): Three axis definitions defining order.
            Each entry can be a label ("x"/"y"/"z") or a 3D axis vector.
            Defaults to canonical Cartesian vectors (x, y, z).

    Notes:
        - Exactly one of `angle_deg` or `angle_rad` must be provided.
        - For backward compatibility, a single scalar angle is still accepted and
          interpreted as a rotation about +Y only (legacy behavior).
    """
    if (angle_deg is None) == (angle_rad is None):
        raise ValueError("Provide exactly one of `angle_deg` or `angle_rad`.")

    # Backward compatibility with previous API that used one Y-axis angle.
    if angle_deg is not None and np.isscalar(angle_deg):
        angle_vec_rad = np.radians(np.array([0.0, float(angle_deg), 0.0]))
        axes_norm = (
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0]),
            np.array([0.0, 0.0, 1.0]),
        )
    elif angle_rad is not None and np.isscalar(angle_rad):
        angle_vec_rad = np.array([0.0, float(angle_rad), 0.0], dtype=float)
        axes_norm = (
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0]),
            np.array([0.0, 0.0, 1.0]),
        )
    else:
        if angle_deg is not None:
            angle_vec_rad = np.radians(_to_3_vector(angle_deg, "angle_deg"))
        else:
            angle_vec_rad = _to_3_vector(angle_rad, "angle_rad")

        axes_norm = _normalize_axes(axes)

    pivot = _to_3_vector(point, "point")
    nodes = np.asarray(struc_nodes, dtype=float)
    if nodes.ndim != 2 or nodes.shape[1] != 3:
        raise ValueError(
            f"`struc_nodes` must have shape (n_nodes, 3). Got shape {nodes.shape}."
        )

    rotated = nodes - pivot
    for ax, ang in zip(axes_norm, angle_vec_rad):
        R = _rotation_matrix_from_axis(ax, ang)
        rotated = rotated @ R.T

    return rotated + pivot
