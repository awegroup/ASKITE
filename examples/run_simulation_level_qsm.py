import numpy as np
from pathlib import Path
import argparse
import subprocess
import sys
import csv
import copy
from collections.abc import Mapping
from kitesim.logging_config import *
from datetime import datetime
from kitesim.utils import (
    load_and_save_config_files,
    load_yaml,
    load_sim_output,
    save_results,
    printing_rest_lengths,
    rotate_geometry,
)
from kitesim import (
    aero2struc_level_1,
    aerodynamic_vsm,
    aerostructural_coupled_solver_qsm,
    read_struc_geometry_yaml_level_1,
    structural_pss,
)
from awetrim.system.system_model import SystemModel


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run ASKITE QSM simulation")
    parser.add_argument(
        "--steering-final-extension",
        type=float,
        default=None,
        help="Override steering_tape_final_extension [m] for this run.",
    )
    parser.add_argument(
        "--steering-sweep-start",
        type=float,
        default=None,
        help="Sweep start for steering_tape_final_extension [m].",
    )
    parser.add_argument(
        "--steering-sweep-end",
        type=float,
        default=None,
        help="Sweep end for steering_tape_final_extension [m].",
    )
    parser.add_argument(
        "--steering-sweep-step",
        type=float,
        default=None,
        help="Sweep step for steering_tape_final_extension [m].",
    )
    return parser


def _run_steering_sweep(args):
    """Launch one process per steering setting to keep runs isolated."""
    if args.steering_sweep_step is None or args.steering_sweep_step <= 0:
        raise ValueError("--steering-sweep-step must be > 0")
    if args.steering_sweep_end < args.steering_sweep_start:
        raise ValueError("--steering-sweep-end must be >= --steering-sweep-start")

    values = np.arange(
        args.steering_sweep_start,
        args.steering_sweep_end + 0.5 * args.steering_sweep_step,
        args.steering_sweep_step,
    )

    script_path = Path(__file__).resolve()
    for idx, value in enumerate(values, start=1):
        print(
            f"\n=== Steering sweep {idx}/{len(values)}: steering_tape_final_extension={value:.4f} m ==="
        )
        cmd = [
            sys.executable,
            str(script_path),
            "--steering-final-extension",
            f"{float(value):.10g}",
        ]
        completed = subprocess.run(cmd, check=False)
        if completed.returncode != 0:
            raise RuntimeError(
                f"Sweep aborted at steering_tape_final_extension={value:.4f} m "
                f"(exit code {completed.returncode})."
            )


def _format_length_tag(value_m: float) -> str:
    """Format a signed length [m] into a filesystem-friendly tag."""
    sign = "p" if value_m >= 0 else "m"
    milli = int(round(abs(float(value_m)) * 1000.0))
    return f"{sign}{milli:04d}mm"


def _build_actuation_case_folder(config):
    """Build case folder name from final depower/steering actuation settings."""
    depower_tag = _format_length_tag(config.get("power_tape_final_extension", 0.0))
    steering_tag = _format_length_tag(config.get("steering_tape_final_extension", 0.0))
    return f"depower_{depower_tag}_steer_{steering_tag}"


def _resolve_starting_struc_nodes(
    config,
    project_dir,
    kite_name,
    struc_nodes_default,
):
    """
    Optionally override start nodes from a previous simulation result folder.

    Priority:
      1) config["starting_from_sim_subdir"] (new)
      2) config["starting_from_sim_of_date"] (legacy)

    The value is treated as a subdir under results/<kite_name>/, e.g.
    depower_p0100mm_steer_m0020mm/run_003.

    If both keys are empty, return struc_nodes_default.
    """
    sim_subdir = str(config.get("starting_from_sim_subdir", "")).strip()
    if sim_subdir == "":
        sim_subdir = str(config.get("starting_from_sim_of_date", "")).strip()

    if sim_subdir == "":
        return struc_nodes_default

    base_results_dir = Path(project_dir) / "results" / kite_name

    # Candidate 1: exact path from config
    candidates = [base_results_dir / sim_subdir]

    # Candidate 2/3: tolerate zero-sign naming mismatch, e.g. m0000mm vs p0000mm
    sim_subdir_m_to_p = sim_subdir.replace("m0000mm", "p0000mm")
    sim_subdir_p_to_m = sim_subdir.replace("p0000mm", "m0000mm")
    if sim_subdir_m_to_p != sim_subdir:
        candidates.append(base_results_dir / sim_subdir_m_to_p)
    if sim_subdir_p_to_m != sim_subdir:
        candidates.append(base_results_dir / sim_subdir_p_to_m)

    start_dir = None
    for cand in candidates:
        if cand.exists() and cand.is_dir():
            start_dir = cand
            break

    if start_dir is None:
        raise FileNotFoundError(
            "Configured starting simulation directory does not exist. "
            f"Tried: {', '.join(str(c) for c in candidates)}"
        )

    # Preferred: direct case-folder storage (sim_output.h5 inside start_dir).
    h5_path = start_dir / "sim_output.h5"
    # Backward compatibility: if not found, try legacy run_XXX layout.
    if not h5_path.exists():
        run_dirs = [
            d
            for d in start_dir.iterdir()
            if d.is_dir() and d.name.startswith("run_") and d.name[4:].isdigit()
        ]
        if len(run_dirs) > 0:
            start_dir = sorted(run_dirs, key=lambda p: int(p.name[4:]))[-1]
            logging.info(
                f"Using latest legacy run folder inside case folder: {start_dir.name}"
            )
            h5_path = start_dir / "sim_output.h5"

    if not h5_path.exists():
        raise FileNotFoundError(
            f"Configured starting simulation has no sim_output.h5: {h5_path}"
        )

    _, tracking_data = load_sim_output(h5_path)
    if "positions" not in tracking_data:
        raise KeyError(f"Expected 'positions' dataset in: {h5_path}")

    positions = np.asarray(tracking_data["positions"])
    if positions.ndim != 3 or positions.shape[2] != 3:
        raise ValueError(
            f"Invalid positions shape in {h5_path}: {positions.shape}. Expected (nt, n_nodes, 3)."
        )

    struc_nodes_loaded = np.array(positions[-1], dtype=float)
    if struc_nodes_loaded.shape != np.asarray(struc_nodes_default).shape:
        raise ValueError(
            "Loaded node shape does not match current geometry. "
            f"loaded={struc_nodes_loaded.shape}, current={np.asarray(struc_nodes_default).shape}"
        )

    logging.info(
        f"Starting from previous simulation final nodes: {start_dir} (n_nodes={len(struc_nodes_loaded)})"
    )
    return struc_nodes_loaded


def _resolve_starting_rest_lengths(
    config,
    project_dir,
    kite_name,
    l0_arr_default,
):
    """
    Optionally load final rest_lengths from a previous simulation result.

    If config["starting_from_sim_subdir"] is set, load the rest_lengths from the
    corresponding H5 file. Otherwise return l0_arr_default.

    Args:
        config: Configuration dictionary
        project_dir: Path to project root
        kite_name: Name of the kite
        l0_arr_default: Default rest length array from YAML geometry

    Returns:
        np.ndarray: Updated rest lengths (or defaults if not recovering)
    """
    sim_subdir = str(config.get("starting_from_sim_subdir", "")).strip()
    if sim_subdir == "":
        sim_subdir = str(config.get("starting_from_sim_of_date", "")).strip()

    if sim_subdir == "":
        return l0_arr_default

    base_results_dir = Path(project_dir) / "results" / kite_name

    # Candidate 1: exact path from config
    candidates = [base_results_dir / sim_subdir]

    # Candidate 2/3: tolerate zero-sign naming mismatch
    sim_subdir_m_to_p = sim_subdir.replace("m0000mm", "p0000mm")
    sim_subdir_p_to_m = sim_subdir.replace("p0000mm", "m0000mm")
    if sim_subdir_m_to_p != sim_subdir:
        candidates.append(base_results_dir / sim_subdir_m_to_p)
    if sim_subdir_p_to_m != sim_subdir:
        candidates.append(base_results_dir / sim_subdir_p_to_m)

    start_dir = None
    for cand in candidates:
        if cand.exists() and cand.is_dir():
            start_dir = cand
            break

    if start_dir is None:
        # No previous sim found, return defaults
        return l0_arr_default

    # Try to find sim_output.h5
    h5_path = start_dir / "sim_output.h5"
    if not h5_path.exists():
        run_dirs = [
            d
            for d in start_dir.iterdir()
            if d.is_dir() and d.name.startswith("run_") and d.name[4:].isdigit()
        ]
        if len(run_dirs) > 0:
            start_dir = sorted(run_dirs, key=lambda p: int(p.name[4:]))[-1]
            h5_path = start_dir / "sim_output.h5"

    if not h5_path.exists():
        logging.warning(
            f"No sim_output.h5 found in {start_dir}, using default rest lengths"
        )
        return l0_arr_default

    # Load rest_lengths from metadata
    try:
        metadata, _ = load_sim_output(h5_path)
        if "rest_lengths" in metadata:
            rest_lengths_loaded = np.asarray(metadata["rest_lengths"], dtype=float)
            if rest_lengths_loaded.shape == np.asarray(l0_arr_default).shape:
                logging.info(
                    f"Loaded final rest lengths from previous simulation: {start_dir.name}"
                )
                return rest_lengths_loaded
            else:
                logging.warning(
                    f"Loaded rest_lengths shape {rest_lengths_loaded.shape} "
                    f"does not match current geometry {np.asarray(l0_arr_default).shape}, "
                    f"using defaults"
                )
                return l0_arr_default
        else:
            logging.warning(
                f"No 'rest_lengths' in {h5_path} metadata, using default rest lengths"
            )
            return l0_arr_default
    except Exception as e:
        logging.warning(
            f"Error loading rest_lengths from {h5_path}: {e}, using defaults"
        )
        return l0_arr_default


def _resolve_initial_geometry_rotation_kwargs(config):
    """
    Build rotate_geometry kwargs from config while keeping legacy behavior.

    New keys (optional):
      - initial_geometry_rotation_angles_deg: [ax, ay, az]
      - initial_geometry_rotation_angles_rad: [ax, ay, az]
      - initial_geometry_rotation_point: [px, py, pz]
      - initial_geometry_rotation_axes: [a1, a2, a3]

    Legacy key fallback:
      - initial_geometry_rotation_deg (single Y-axis angle in degrees)
    """
    angle_deg = config.get("initial_geometry_rotation_angles_deg")
    angle_rad = config.get("initial_geometry_rotation_angles_rad")
    if angle_deg is not None and angle_rad is not None:
        raise ValueError(
            "Provide only one of `initial_geometry_rotation_angles_deg` or "
            "`initial_geometry_rotation_angles_rad`."
        )

    if angle_deg is None and angle_rad is None:
        angle_deg = [0.0, float(config.get("initial_geometry_rotation_deg", 0.0)), 0.0]

    return {
        "angle_deg": angle_deg,
        "angle_rad": angle_rad,
        "point": config.get("initial_geometry_rotation_point", [0.0, 0.0, 0.0]),
        "axes": config.get("initial_geometry_rotation_axes", ["x", "y", "z"]),
    }


def _configure_system_model_from_config(system_model, config):
    """Populate the AWETRIM SystemModel from config values."""
    system_model.angle_elevation = np.deg2rad(
        float(config.get("angle_elevation_deg", 30.0))
    )
    system_model.angle_azimuth = np.deg2rad(
        float(config.get("angle_azimuth_deg", 20.0))
    )
    system_model.angle_course = np.deg2rad(float(config.get("angle_course_deg", 90.0)))
    system_model.speed_radial = float(config.get("speed_radial", 0.0))
    system_model.distance_radial = float(config.get("distance_radial", 200.0))
    system_model.wind.speed_wind_ref = float(config.get("wind_speed_wind_ref", 4.0))
    system_model.timeder_speed_tangential = float(
        config.get("timeder_speed_tangential", 0.0)
    )
    system_model.timeder_speed_radial = float(config.get("timeder_speed_radial", 0.0))


def _build_qsm_csv_row(
    config,
    results,
    case_folder,
    results_dir,
    power_tape_index=None,
    steering_tape_indices=None,
):
    """
    Flatten one simulation run into a CSV row.

    Args:
        config: Configuration dictionary
        results: Metadata dictionary from solver (includes rest_lengths)
        case_folder: Case folder name
        results_dir: Results directory path
        power_tape_index: Index of power tape in rest_lengths array (optional)
        steering_tape_indices: List of [left_idx, right_idx] for steering tapes (optional)
    """
    opt_x = np.asarray(results.get("opt_x", []), dtype=float).reshape(-1)
    opt_names = [
        "kite_speed",
        "roll_deg",
        "pitch_deg",
        "yaw_deg",
        "course_rate_body",
    ]
    row = {
        "case_folder": case_folder,
        "results_dir": str(results_dir),
        "is_with_gravity": bool(config.get("is_with_gravity", False)),
        "is_with_aero_bridle": bool(config.get("is_with_aero_bridle", False)),
        "angle_elevation_deg": float(config.get("angle_elevation_deg", 30.0)),
        "angle_azimuth_deg": float(config.get("angle_azimuth_deg", 20.0)),
        "angle_course_deg": float(config.get("angle_course_deg", 90.0)),
        "speed_radial": float(config.get("speed_radial", 0.0)),
        "distance_radial": float(config.get("distance_radial", 200.0)),
        "wind_speed_wind_ref": float(config.get("wind_speed_wind_ref", 4.0)),
        "timeder_speed_tangential": float(config.get("timeder_speed_tangential", 0.0)),
        "timeder_speed_radial": float(config.get("timeder_speed_radial", 0.0)),
        "aero_roll_deg": float(results.get("aero_roll_deg", np.nan)),
        "aoa_deg": float(results.get("aoa_deg", np.nan)),
        "side_slip_deg": float(results.get("side_slip_deg", np.nan)),
    }

    # Add final actual rest_lengths instead of input extensions
    rest_lengths = np.asarray(results.get("rest_lengths", []), dtype=float)
    if power_tape_index is not None and rest_lengths.size > power_tape_index:
        row["power_tape_final_length_m"] = float(rest_lengths[power_tape_index])
    else:
        row["power_tape_final_length_m"] = np.nan

    if steering_tape_indices is not None and len(steering_tape_indices) >= 2:
        left_idx = int(steering_tape_indices[0])
        right_idx = int(steering_tape_indices[1])
        if rest_lengths.size > max(left_idx, right_idx):
            row["steering_tape_left_final_length_m"] = float(rest_lengths[left_idx])
            row["steering_tape_right_final_length_m"] = float(rest_lengths[right_idx])
        else:
            row["steering_tape_left_final_length_m"] = np.nan
            row["steering_tape_right_final_length_m"] = np.nan
    else:
        row["steering_tape_left_final_length_m"] = np.nan
        row["steering_tape_right_final_length_m"] = np.nan

    for idx, name in enumerate(opt_names):
        row[f"opt_{name}"] = float(opt_x[idx]) if idx < opt_x.size else np.nan

    return row


def _append_row_to_csv(csv_path, row):
    """Append one row to a CSV file, creating the header if needed."""
    csv_path = Path(csv_path)
    csv_path.parent.mkdir(parents=True, exist_ok=True)

    file_exists = csv_path.exists()
    fieldnames = list(row.keys())
    with csv_path.open("a", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)


# Import modules
def main():
    """Main function"""
    args = _build_arg_parser().parse_args()

    is_sweep_requested = (
        args.steering_sweep_start is not None
        or args.steering_sweep_end is not None
        or args.steering_sweep_step is not None
    )
    if is_sweep_requested:
        if (
            args.steering_sweep_start is None
            or args.steering_sweep_end is None
            or args.steering_sweep_step is None
        ):
            raise ValueError(
                "Provide all sweep args: --steering-sweep-start, --steering-sweep-end, --steering-sweep-step"
            )
        _run_steering_sweep(args)
        return

    PROJECT_DIR = Path(__file__).resolve().parents[1]
    kite_name = "TUDELFT_V3_KITE"  # the dir name with the relevant .yaml files
    # kite_name = "3plate_kite"  # the dir name with the relevant .yaml files
    # load config.yaml & geometry.yaml, save both, and return them as dicts
    config_path = Path(PROJECT_DIR) / "data" / f"{kite_name}" / "config.yaml"
    struc_geometry_path = (
        Path(PROJECT_DIR)
        / "data"
        / f"{kite_name}"
        # / "struc_geometry_level_1_manual.yaml"
        / "struc_geometry_level_1_manual_JULIA.yaml"
        # / "struc_geometry_level_1_converged.yaml"
    )
    aero_geometry_path = (
        Path(PROJECT_DIR) / "data" / f"{kite_name}" / "aero_geometry.yaml"
    )

    config_for_case = load_yaml(config_path)
    if args.steering_final_extension is not None:
        config_for_case["steering_tape_final_extension"] = float(
            args.steering_final_extension
        )

    case_folder = _build_actuation_case_folder(config=config_for_case)
    case_dir = Path(PROJECT_DIR) / "results" / f"{kite_name}" / case_folder
    results_dir = case_dir
    config, struc_geometry, aero_geometry, results_dir = load_and_save_config_files(
        config_path, struc_geometry_path, aero_geometry_path, results_dir
    )
    if args.steering_final_extension is not None:
        config["steering_tape_final_extension"] = float(args.steering_final_extension)

    logging.info(f"config files saved in {results_dir}\n")

    ###################
    ### AERODYNAMIC ###
    ###################
    n_wing_struc_nodes = len(struc_geometry["wing_particles"]["data"])
    n_struc_ribs = n_wing_struc_nodes / 2
    n_panels_aero = (n_struc_ribs - 1) * config["aerodynamic"][
        "n_aero_panels_per_struc_section"
    ]
    bridle_path = (
        struc_geometry_path if config.get("is_with_aero_bridle", False) else None
    )
    body_aero, vsm_solver, vel_app, initial_polar_data = aerodynamic_vsm.initialize(
        aero_geometry_path,
        config,
        n_panels_aero,
        bridle_path=bridle_path,
    )

    ##################
    ### STRUCTURAL ###
    ##################
    (
        # node level
        struc_nodes,
        m_arr,
        struc_node_le_indices,
        struc_node_te_indices,
        power_tape_index,
        steering_tape_indices,
        pulley_node_indices,
        # element level
        kite_connectivity_arr,
        bridle_connectivity_arr,
        bridle_diameter_arr,
        l0_arr,
        k_arr,
        c_arr,
        linktype_arr,
        pulley_line_indices,
        pulley_line_to_other_node_pair_dict,
    ) = read_struc_geometry_yaml_level_1.main(struc_geometry)

    #####################################################
    ### rotating the initial geometry by some angle,
    ### to enable the wind to be horizontal
    #####################################################
    struc_nodes = rotate_geometry(
        struc_nodes,
        **_resolve_initial_geometry_rotation_kwargs(config),
    )
    struc_nodes = _resolve_starting_struc_nodes(
        config=config,
        project_dir=PROJECT_DIR,
        kite_name=kite_name,
        struc_nodes_default=struc_nodes,
    )
    # Also recover the final rest_lengths (element l0 values) from previous simulation if available
    l0_arr = _resolve_starting_rest_lengths(
        config=config,
        project_dir=PROJECT_DIR,
        kite_name=kite_name,
        l0_arr_default=l0_arr,
    )

    # logging initial conditions
    logging.info(f"\n\nINITIAL CONDITIONS, NODES \n")
    for idx, (node_i, m_i) in enumerate(zip(struc_nodes, m_arr)):
        logging.info(f"node_idx: {idx}: node: {node_i}, mass: {m_i}")

    logging.info(f"\n\nINITIAL CONDITIONS, ELEMENTS \n")
    for idx, conn in enumerate(kite_connectivity_arr):
        logging.info(
            f"conn_idx: {idx}: conn: {conn}, l0: {l0_arr[idx]}, k: {k_arr[idx]}, c: {c_arr[idx]}, linktype: {linktype_arr[idx]}"
        )

    if config["structural_solver"] == "pss":
        ## pss -- https://github.com/awegroup/Particle_System_Simulator
        ##TODO: Fix the comment below, it SHOULD read l0
        # Note: ParticleSystem doesn’t read l0_arr. SpringDamper sets l0
        # from the initial particle positions.
        # So l0_arr is a bookkeeping array for you, not used at instantiation.
        (psystem, pss_initial_conditions, pss_params, struc_nodes_initial) = (
            structural_pss.instantiate(
                # yaml files
                config,
                # node level
                struc_nodes,
                m_arr,
                # element_level
                kite_connectivity_arr,
                l0_arr,
                k_arr,
                c_arr,
                linktype_arr,
                pulley_line_to_other_node_pair_dict,
            )
        )
        if config["is_with_initial_structure_plot"]:
            structural_pss.plot_3d_kite_structure(
                struc_nodes,
                kite_connectivity_arr,
                power_tape_index,
                k_arr=k_arr,
                c_arr=c_arr,
                linktype_arr=linktype_arr,
                pulley_nodes=pulley_node_indices,
            )
        # setting kite_fem related output to None
        kite_fem_structure = None
    else:
        raise ValueError(
            "Invalid structural solver specified, only pss is available in qsm couplded simulation for now."
        )

    ##################
    ### AERO2STRUC ###
    ##################
    aero2struc_mapping = aero2struc_level_1.initialize_mapping(
        body_aero.panels,
        struc_nodes,
        struc_node_le_indices,
        struc_node_te_indices,
    )

    #################
    ### ACTUATION ###
    #################
    initial_length_power_tape = l0_arr[power_tape_index]
    power_tape_extension_step = config["power_tape_extension_step"]
    power_tape_final_extension = config["power_tape_final_extension"]
    if power_tape_extension_step != 0:
        n_power_tape_steps = int(power_tape_final_extension / power_tape_extension_step)
    else:
        n_power_tape_steps = 0
    logging.info(f"Initial depower tape length: {l0_arr[power_tape_index]:.3f}m")
    logging.info(
        f"Desired depower tape length: {initial_length_power_tape + power_tape_final_extension:.3f}m"
    )

    initial_length_steering_left = l0_arr[steering_tape_indices[0]]
    initial_length_steering_right = l0_arr[steering_tape_indices[1]]
    steering_tape_extension_step = config["steering_tape_extension_step"]
    steering_tape_final_extension = config["steering_tape_final_extension"]
    # Apply steering if final extension != 0 (step=0 means single application, step>0 means gradual)
    if abs(float(steering_tape_final_extension)) > 1e-9:
        # Use step size from config, or final extension for single application (step=0)
        effective_step = (
            float(steering_tape_extension_step)
            if steering_tape_extension_step != 0
            else steering_tape_final_extension
        )
        aerostructural_coupled_solver_qsm.update_steering_tape_actuation(
            config=config,
            psystem=psystem,
            kite_fem_structure=kite_fem_structure,
            kite_connectivity_arr=kite_connectivity_arr,
            steering_tape_indices=steering_tape_indices,
            steering_tape_extension_step=effective_step,
            initial_length_steering_left=initial_length_steering_left,
            initial_length_steering_right=initial_length_steering_right,
            steering_tape_final_extension=steering_tape_final_extension,
        )

    ########################################
    # AWETRIM SYSTEM MODEL
    ########################################
    from awetrim.system.tether import RigidLumpedTether

    tether = RigidLumpedTether(diameter=config["tether"]["diameter"])
    system_model = SystemModel(tether=tether)
    system_model.mass_wing = float(np.sum(m_arr))
    print(
        f"Total mass of the wing (sum of particle masses): {system_model.mass_wing:.3f} kg"
    )
    _configure_system_model_from_config(system_model, config)

    ########################################
    ### AEROSTUCTURAL COUPLED SIMULATION ###
    ########################################
    tracking_data, meta = aerostructural_coupled_solver_qsm.main(
        m_arr=m_arr,
        struc_nodes=struc_nodes,
        struc_nodes_initial=struc_nodes_initial,
        system_model=system_model,
        config=config,
        ### ACTUATION
        initial_length_power_tape=initial_length_power_tape,
        n_power_tape_steps=n_power_tape_steps,
        power_tape_final_extension=power_tape_final_extension,
        power_tape_extension_step=power_tape_extension_step,
        ### CONNECTIVITY
        kite_connectivity_arr=kite_connectivity_arr,
        bridle_connectivity_arr=bridle_connectivity_arr,
        pulley_line_indices=pulley_line_indices,
        pulley_line_to_other_node_pair_dict=pulley_line_to_other_node_pair_dict,
        ### STRUC --> AERO
        struc_node_le_indices=struc_node_le_indices,
        struc_node_te_indices=struc_node_te_indices,
        ### AERO
        body_aero=copy.deepcopy(body_aero),
        vsm_solver=copy.deepcopy(vsm_solver),
        vel_app=vel_app,
        initial_polar_data=copy.deepcopy(initial_polar_data),
        bridle_diameter_arr=bridle_diameter_arr,
        ### AERO --> STRUC
        aero2struc_mapping=aero2struc_mapping,
        power_tape_index=power_tape_index,
        ### STRUC
        psystem=psystem,
        kite_fem_structure=kite_fem_structure,
    )

    # Save results
    h5_path = Path(results_dir) / "sim_output.h5"
    save_results(tracking_data, meta, h5_path)

    summary_csv_name = config.get("qsm_summary_csv_name", "qsm_summary.csv")
    summary_csv_path = Path(PROJECT_DIR) / "results" / kite_name / summary_csv_name
    summary_row = _build_qsm_csv_row(
        config=config,
        results=meta,
        case_folder=case_folder,
        results_dir=results_dir,
        power_tape_index=power_tape_index,
        steering_tape_indices=steering_tape_indices,
    )
    _append_row_to_csv(summary_csv_path, summary_row)

    # Load results
    meta_data_dict, tracking_data = load_sim_output(h5_path)

    # logging.info(f"meta_data: {meta_data_dict}")
    # - here you could add functions to plot the tracking of f_int, f_ext and f_residual over the iterations
    # - functions that make an animation of the kite going through the iterations
    # - etc.
    f_residual = tracking_data["f_int"] - tracking_data["f_ext"]

    printing_rest_lengths(tracking_data, struc_geometry)


if __name__ == "__main__":
    main()
