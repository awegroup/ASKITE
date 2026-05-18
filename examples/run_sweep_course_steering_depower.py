"""
3-D parameter sweep: course angle × steering tape extension × depower tape extension.

Configurable sweep parameters: easily adjust COURSE_ANGLES_DEG, STEERING_*, and DEPOWER_*
at the top of this file to change sweep scope without rewriting code.

All combinations are run in-process and every result row is appended to a CSV inside
results/<kite_name>/. The CSV filename includes the sweep counts (e.g.,
sweep_c1_s6_d1.csv for 1 course × 6 steering × 1 depower = 6 total cases).

The CSV includes:
  - course angle, steering, depower actuation
  - angle of attack (aoa_deg), side-slip (side_slip_deg)
  - lift and drag coefficients (cl, cd)
  - quasi-steady optimised state (kite speed, roll, pitch, yaw, course rate)
  - convergence flag

Examples:
  - Single depower only:  Set DEPOWER_START_M = DEPOWER_END_M = 0.0
  - 5-point steering sweep: Set STEERING_START_M = 0.0, STEERING_END_M = 0.20, STEERING_STEP_M = 0.05
  - Single course angle: Set COURSE_ANGLES_DEG = [90.0]
"""

import copy
import csv
import logging
from pathlib import Path

import numpy as np

from kitesim.logging_config import *  # noqa: F401,F403  (sets up root logger)
from kitesim.utils import (
    load_and_save_config_files,
    load_yaml,
    load_sim_output,
    printing_rest_lengths,
    rotate_geometry,
    save_results,
)
from kitesim import (
    aero2struc_level_1,
    aerodynamic_vsm,
    aerostructural_coupled_solver_qsm,
    read_struc_geometry_yaml_level_1,
    structural_pss,
)
from awetrim.system.system_model import SystemModel
from awetrim.system.tether import RigidLumpedTether

# ── Sweep parameters ─────────────────────────────────────────────────────────
# Easily configure each parameter:
#   - COURSE_ANGLES_DEG: List of course angles (degrees)
#   - STEERING: start, end, step (meters)
#   - DEPOWER: start, end, step (meters)

# Example configurations:
#   Single depower, multiple steering: set DEPOWER_START_M == DEPOWER_END_M
#   Single steering, multiple depower: set STEERING_START_M == STEERING_END_M
#   Single course: COURSE_ANGLES_DEG = [90.0]

COURSE_ANGLES_DEG: list = [90.0]  # Course angle in degrees
STEERING_START_M: float = 0.0  # m
STEERING_END_M: float = 0.5  # m
STEERING_STEP_M: float = 0.05  # m
DEPOWER_START_M: float = 0.0  # m|
DEPOWER_END_M: float = 0.8  # m (only depower=0.0)
DEPOWER_STEP_M: float = 0.2  # m (step size; ignored if START==END)

# Calculate total cases (will be shown at runtime)
_n_course = len(COURSE_ANGLES_DEG)
_steering_values = np.arange(
    STEERING_START_M, STEERING_END_M + 0.5 * STEERING_STEP_M, STEERING_STEP_M
)
_depower_values = (
    np.arange(DEPOWER_START_M, DEPOWER_END_M + 0.5 * DEPOWER_STEP_M, DEPOWER_STEP_M)
    if DEPOWER_STEP_M > 0
    else [DEPOWER_START_M]
)
_n_steering = len(_steering_values)
_n_depower = len(_depower_values)
_total_cases = _n_course * _n_steering * _n_depower

# CSV filename with sweep parameters (editable for custom names)
SUMMARY_CSV_NAME: str = f"sweep_c{_n_course}_s{_n_steering}_d{_n_depower}.csv"

# Fixed CSV fieldnames (order matters for alignment)
CSV_FIELDNAMES: list = [
    "wind_speed_wind_ref_ms",
    "depower_tape_final_extension_m",
    "steering_tape_final_extension_m",
    "case_folder",
    "results_dir",
    "is_with_gravity",
    "is_with_aero_bridle",
    "angle_elevation_deg",
    "angle_azimuth_deg",
    "angle_course_deg",
    "speed_radial",
    "distance_radial",
    "timeder_speed_tangential",
    "timeder_speed_radial",
    "aoa_deg",
    "side_slip_deg",
    "aero_roll_deg",
    "cl",
    "cd",
    "converged",
    "depower_tape_final_length_m",
    "steering_tape_left_final_length_m",
    "steering_tape_right_final_length_m",
    "steering_tape_avg_length_m",
    "steering_tape_asymmetry_m",
    "opt_kite_speed",
    "opt_roll_deg",
    "opt_pitch_deg",
    "opt_yaw_deg",
    "opt_course_rate_body",
    "va",
    "tether_force",
]
# ─────────────────────────────────────────────────────────────────────────────


# ── Helpers (self-contained copies) ───────────────────────────────────────────


def _format_length_tag(value_m: float) -> str:
    """Format a signed length [m] as a filesystem-friendly tag, e.g. p0150mm."""
    sign = "p" if value_m >= 0 else "m"
    milli = int(round(abs(float(value_m)) * 1000.0))
    return f"{sign}{milli:04d}mm"


def _parse_length_tag_mm(tag: str) -> float:
    """Parse signed mm tag like p0200mm/m0050mm into meters."""
    clean = str(tag).strip()
    if len(clean) < 4 or not clean.endswith("mm"):
        raise ValueError(f"Invalid length tag format: {tag}")

    sign_char = clean[0]
    if sign_char == "p":
        sign = 1.0
    elif sign_char == "m":
        sign = -1.0
    else:
        raise ValueError(f"Invalid length tag sign in: {tag}")

    milli = int(clean[1:-2])
    return sign * (milli / 1000.0)


def _infer_source_extensions_from_subdir(sim_subdir: str) -> tuple[float, float]:
    """Infer (depower, steering) source extensions [m] from folder naming."""
    if str(sim_subdir).strip() == "":
        return 0.0, 0.0

    name = str(sim_subdir).replace("\\", "/").strip("/").split("/")[-1]
    depower_marker = "depower_"
    steer_marker = "_steer_"

    if depower_marker not in name or steer_marker not in name:
        logging.warning(
            "Could not infer source extensions from starting_from_sim_subdir='%s'; assuming 0.0m.",
            sim_subdir,
        )
        return 0.0, 0.0

    depower_start = name.index(depower_marker) + len(depower_marker)
    steer_start = name.index(steer_marker, depower_start)
    depower_tag = name[depower_start:steer_start]
    steering_tag = name[steer_start + len(steer_marker) :]

    try:
        return _parse_length_tag_mm(depower_tag), _parse_length_tag_mm(steering_tag)
    except Exception:
        logging.warning(
            "Failed parsing source extensions from starting folder '%s'; assuming 0.0m.",
            name,
        )
        return 0.0, 0.0


def _build_actuation_case_folder(config: dict) -> str:
    depower_tag = _format_length_tag(config.get("power_tape_final_extension", 0.0))
    steering_tag = _format_length_tag(config.get("steering_tape_final_extension", 0.0))
    course_tag = f"course_{int(config.get('angle_course_deg', 90)):03d}deg"
    return f"{course_tag}_depower_{depower_tag}_steer_{steering_tag}"


def _resolve_initial_geometry_rotation_kwargs(config: dict) -> dict:
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


def _configure_system_model(system_model: SystemModel, config: dict) -> None:
    """Populate a SystemModel from config values, including wind speed."""
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


def _append_csv_row(csv_path: Path, row: dict) -> None:
    """Append *row* to *csv_path*, writing the header on first write.

    Uses fixed CSV_FIELDNAMES to ensure consistent column alignment across all rows.
    """
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    file_exists = csv_path.exists()
    with csv_path.open("a", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDNAMES)
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)


def _build_csv_row(
    wind_speed: float,
    steering: float,
    config: dict,
    meta: dict,
    case_folder: str,
    results_dir: Path,
    power_tape_index: int = None,
    steering_tape_indices: list = None,
) -> dict:
    """Flatten one simulation result into a CSV row."""
    opt_x = np.asarray(meta.get("opt_x", []), dtype=float).reshape(-1)
    opt_names = ["kite_speed", "roll_deg", "pitch_deg", "yaw_deg", "course_rate_body"]

    row = {
        "wind_speed_wind_ref_ms": wind_speed,
        "depower_tape_final_extension_m": float(
            config.get("power_tape_final_extension", 0.0)
        ),
        "steering_tape_final_extension_m": float(steering),
        "case_folder": case_folder,
        "results_dir": str(results_dir),
        "is_with_gravity": bool(config.get("is_with_gravity", False)),
        "is_with_aero_bridle": bool(config.get("is_with_aero_bridle", False)),
        "angle_elevation_deg": float(config.get("angle_elevation_deg", 30.0)),
        "angle_azimuth_deg": float(config.get("angle_azimuth_deg", 20.0)),
        "angle_course_deg": float(config.get("angle_course_deg", 90.0)),
        "speed_radial": float(config.get("speed_radial", 0.0)),
        "distance_radial": float(config.get("distance_radial", 200.0)),
        "timeder_speed_tangential": float(config.get("timeder_speed_tangential", 0.0)),
        "timeder_speed_radial": float(config.get("timeder_speed_radial", 0.0)),
        # Aerodynamic solution
        "aoa_deg": float(meta.get("aoa_deg", np.nan)),
        "side_slip_deg": float(meta.get("side_slip_deg", np.nan)),
        "aero_roll_deg": float(meta.get("aero_roll_deg", np.nan)),
        "cl": float(meta.get("cl", np.nan)),
        "cd": float(meta.get("cd", np.nan)),
        "converged": bool(meta.get("converged", False)),
    }

    # Add final actual rest_lengths
    rest_lengths = np.asarray(meta.get("rest_lengths", []), dtype=float)
    if power_tape_index is not None and rest_lengths.size > power_tape_index:
        row["depower_tape_final_length_m"] = float(rest_lengths[power_tape_index])
    else:
        row["depower_tape_final_length_m"] = np.nan

    if steering_tape_indices is not None and len(steering_tape_indices) >= 2:
        left_idx = int(steering_tape_indices[0])
        right_idx = int(steering_tape_indices[1])
        if rest_lengths.size > max(left_idx, right_idx):
            left_length = float(rest_lengths[left_idx])
            right_length = float(rest_lengths[right_idx])
            row["steering_tape_left_final_length_m"] = left_length
            row["steering_tape_right_final_length_m"] = right_length
            row["steering_tape_avg_length_m"] = (left_length + right_length) / 2.0
            row["steering_tape_asymmetry_m"] = (left_length - right_length) / 2.0
        else:
            row["steering_tape_left_final_length_m"] = np.nan
            row["steering_tape_right_final_length_m"] = np.nan
            row["steering_tape_avg_length_m"] = np.nan
            row["steering_tape_asymmetry_m"] = np.nan
    else:
        row["steering_tape_left_final_length_m"] = np.nan
        row["steering_tape_right_final_length_m"] = np.nan
        row["steering_tape_avg_length_m"] = np.nan
        row["steering_tape_asymmetry_m"] = np.nan

    for idx, name in enumerate(opt_names):
        row[f"opt_{name}"] = float(opt_x[idx]) if idx < opt_x.size else np.nan

    row["va"] = float(meta.get("va", np.nan))
    row["tether_force"] = float(meta.get("tether_force", np.nan))

    return row


def _resolve_starting_struc_nodes(case_dir: Path, l0_arr_default):
    """Load final struc_nodes from a previous simulation in this case folder."""
    h5_path = case_dir / "sim_output.h5"
    if not h5_path.exists():
        return None

    try:
        _, tracking_data = load_sim_output(h5_path)
        if "positions" not in tracking_data:
            return None
        positions = np.asarray(tracking_data["positions"])
        if positions.ndim != 3 or positions.shape[2] != 3:
            return None
        struc_nodes_loaded = np.array(positions[-1], dtype=float)
        logging.info(f"Recovered struc_nodes from previous run: {h5_path.parent.name}")
        return struc_nodes_loaded
    except Exception as e:
        logging.warning(f"Could not load struc_nodes from {h5_path}: {e}")
        return None


def _resolve_starting_rest_lengths(case_dir: Path, l0_arr_default):
    """Load final rest_lengths from a previous simulation in this case folder."""
    h5_path = case_dir / "sim_output.h5"
    if not h5_path.exists():
        return l0_arr_default

    try:
        metadata, _ = load_sim_output(h5_path)
        if "rest_lengths" in metadata:
            rest_lengths_loaded = np.asarray(metadata["rest_lengths"], dtype=float)
            if rest_lengths_loaded.shape == np.asarray(l0_arr_default).shape:
                logging.info(
                    f"Recovered rest_lengths from previous run: {h5_path.parent.name}"
                )
                return rest_lengths_loaded
        return l0_arr_default
    except Exception as e:
        logging.warning(f"Could not load rest_lengths from {h5_path}: {e}")
        return l0_arr_default


def _select_nearest_case(
    target_course: float,
    target_steering: float,
    target_depower: float,
    simulated_cases: list,
    course_scale: float,
    steering_scale: float,
    depower_scale: float,
):
    """Return nearest simulated case in normalized (course, steering, depower) space."""
    if not simulated_cases:
        return None

    def _dist2(case):
        dc = (float(target_course) - float(case["course"])) / course_scale
        ds = (float(target_steering) - float(case["steering"])) / steering_scale
        dd = (float(target_depower) - float(case["depower"])) / depower_scale
        return dc * dc + ds * ds + dd * dd

    return min(simulated_cases, key=_dist2)


# ── Main sweep ────────────────────────────────────────────────────────────────


def main() -> None:
    PROJECT_DIR = Path(__file__).resolve().parents[1]
    kite_name = "TUDELFT_V3_KITE"

    config_path = PROJECT_DIR / "data" / kite_name / "config.yaml"
    struc_geometry_path = (
        PROJECT_DIR / "data" / kite_name / "struc_geometry_level_1_manual_JULIA.yaml"
    )
    aero_geometry_path = PROJECT_DIR / "data" / kite_name / "aero_geometry.yaml"

    # Build parameter ranges
    course_values = COURSE_ANGLES_DEG
    steering_values = np.arange(
        STEERING_START_M,
        STEERING_END_M + 0.5 * STEERING_STEP_M,
        STEERING_STEP_M,
    )
    # Handle depower: if START == END, use single value; otherwise sweep
    if DEPOWER_START_M == DEPOWER_END_M:
        depower_values = np.array([DEPOWER_START_M])
    else:
        depower_values = np.arange(
            DEPOWER_START_M,
            DEPOWER_END_M + 0.5 * DEPOWER_STEP_M,
            DEPOWER_STEP_M,
        )

    summary_csv_path = PROJECT_DIR / "results" / kite_name / SUMMARY_CSV_NAME

    total_runs = len(course_values) * len(steering_values) * len(depower_values)

    # Log sweep configuration
    logging.info("\n" + "=" * 80)
    logging.info("SWEEP CONFIGURATION")
    logging.info("=" * 80)
    logging.info(f"Course angles:    {course_values} ({len(course_values)} values)")
    logging.info(
        f"Steering range:   {STEERING_START_M:.4f} to {STEERING_END_M:.4f} m (step: {STEERING_STEP_M:.4f} m) → {len(steering_values)} values"
    )
    logging.info(
        f"Depower range:    {DEPOWER_START_M:.4f} to {DEPOWER_END_M:.4f} m (step: {DEPOWER_STEP_M:.4f} m) → {len(depower_values)} values"
    )
    logging.info(f"Total runs:       {total_runs}")
    logging.info(f"Output CSV:       {summary_csv_path}")
    logging.info("=" * 80 + "\n")

    run_idx = 0

    # ── One-time: load base config and geometry ───────────────────────────────
    base_config = load_yaml(config_path)
    struc_geometry = load_yaml(struc_geometry_path)

    # Structural arrays from geometry
    (
        struc_nodes_base,
        m_arr,
        struc_node_le_indices,
        struc_node_te_indices,
        power_tape_index,
        steering_tape_indices,
        pulley_node_indices,
        kite_connectivity_arr,
        bridle_connectivity_arr,
        bridle_diameter_arr,
        l0_arr,
        k_arr,
        c_arr,
        linktype_arr,
        pulley_line_indices,
        pulley_line_to_other_node_pair_dict,
    ) = read_struc_geometry_yaml_level_1.main(struc_geometry, config=base_config)

    # Apply initial geometry rotation once
    struc_nodes_base = rotate_geometry(
        struc_nodes_base,
        **_resolve_initial_geometry_rotation_kwargs(base_config),
    )

    # ── One-time: initialise aerodynamic solver ───────────────────────────────
    n_wing_struc_nodes = len(struc_geometry["wing_particles"]["data"])
    n_struc_ribs = n_wing_struc_nodes / 2
    n_panels_aero = (n_struc_ribs - 1) * base_config["aerodynamic"][
        "n_aero_panels_per_struc_section"
    ]
    bridle_path = (
        struc_geometry_path if base_config.get("is_with_aero_bridle", False) else None
    )
    body_aero_init, vsm_solver_init, vel_app, initial_polar_data_init = (
        aerodynamic_vsm.initialize(
            aero_geometry_path,
            base_config,
            n_panels_aero,
            bridle_path=bridle_path,
        )
    )

    # Extract wind speed from config (single wind speed for entire sweep)
    wind_speed = base_config.get("wind_speed_wind_ref", 6.0)

    # Resolve optional fixed restart source from config (same behavior as single-run).
    start_subdir = str(base_config.get("starting_from_sim_subdir", "")).strip()
    if start_subdir == "":
        start_subdir = str(base_config.get("starting_from_sim_of_date", "")).strip()

    source_depower_from_start, source_steering_from_start = (
        _infer_source_extensions_from_subdir(start_subdir)
    )

    start_case_dir = None
    if start_subdir != "":
        base_results_dir = PROJECT_DIR / "results" / kite_name
        start_subdir_m_to_p = start_subdir.replace("m0000mm", "p0000mm")
        start_subdir_p_to_m = start_subdir.replace("p0000mm", "m0000mm")
        candidates = [base_results_dir / start_subdir]
        if start_subdir_m_to_p != start_subdir:
            candidates.append(base_results_dir / start_subdir_m_to_p)
        if start_subdir_p_to_m != start_subdir:
            candidates.append(base_results_dir / start_subdir_p_to_m)

        for cand in candidates:
            if cand.exists() and cand.is_dir():
                start_case_dir = cand
                break

        if start_case_dir is None:
            raise FileNotFoundError(
                "Configured starting simulation directory does not exist for sweep. "
                f"Tried: {', '.join(str(c) for c in candidates)}"
            )
        logging.info(
            "Sweep initialization mode: always start from configured simulation folder: %s",
            start_case_dir,
        )
        logging.info(
            "Source extensions from configured start folder: depower=%.4fm, steering=%.4fm",
            source_depower_from_start,
            source_steering_from_start,
        )
    else:
        logging.info(
            "Sweep initialization mode: always start from baseline initial state"
        )

    # ── Loop: THREE nested loops for course × steering × depower ──────────────
    for course_angle in course_values:
        for steering in steering_values:
            for depower in depower_values:
                run_idx += 1
                logging.info(
                    "\n=== Sweep run %d/%d: "
                    "course=%.1f deg  steering=%.4f m  depower=%.4f m  wind=%.2f m/s ===",
                    run_idx,
                    total_runs,
                    course_angle,
                    steering,
                    depower,
                    wind_speed,
                )

                # Build per-run config with current sweep parameters
                cfg = copy.deepcopy(base_config)
                cfg["wind_speed_wind_ref"] = wind_speed
                cfg["angle_course_deg"] = float(course_angle)
                cfg["steering_tape_final_extension"] = float(steering)
                cfg["power_tape_final_extension"] = float(depower)
                # Keep depower actuation step from config.
                # DEPOWER_STEP_M defines sweep grid spacing, not in-simulation actuation step.
                cfg["power_tape_extension_step"] = float(
                    base_config.get("power_tape_extension_step", 0.0)
                )

                # Determine output directory for this run
                case_folder = _build_actuation_case_folder(cfg)
                case_dir = PROJECT_DIR / "results" / kite_name / case_folder
                case_dir.mkdir(parents=True, exist_ok=True)

                # Persist config snapshot
                _, _, _, results_dir = load_and_save_config_files(
                    config_path, struc_geometry_path, aero_geometry_path, case_dir
                )

                # Always initialize from the fixed config-selected source.
                source_steering = float(source_steering_from_start)
                source_depower = float(source_depower_from_start)

                l0_arr_active = l0_arr
                if start_case_dir is not None:
                    struc_nodes_recovered = _resolve_starting_struc_nodes(
                        start_case_dir, struc_nodes_base
                    )
                    l0_arr_active = _resolve_starting_rest_lengths(
                        start_case_dir, l0_arr
                    )
                    if struc_nodes_recovered is None:
                        raise FileNotFoundError(
                            f"Configured start case has no valid sim_output positions: {start_case_dir}"
                        )
                    struc_nodes = struc_nodes_recovered.copy()
                    logging.info(
                        "Starting from configured simulation state: %s",
                        start_case_dir,
                    )
                else:
                    struc_nodes = struc_nodes_base.copy()
                    logging.info(
                        "Starting from baseline initial geometry and rest lengths"
                    )

                initial_length_power_tape = l0_arr_active[power_tape_index]

                (psystem, pss_initial_conditions, pss_params, struc_nodes_initial) = (
                    structural_pss.instantiate(
                        cfg,
                        struc_nodes,
                        m_arr,
                        kite_connectivity_arr,
                        l0_arr_active,
                        k_arr,
                        c_arr,
                        linktype_arr,
                        pulley_line_to_other_node_pair_dict,
                    )
                )

                if cfg["is_with_initial_structure_plot"]:
                    structural_pss.plot_3d_kite_structure(
                        struc_nodes,
                        kite_connectivity_arr,
                        power_tape_index,
                        k_arr=k_arr,
                        c_arr=c_arr,
                        linktype_arr=linktype_arr,
                        pulley_nodes=pulley_node_indices,
                    )

                # Apply steering actuation
                steering_tape_extension_step = cfg.get(
                    "steering_tape_extension_step", 0.0
                )
                steering_tape_final_extension = cfg.get(
                    "steering_tape_final_extension", 0.0
                )

                steering_to_apply = (
                    float(steering_tape_final_extension) - source_steering
                )
                logging.info(
                    "Steering delta to apply: target %.4fm - source %.4fm = %.4fm",
                    float(steering_tape_final_extension),
                    source_steering,
                    float(steering_to_apply),
                )
                initial_length_steering_left = float(
                    l0_arr_active[steering_tape_indices[0]]
                )
                initial_length_steering_right = float(
                    l0_arr_active[steering_tape_indices[1]]
                )

                # Power-tape actuation
                power_tape_extension_step = cfg.get("power_tape_extension_step", 0.0)
                power_tape_final_extension = cfg.get("power_tape_final_extension", 0.0)

                depower_to_apply = float(power_tape_final_extension) - source_depower
                desired_length_power_tape = float(initial_length_power_tape) + float(
                    depower_to_apply
                )
                logging.info(
                    "Depower delta to apply: target %.4fm - source %.4fm = %.4fm | length %.3fm -> %.3fm",
                    float(power_tape_final_extension),
                    source_depower,
                    float(depower_to_apply),
                    float(initial_length_power_tape),
                    float(desired_length_power_tape),
                )

                # Two-level depower behavior:
                # 1) Sweep target (this file): DEPOWER_START/END/STEP define case grid.
                # 2) In-simulation progression (config): power_tape_extension_step ramps toward target.
                if (
                    abs(float(depower_to_apply)) > 1e-9
                    and abs(float(power_tape_extension_step)) <= 1e-12
                ):
                    logging.warning(
                        "Depower target is non-zero (%.4fm) but power_tape_extension_step in config is 0.0; "
                        "internal progressive depower will not move.",
                        float(depower_to_apply),
                    )

                n_power_tape_steps = (
                    int(np.ceil(abs(depower_to_apply) / abs(power_tape_extension_step)))
                    if power_tape_extension_step != 0
                    else 0
                )

                # ── SystemModel ───────────────────────────────────────────────────────
                from awetrim.system.tether import Tether

                tether = RigidLumpedTether(diameter=cfg["tether"]["diameter"])
                system_model = SystemModel(tether=tether)
                system_model.mass_wing = float(np.sum(m_arr))
                _configure_system_model(system_model, cfg)

                # ── Aero–structure mapping (created per-run) ───────────────────────────
                aero2struc_mapping = aero2struc_level_1.initialize_mapping(
                    body_aero_init.panels,
                    struc_nodes,
                    struc_node_le_indices,
                    struc_node_te_indices,
                )

                # ── Run the coupled solver ─────────────────────────────────────────────
                converged = False
                tracking_data = None
                try:
                    tracking_data, meta = aerostructural_coupled_solver_qsm.main(
                        m_arr=m_arr,
                        struc_nodes=struc_nodes,
                        struc_nodes_initial=struc_nodes_initial,
                        system_model=system_model,
                        config=cfg,
                        # Actuation
                        initial_length_power_tape=initial_length_power_tape,
                        n_power_tape_steps=n_power_tape_steps,
                        power_tape_final_extension=depower_to_apply,
                        power_tape_extension_step=power_tape_extension_step,
                        initial_length_steering_left=initial_length_steering_left,
                        initial_length_steering_right=initial_length_steering_right,
                        steering_tape_indices=steering_tape_indices,
                        steering_tape_final_extension=steering_to_apply,
                        steering_tape_extension_step=steering_tape_extension_step,
                        # Connectivity
                        kite_connectivity_arr=kite_connectivity_arr,
                        bridle_connectivity_arr=bridle_connectivity_arr,
                        pulley_line_indices=pulley_line_indices,
                        pulley_line_to_other_node_pair_dict=pulley_line_to_other_node_pair_dict,
                        # Struc → Aero
                        struc_node_le_indices=struc_node_le_indices,
                        struc_node_te_indices=struc_node_te_indices,
                        # Aero
                        body_aero=copy.deepcopy(body_aero_init),
                        vsm_solver=copy.deepcopy(vsm_solver_init),
                        vel_app=vel_app,
                        initial_polar_data=copy.deepcopy(initial_polar_data_init),
                        bridle_diameter_arr=bridle_diameter_arr,
                        # Aero → Struc
                        aero2struc_mapping=aero2struc_mapping,
                        power_tape_index=power_tape_index,
                        # Struc
                        psystem=psystem,
                        kite_fem_structure=None,
                    )
                    converged = True
                except Exception as e:
                    logging.error(f"Solver failed: {e}")
                    meta = {}

                # Save run state for post-processing.
                if tracking_data is not None:
                    h5_path = Path(results_dir) / "sim_output.h5"
                    save_results(tracking_data, meta, h5_path)

                # ── Append result to CSV ───────────────────────────────────────────────
                csv_row = _build_csv_row(
                    wind_speed=wind_speed,
                    steering=steering,
                    config=cfg,
                    meta=meta,
                    case_folder=case_folder,
                    results_dir=results_dir,
                    power_tape_index=power_tape_index,
                    steering_tape_indices=steering_tape_indices,
                )
                _append_csv_row(summary_csv_path, csv_row)
                logging.info(f"Appended to CSV: {summary_csv_path}")

    logging.info(f"\n=== Sweep complete! Total runs: {total_runs} ===")
    logging.info(f"Summary CSV: {summary_csv_path}")


if __name__ == "__main__":
    main()
