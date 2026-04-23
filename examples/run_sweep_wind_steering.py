"""
2-D parameter sweep: wind reference speed × steering tape extension.

Wind speeds are defined in WIND_SPEEDS_MS (m/s).
Steering range is swept from STEERING_START_M to STEERING_END_M in steps of
STEERING_STEP_M.  All combinations are run in-process and every result row
is appended to SUMMARY_CSV_NAME inside results/<kite_name>/.

The CSV includes:
  - wind_speed_wind_ref, steering, depower actuation
  - angle of attack (aoa_deg), side-slip (side_slip_deg)
  - lift and drag coefficients (cl, cd)
  - quasi-steady optimised state (kite speed, roll, pitch, yaw, course rate)
  - convergence flag

NOTE on AOA calculation:
  The AOA in results is `aoa_deg` = atan2(v_z, v_x) in the course frame, where:
    - v_x is apparent velocity along course direction (horizontal)
    - v_z is apparent velocity along radial direction (vertical/lift)
  This gives the pitch angle of the apparent wind relative to the kite.
  This value accounts for the kite's trim attitude rotations and includes
  induced velocity effects from the VSM solver.
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
# Wind speed is taken from config file (wind_speed_wind_ref)

STEERING_START_M: float = 0.0  # m
STEERING_END_M: float = 0.25  # m
STEERING_STEP_M: float = 0.05  # m

# Depower tape extension applied to every run in the sweep.
# Set DEPOWER_STEP_M = 0 to jump in a single step (recommended for sweep use).
DEPOWER_M: float = 0.0  # m  — final power-tape extension
DEPOWER_STEP_M: float = 0.0  # m  — actuation step size (0 = single step)

SUMMARY_CSV_NAME: str = "sweep_wind_steering.csv"

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


# ── Helpers (self-contained copies so this file is independent) ───────────────


def _format_length_tag(value_m: float) -> str:
    """Format a signed length [m] as a filesystem-friendly tag, e.g. p0150mm."""
    sign = "p" if value_m >= 0 else "m"
    milli = int(round(abs(float(value_m)) * 1000.0))
    return f"{sign}{milli:04d}mm"


def _build_actuation_case_folder(config: dict) -> str:
    depower_tag = _format_length_tag(config.get("power_tape_final_extension", 0.0))
    steering_tag = _format_length_tag(config.get("steering_tape_final_extension", 0.0))
    return f"depower_{depower_tag}_steer_{steering_tag}"


def _parse_steering_from_case_folder(case_dir_or_name: str) -> float:
    """
    Extract steering amount (in meters) from case folder name.
    E.g., "depower_p0000mm_steer_p0150mm" -> 0.150 m
    """
    case_name = str(case_dir_or_name).split("/")[-1]  # Handle both Path and string
    if "steer_" not in case_name:
        return 0.0
    try:
        steer_part = case_name.split("steer_")[1]  # e.g., "p0150mm"
        sign_char = steer_part[0]  # 'p' or 'm'
        milli_str = steer_part[1:5]  # e.g., "0150"
        milli = int(milli_str)
        value_m = milli / 1000.0
        if sign_char == "m":
            value_m = -value_m
        return value_m
    except (IndexError, ValueError):
        return 0.0


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
    """
    Flatten one simulation result into a CSV row.

    Args:
        wind_speed: Wind speed in m/s
        steering: Steering extension in m
        config: Configuration dictionary
        meta: Metadata dictionary from solver (includes rest_lengths)
        case_folder: Case folder name
        results_dir: Results directory path
        power_tape_index: Index of power tape in rest_lengths array (optional)
        steering_tape_indices: List of [left_idx, right_idx] for steering tapes (optional)
    """
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

    # Add final actual rest_lengths instead of input extensions
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
            # Add computed columns: average length and asymmetry
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

    # Apparent wind speed magnitude
    row["va"] = float(meta.get("va", np.nan))

    # Tether force: resultant force in radial (Z) direction
    # This is computed from aero, gravity, and inertial forces at equilibrium
    row["tether_force"] = float(meta.get("tether_force", np.nan))

    return row


def _resolve_starting_struc_nodes(
    case_dir: Path,
    l0_arr_default,
):
    """
    Load final struc_nodes from a previous simulation in this case folder.
    If not found, return default.
    """
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


def _resolve_starting_rest_lengths(
    case_dir: Path,
    l0_arr_default,
):
    """
    Load final rest_lengths from a previous simulation in this case folder.
    If not found, return defaults.
    """
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


# ── Main sweep ────────────────────────────────────────────────────────────────


def main() -> None:
    PROJECT_DIR = Path(__file__).resolve().parents[1]
    kite_name = "TUDELFT_V3_KITE"

    config_path = PROJECT_DIR / "data" / kite_name / "config.yaml"
    struc_geometry_path = (
        PROJECT_DIR / "data" / kite_name / "struc_geometry_level_1_manual_JULIA.yaml"
    )
    aero_geometry_path = PROJECT_DIR / "data" / kite_name / "aero_geometry.yaml"

    # Build steering values
    steering_values = np.arange(
        STEERING_START_M,
        STEERING_END_M + 0.5 * STEERING_STEP_M,
        STEERING_STEP_M,
    )

    summary_csv_path = PROJECT_DIR / "results" / kite_name / SUMMARY_CSV_NAME

    total_runs = len(steering_values)
    run_idx = 0

    # ── One-time: load base config and geometry ───────────────────────────────
    base_config = load_yaml(config_path)
    struc_geometry = load_yaml(struc_geometry_path)

    # Structural arrays from geometry (shared – steering actuation modifies
    # spring lengths inside psystem at runtime, not these arrays directly)
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

    # ── Aero–structure mapping (created per-run with current geometry) ────────────
    # Note: Must be created inside loop to match run_simulation_level_qsm.py behavior
    # which creates it after geometry setup/recovery

    # ── Track previous case for recovery sequence ──────────────────────────────
    previous_case_dir = None

    # Extract wind speed from config (single wind speed per sweep)
    wind_speed = base_config.get("wind_speed_wind_ref", 6.0)

    # ── Loop: steering values ─────────────────────────────────────────────────
    for steering in steering_values:
        run_idx += 1
        logging.info(
            "\n=== Sweep run %d/%d: wind=%.2f m/s  depower=%.4f m  steering=%.4f m ===",
            run_idx,
            total_runs,
            wind_speed,
            DEPOWER_M,
            steering,
        )

        # Build a per-run config copy with the current wind/steering/depower values
        cfg = copy.deepcopy(base_config)
        cfg["wind_speed_wind_ref"] = wind_speed
        cfg["steering_tape_final_extension"] = float(steering)
        cfg["power_tape_final_extension"] = float(DEPOWER_M)
        # Step size: use DEPOWER_STEP_M if > 0, otherwise a single step equal
        # to the full extension (avoids divide-by-zero in the coupled solver).
        cfg["power_tape_extension_step"] = (
            float(DEPOWER_STEP_M)
            if DEPOWER_STEP_M > 0
            else float(DEPOWER_M) if DEPOWER_M != 0 else 0.0
        )

        # Determine output directory for this run
        case_folder = _build_actuation_case_folder(cfg)
        case_dir = PROJECT_DIR / "results" / kite_name / case_folder
        case_dir.mkdir(parents=True, exist_ok=True)

        # Persist config snapshot alongside the run results
        _, _, _, results_dir = load_and_save_config_files(
            config_path, struc_geometry_path, aero_geometry_path, case_dir
        )

        # Recover final struc_nodes and rest_lengths from PREVIOUS run in sequence
        # (not from current case folder; only after first run)
        if run_idx > 1 and previous_case_dir is not None:
            struc_nodes_recovered = _resolve_starting_struc_nodes(
                previous_case_dir, struc_nodes_base
            )
            l0_arr_active = _resolve_starting_rest_lengths(previous_case_dir, l0_arr)
        else:
            # First run: use defaults
            struc_nodes_recovered = None
            l0_arr_active = l0_arr
        initial_length_power_tape = l0_arr_active[power_tape_index]

        # Use recovered nodes if available, else start fresh
        struc_nodes = (
            struc_nodes_recovered.copy()
            if struc_nodes_recovered is not None
            else struc_nodes_base.copy()
        )
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
        steering_tape_extension_step = cfg.get("steering_tape_extension_step", 0.0)
        steering_tape_final_extension = cfg.get("steering_tape_final_extension", 0.0)

        # When recovering from a previous run, apply INCREMENTAL steering
        # (not the absolute value from config)
        steering_to_apply = steering_tape_final_extension
        if run_idx > 1 and previous_case_dir is not None:
            previous_steering = _parse_steering_from_case_folder(previous_case_dir.name)
            steering_to_apply = steering_tape_final_extension - previous_steering
            logging.info(
                f"Steering adjustment: current={steering_tape_final_extension:.4f}m - "
                f"previous={previous_steering:.4f}m = incremental={steering_to_apply:.4f}m"
            )

        # Apply steering if adjustment != 0 (step=0 means single application, step>0 means gradual)
        if abs(float(steering_to_apply)) > 1e-9:
            # Use step size from config, or incremental extension for single application (step=0)
            effective_step = (
                float(steering_tape_extension_step)
                if steering_tape_extension_step != 0
                else steering_to_apply
            )
            aerostructural_coupled_solver_qsm.update_steering_tape_actuation(
                config=cfg,
                psystem=psystem,
                kite_fem_structure=None,
                kite_connectivity_arr=kite_connectivity_arr,
                steering_tape_indices=steering_tape_indices,
                steering_tape_extension_step=effective_step,
                initial_length_steering_left=l0_arr_active[steering_tape_indices[0]],
                initial_length_steering_right=l0_arr_active[steering_tape_indices[1]],
                steering_tape_final_extension=steering_to_apply,
            )

        # Power-tape actuation
        power_tape_extension_step = cfg.get("power_tape_extension_step", 0.0)
        power_tape_final_extension = cfg.get("power_tape_final_extension", 0.0)
        n_power_tape_steps = (
            int(power_tape_final_extension / power_tape_extension_step)
            if power_tape_extension_step != 0
            else 0
        )

        # ── SystemModel ───────────────────────────────────────────────────
        from awetrim.system.tether import Tether

        tether = RigidLumpedTether(diameter=cfg["tether"]["diameter"])
        system_model = SystemModel(tether=tether)
        system_model.mass_wing = float(np.sum(m_arr))
        _configure_system_model(system_model, cfg)

        # ── Aero–structure mapping (created per-run with current geometry) ────────
        aero2struc_mapping = aero2struc_level_1.initialize_mapping(
            body_aero_init.panels,
            struc_nodes,  # Use current struc_nodes (may be recovered)
            struc_node_le_indices,
            struc_node_te_indices,
        )

        # ── Run coupled solver ────────────────────────────────────────────
        # Deepcopy aero objects so each run starts from a clean initial state.
        tracking_data, meta = aerostructural_coupled_solver_qsm.main(
            m_arr=m_arr,
            struc_nodes=struc_nodes,
            struc_nodes_initial=struc_nodes_initial,
            system_model=system_model,
            config=cfg,
            # Actuation
            initial_length_power_tape=initial_length_power_tape,
            n_power_tape_steps=n_power_tape_steps,
            power_tape_final_extension=power_tape_final_extension,
            power_tape_extension_step=power_tape_extension_step,
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

        # Save simulation output
        h5_path = results_dir / "sim_output.h5"
        save_results(tracking_data, meta, h5_path)

        # Structural warm-start: PSS particle positions (pre-rotation) for next step.
        final_nodes = meta.get("final_struc_nodes")
        if final_nodes is not None and np.all(np.isfinite(final_nodes)):
            struc_nodes_warm = np.asarray(final_nodes, dtype=float)

        # QSM warm-start: carry converged opt_x to the next steering step.
        opt_x = np.asarray(meta.get("opt_x", []), dtype=float)

        logging.info(
            "Run complete: wind=%.2f m/s  steering=%.4f m  opt_x=%s",
            wind_speed,
            steering,
            opt_x,
        )

        # Append to sweep summary CSV
        row = _build_csv_row(
            wind_speed=wind_speed,
            steering=float(steering),
            config=cfg,
            meta=meta,
            case_folder=case_folder,
            results_dir=results_dir,
            power_tape_index=power_tape_index,
            steering_tape_indices=steering_tape_indices,
        )
        _append_csv_row(summary_csv_path, row)

        logging.info(
            "Run %d/%d done | depower=%.4f m  steering=%.4f m | "
            "aoa=%.2f deg  side_slip=%.2f deg  converged=%s",
            run_idx,
            total_runs,
            DEPOWER_M,
            float(steering),
            row["aoa_deg"],
            row["side_slip_deg"],
            row["converged"],
        )

        # Track this case for next run's recovery
        previous_case_dir = case_dir

    logging.info(
        "\nSweep complete (wind=%.2f m/s). Summary CSV: %s",
        wind_speed,
        summary_csv_path,
    )


if __name__ == "__main__":
    main()
