import copy
import json
import math
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

from kitesim import (
    aero2struc_level_1,
    aerodynamic_vsm,
    aerostructural_coupled_solver_level_1,
    aerostructural_coupled_solver_qsm,
    read_struc_geometry_yaml_level_1,
    structural_kite_fem_level_1,
    structural_pss,
)
from kitesim.analysis_metrics import compute_geometry_metrics
from kitesim.utils import load_sim_output, load_yaml, rotate_geometry, save_results


PROJECT_DIR = Path(__file__).resolve().parents[2]
KITE_NAME = "ch9"
DEFAULT_DATA_DIR = PROJECT_DIR / "data" / "ch9"
DEFAULT_CONFIG = DEFAULT_DATA_DIR / "config.yaml"
DEFAULT_STRUC_GEOMETRY = DEFAULT_DATA_DIR / "struc_geometry_PSM_reduced.yaml"
DEFAULT_AERO_GEOMETRY = DEFAULT_DATA_DIR / "aero_geometry.yaml"
UDP_DEPOWER_OFFSET_M = 0.2
UDP_DEPOWER_SCALE_M = 5.0


def udp_to_depower_tape_length_m(udp):
    return UDP_DEPOWER_OFFSET_M + UDP_DEPOWER_SCALE_M * float(udp)


def depower_tape_length_m_to_udp(length_m):
    return (float(length_m) - UDP_DEPOWER_OFFSET_M) / UDP_DEPOWER_SCALE_M


def parse_float_list(text):
    if text is None or str(text).strip() == "":
        return []
    return [float(item.strip()) for item in str(text).split(",") if item.strip()]


def parse_formats(text):
    if text is None:
        return ["pdf"]
    return [item.strip().lower() for item in str(text).split(",") if item.strip()]


def json_default(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return str(value)


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2, default=json_default)


def write_yaml(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        yaml.dump(data, fh, sort_keys=False)


def write_csv(path, rows, fieldnames=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    if fieldnames is None:
        fieldnames = []
        for row in rows:
            for key in row.keys():
                if key not in fieldnames:
                    fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_markdown_table(path, rows, fieldnames=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", encoding="utf-8") as fh:
        fh.write("| " + " | ".join(fieldnames) + " |\n")
        fh.write("| " + " | ".join(["---"] * len(fieldnames)) + " |\n")
        for row in rows:
            fh.write(
                "| "
                + " | ".join(_format_md_value(row.get(key, "")) for key in fieldnames)
                + " |\n"
            )


def _format_md_value(value):
    if value is None:
        return ""
    if isinstance(value, float):
        if not np.isfinite(value):
            return ""
        return f"{value:.5g}"
    return str(value)


def save_figure(fig, output_dir, stem, formats):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for fmt in formats:
        out = output_dir / f"{stem}.{fmt}"
        fig.savefig(out, bbox_inches="tight")
        paths.append(out)
    plt.close(fig)
    return paths


def final_valid_iteration(meta, tracking_data):
    if "positions" not in tracking_data:
        raise KeyError("tracking_data does not contain positions")
    n_rows = int(np.asarray(tracking_data["positions"]).shape[0])
    n_iter = int(meta.get("n_iter", n_rows))
    return max(0, min(n_iter - 1, n_rows - 1))


def load_case(case_dir):
    case_dir = Path(case_dir)
    meta, tracking_data = load_sim_output(case_dir / "sim_output.h5")
    return meta, tracking_data, final_valid_iteration(meta, tracking_data)


def finite_or_nan(value):
    try:
        value = float(value)
    except (TypeError, ValueError):
        return np.nan
    return value if np.isfinite(value) else np.nan


def get_first(row, names, default=np.nan):
    for name in names:
        if name in row and pd.notna(row[name]):
            return row[name]
    return default


def infer_column(df, candidates, required=False):
    for candidate in candidates:
        if candidate in df.columns:
            return candidate
    if required:
        raise KeyError(f"None of the required columns are present: {candidates}")
    return None


def _resolve_initial_geometry_rotation_kwargs(config):
    angle_deg = config.get("initial_geometry_rotation_angles_deg")
    angle_rad = config.get("initial_geometry_rotation_angles_rad")
    if angle_deg is not None and angle_rad is not None:
        raise ValueError(
            "Provide only one of initial_geometry_rotation_angles_deg or "
            "initial_geometry_rotation_angles_rad."
        )
    if angle_deg is None and angle_rad is None:
        angle_deg = [0.0, float(config.get("initial_geometry_rotation_deg", 0.0)), 0.0]
    return {
        "angle_deg": angle_deg,
        "angle_rad": angle_rad,
        "point": config.get("initial_geometry_rotation_point", [0.0, 0.0, 0.0]),
        "axes": config.get("initial_geometry_rotation_axes", ["x", "y", "z"]),
    }


def _configure_system_model(system_model, config):
    system_model.angle_elevation = np.deg2rad(float(config.get("angle_elevation_deg", 0.0)))
    system_model.angle_azimuth = np.deg2rad(float(config.get("angle_azimuth_deg", 0.0)))
    system_model.angle_course = np.deg2rad(float(config.get("angle_course_deg", 90.0)))
    system_model.speed_radial = float(config.get("speed_radial", 0.0))
    system_model.distance_radial = float(config.get("distance_radial", 200.0))
    system_model.wind.speed_wind_ref = float(config.get("wind_speed_wind_ref", 6.0))
    system_model.timeder_speed_tangential = float(
        config.get("timeder_speed_tangential", 0.0)
    )
    system_model.timeder_speed_radial = float(config.get("timeder_speed_radial", 0.0))


def _instantiate_structure(
    config,
    struc_geometry,
    struc_nodes,
    m_arr,
    kite_connectivity_arr,
    l0_arr,
    k_arr,
    c_arr,
    linktype_arr,
    pulley_line_to_other_node_pair_dict,
    power_tape_index,
    pulley_node_indices,
):
    if config["structural_solver"] == "pss":
        psystem, _, _, struc_nodes_initial = structural_pss.instantiate(
            config,
            struc_nodes,
            m_arr,
            kite_connectivity_arr,
            l0_arr,
            k_arr,
            c_arr,
            linktype_arr,
            pulley_line_to_other_node_pair_dict,
        )
        kite_fem_structure = None
    elif config["structural_solver"] == "kite_fem":
        (
            kite_fem_structure,
            _,
            _,
            _,
            struc_nodes_initial,
        ) = structural_kite_fem_level_1.instantiate(
            config,
            struc_geometry,
            struc_nodes,
            kite_connectivity_arr,
            l0_arr,
            k_arr,
            c_arr,
            m_arr,
            linktype_arr,
            pulley_line_to_other_node_pair_dict,
        )
        struc_nodes = struc_nodes_initial.copy()
        psystem = None
    else:
        raise ValueError("structural_solver must be pss or kite_fem")

    if config.get("is_with_initial_structure_plot", False):
        structural_pss.plot_3d_kite_structure(
            struc_nodes,
            kite_connectivity_arr,
            power_tape_index,
            k_arr=k_arr,
            c_arr=c_arr,
            linktype_arr=linktype_arr,
            pulley_nodes=pulley_node_indices,
        )
    return psystem, kite_fem_structure, struc_nodes, struc_nodes_initial


def run_coupled_case(
    case_dir,
    config_path=DEFAULT_CONFIG,
    struc_geometry_path=DEFAULT_STRUC_GEOMETRY,
    aero_geometry_path=DEFAULT_AERO_GEOMETRY,
    solver_mode="level_1",
    requested_va_ms=None,
    target_depower_tape_length_m=None,
    power_tape_final_extension_m=None,
    steering_tape_final_extension_m=0.0,
    config_overrides=None,
    max_iter=None,
):
    """
    Run one ASKITE case and save config/geometry, sim_output.h5, and JSON inputs.
    """
    case_dir = Path(case_dir)
    case_dir.mkdir(parents=True, exist_ok=True)
    config = copy.deepcopy(load_yaml(Path(config_path)))
    struc_geometry = load_yaml(Path(struc_geometry_path))
    aero_geometry = load_yaml(Path(aero_geometry_path))

    config["starting_from_sim_subdir"] = ""
    config["starting_from_sim_of_date"] = ""
    config["steering_tape_final_extension"] = float(steering_tape_final_extension_m)
    config["steering_tape_extension_step"] = 0.0
    config["is_with_initial_structure_plot"] = False
    config["is_with_coupling_plot_per_iteration"] = False
    config["is_with_struc_plot_per_iteration"] = False
    config["is_with_aero_plot_per_iteration"] = False
    config["is_with_final_plot"] = False
    config.setdefault("aero_structural_solver", {})
    config["aero_structural_solver"]["aero_solve_type"] = solver_mode
    if max_iter is not None:
        config["aero_structural_solver"]["max_iter"] = int(max_iter)
    if requested_va_ms is not None:
        config["wind_speed_wind_ref"] = float(requested_va_ms)
        config["vel_wind"] = [float(requested_va_ms), 0.0, 0.0]
        config["vel_kite"] = [0.0, 0.0, 0.0]
        config["acc_kite"] = [0.0, 0.0, 0.0]

    if config_overrides:
        for key, value in config_overrides.items():
            if "." in key:
                first, rest = key.split(".", 1)
                config.setdefault(first, {})[rest] = value
            else:
                config[key] = value

    (
        struc_nodes,
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
    ) = read_struc_geometry_yaml_level_1.main(struc_geometry, config=config)

    initial_length_power_tape = float(l0_arr[power_tape_index])
    if target_depower_tape_length_m is not None:
        config["power_tape_final_extension"] = (
            float(target_depower_tape_length_m) - initial_length_power_tape
        )
    elif power_tape_final_extension_m is not None:
        config["power_tape_final_extension"] = float(power_tape_final_extension_m)
    config.setdefault("power_tape_extension_step", 0.0)

    write_yaml(case_dir / "config.yaml", config)
    write_yaml(case_dir / "struc_geometry.yaml", struc_geometry)
    write_yaml(case_dir / "aero_geometry.yaml", aero_geometry)

    struc_nodes = rotate_geometry(
        struc_nodes,
        **_resolve_initial_geometry_rotation_kwargs(config),
    )

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
    if requested_va_ms is not None:
        vel_app = np.array([float(requested_va_ms), 0.0, 0.0])

    psystem, kite_fem_structure, struc_nodes, struc_nodes_initial = _instantiate_structure(
        config,
        struc_geometry,
        struc_nodes,
        m_arr,
        kite_connectivity_arr,
        l0_arr,
        k_arr,
        c_arr,
        linktype_arr,
        pulley_line_to_other_node_pair_dict,
        power_tape_index,
        pulley_node_indices,
    )

    aero2struc_mapping = aero2struc_level_1.initialize_mapping(
        body_aero.panels,
        struc_nodes,
        struc_node_le_indices,
        struc_node_te_indices,
    )

    power_tape_extension_step = float(config.get("power_tape_extension_step", 0.0))
    power_tape_final_extension = float(config.get("power_tape_final_extension", 0.0))
    n_power_tape_steps = (
        int(np.ceil(abs(power_tape_final_extension) / abs(power_tape_extension_step)))
        if abs(power_tape_extension_step) > 1e-12
        else 0
    )

    case_inputs = {
        "solver_mode": solver_mode,
        "requested_va_ms": requested_va_ms,
        "target_depower_tape_length_m": target_depower_tape_length_m,
        "initial_length_power_tape_m": initial_length_power_tape,
        "initial_u_dp": depower_tape_length_m_to_udp(initial_length_power_tape),
        "power_tape_final_extension_m": power_tape_final_extension,
        "steering_tape_final_extension_m": steering_tape_final_extension_m,
        "config_path": str(config_path),
        "struc_geometry_path": str(struc_geometry_path),
        "aero_geometry_path": str(aero_geometry_path),
    }
    write_json(case_dir / "case_inputs.json", case_inputs)

    if solver_mode == "level_1":
        tracking_data, meta = aerostructural_coupled_solver_level_1.main(
            m_arr=m_arr,
            struc_nodes=struc_nodes,
            struc_nodes_initial=struc_nodes_initial,
            config=config,
            initial_length_power_tape=initial_length_power_tape,
            n_power_tape_steps=n_power_tape_steps,
            power_tape_final_extension=power_tape_final_extension,
            power_tape_extension_step=power_tape_extension_step,
            kite_connectivity_arr=kite_connectivity_arr,
            bridle_connectivity_arr=bridle_connectivity_arr,
            pulley_line_indices=pulley_line_indices,
            pulley_line_to_other_node_pair_dict=pulley_line_to_other_node_pair_dict,
            struc_node_le_indices=struc_node_le_indices,
            struc_node_te_indices=struc_node_te_indices,
            body_aero=body_aero,
            vsm_solver=vsm_solver,
            vel_app=vel_app,
            initial_polar_data=initial_polar_data,
            bridle_diameter_arr=bridle_diameter_arr,
            aero2struc_mapping=aero2struc_mapping,
            power_tape_index=power_tape_index,
            psystem=psystem,
            kite_fem_structure=kite_fem_structure,
        )
    elif solver_mode == "qsm":
        from awetrim.system.system_model import SystemModel
        from awetrim.system.tether import RigidLumpedTether

        tether = RigidLumpedTether(diameter=config.get("tether", {}).get("diameter", 0.01))
        system_model = SystemModel(tether=tether)
        system_model.mass_wing = float(np.sum(m_arr))
        _configure_system_model(system_model, config)
        initial_length_steering_left = float(l0_arr[steering_tape_indices[0]])
        initial_length_steering_right = float(l0_arr[steering_tape_indices[1]])

        tracking_data, meta = aerostructural_coupled_solver_qsm.main(
            m_arr=m_arr,
            struc_nodes=struc_nodes,
            struc_nodes_initial=struc_nodes_initial,
            system_model=system_model,
            config=config,
            initial_length_power_tape=initial_length_power_tape,
            n_power_tape_steps=n_power_tape_steps,
            power_tape_final_extension=power_tape_final_extension,
            power_tape_extension_step=power_tape_extension_step,
            initial_length_steering_left=initial_length_steering_left,
            initial_length_steering_right=initial_length_steering_right,
            steering_tape_indices=steering_tape_indices,
            steering_tape_final_extension=float(steering_tape_final_extension_m),
            steering_tape_extension_step=float(config.get("steering_tape_extension_step", 0.0)),
            kite_connectivity_arr=kite_connectivity_arr,
            bridle_connectivity_arr=bridle_connectivity_arr,
            pulley_line_indices=pulley_line_indices,
            pulley_line_to_other_node_pair_dict=pulley_line_to_other_node_pair_dict,
            struc_node_le_indices=struc_node_le_indices,
            struc_node_te_indices=struc_node_te_indices,
            body_aero=copy.deepcopy(body_aero),
            vsm_solver=copy.deepcopy(vsm_solver),
            vel_app=vel_app,
            initial_polar_data=copy.deepcopy(initial_polar_data),
            bridle_diameter_arr=bridle_diameter_arr,
            aero2struc_mapping=aero2struc_mapping,
            power_tape_index=power_tape_index,
            psystem=psystem,
            kite_fem_structure=kite_fem_structure,
        )
    else:
        raise ValueError("solver_mode must be level_1 or qsm")

    h5_path = case_dir / "sim_output.h5"
    save_results(tracking_data, meta, h5_path)
    return meta, tracking_data, case_inputs


def summary_from_case(
    case_id,
    case_dir,
    meta,
    tracking_data,
    extra=None,
):
    final_idx = final_valid_iteration(meta, tracking_data)
    extra = extra or {}
    wing_force = (
        np.asarray(tracking_data.get("aero_force_wing_total", [[np.nan] * 3]))[
            final_idx
        ]
        if "aero_force_wing_total" in tracking_data
        else np.full(3, np.nan)
    )
    total_force = (
        np.asarray(tracking_data.get("aero_force_total", [[np.nan] * 3]))[final_idx]
        if "aero_force_total" in tracking_data
        else wing_force
    )
    row = {
        "case_id": case_id,
        "mode": extra.get("mode", "level_1"),
        "converged": bool(meta.get("converged", False)),
        "n_iter": int(meta.get("n_iter", final_idx + 1)),
        "results_dir": str(case_dir),
        "requested_V_a_ms": finite_or_nan(extra.get("requested_V_a_ms", np.nan)),
        "solved_V_a_ms": finite_or_nan(meta.get("final_V_a", meta.get("va", np.nan))),
        "requested_u_dp_ch9": finite_or_nan(extra.get("requested_u_dp_ch9", np.nan)),
        "initial_depower_tape_length_m": finite_or_nan(
            extra.get("initial_length_power_tape_m", np.nan)
        ),
        "initial_u_dp": finite_or_nan(extra.get("initial_u_dp", np.nan)),
        "target_depower_tape_length_m": finite_or_nan(
            extra.get("target_depower_tape_length_m", np.nan)
        ),
        "final_depower_tape_length_m": finite_or_nan(
            meta.get("final_depower_tape_length_m", np.nan)
        ),
        "power_tape_final_extension_m": finite_or_nan(
            meta.get("final_power_tape_extension_m", np.nan)
        ),
        "final_residual_norm": finite_or_nan(meta.get("final_residual_norm", np.nan)),
        "sim_tether_or_reaction_force_N": finite_or_nan(
            meta.get("tether_force", np.linalg.norm(total_force))
        ),
        "sim_wing_force_N": finite_or_nan(np.linalg.norm(wing_force)),
        "sim_aero_total_force_N": finite_or_nan(np.linalg.norm(total_force)),
        "sim_CL_wing": finite_or_nan(meta.get("final_C_L_wing", np.nan)),
        "sim_CD_wing": finite_or_nan(meta.get("final_C_D_wing", np.nan)),
        "sim_L_over_D_wing": finite_or_nan(meta.get("final_glide_ratio_wing", np.nan)),
        "sim_CL_total": finite_or_nan(meta.get("final_C_L_total_aero", np.nan)),
        "sim_CD_total": finite_or_nan(meta.get("final_C_D_total_aero", np.nan)),
        "sim_L_over_D_total": finite_or_nan(
            meta.get("final_glide_ratio_total_aero", np.nan)
        ),
        "projected_span_m": finite_or_nan(meta.get("final_projected_span_m", np.nan)),
        "projected_area_m2": finite_or_nan(meta.get("final_projected_area_m2", np.nan)),
        "midspan_pitch_deg": finite_or_nan(meta.get("final_pitch_or_trim_deg", np.nan)),
        "mean_twist_deg": finite_or_nan(meta.get("final_mean_twist_deg", np.nan)),
        "max_abs_twist_deg": finite_or_nan(meta.get("final_max_abs_twist_deg", np.nan)),
    }
    row.update(extra)
    return row


def metric_file_for_shape(case_dir, meta, tracking_data):
    final_idx = final_valid_iteration(meta, tracking_data)
    nodes = np.asarray(tracking_data["positions"])[final_idx]
    le = np.asarray(meta.get("struc_node_le_indices", []), dtype=int)
    te = np.asarray(meta.get("struc_node_te_indices", []), dtype=int)
    return compute_geometry_metrics(nodes, le, te)


def equalize_2d_axes(ax, x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size == 0 or y.size == 0:
        return
    xmin, xmax = np.nanmin(x), np.nanmax(x)
    ymin, ymax = np.nanmin(y), np.nanmax(y)
    span = max(xmax - xmin, ymax - ymin, 1e-6)
    xmid = 0.5 * (xmin + xmax)
    ymid = 0.5 * (ymin + ymax)
    ax.set_xlim(xmid - 0.55 * span, xmid + 0.55 * span)
    ax.set_ylim(ymid - 0.55 * span, ymid + 0.55 * span)
    ax.set_aspect("equal", adjustable="box")
