"""How to run (from ASKITE repo root):

python examples/ch9/force_validation/run_ch9_3_2_vwt_cases.py \
    --cases-csv /home/jellepoland/ownCloud/phd/code/EKF-AWE/data/ch9_3_2_straight_vwt/ch9_3_2_vwt_cases_for_askite.csv \
    --mode apparent_wind_prescribed \
    --output-root results/ch9/force_validation \
    --max-iter 750

Quick wiring check:
python examples/ch9/force_validation/run_ch9_3_2_vwt_cases.py --max-cases 1 --max-iter 1
"""

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd

CH9_DIR = Path(__file__).resolve().parents[1]
if str(CH9_DIR) not in sys.path:
    sys.path.insert(0, str(CH9_DIR))

from ch9_analysis_utils import (
    DEFAULT_AERO_GEOMETRY,
    DEFAULT_CONFIG,
    DEFAULT_STRUC_GEOMETRY,
    PROJECT_DIR,
    depower_tape_length_m_to_udp,
    finite_or_nan,
    get_first,
    run_coupled_case,
    summary_from_case,
    udp_to_depower_tape_length_m,
    write_csv,
    write_json,
    write_markdown_table,
)

DEFAULT_EKF_CASES = (
    Path("/home/jellepoland/ownCloud/phd/code/EKF-AWE")
    / "data"
    / "ch9_3_2_straight_vwt"
    / "ch9_3_2_vwt_cases_for_askite.csv"
)


def build_parser():
    parser = argparse.ArgumentParser(
        description="Run ASKITE VWT simulations for EKF straight-flight validation cases."
    )
    parser.add_argument("--cases-csv", type=Path, default=DEFAULT_EKF_CASES)
    parser.add_argument(
        "--output-root",
        type=Path,
        default=PROJECT_DIR / "results" / "ch9" / "force_validation",
    )
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--struc-geometry", type=Path, default=DEFAULT_STRUC_GEOMETRY)
    parser.add_argument("--aero-geometry", type=Path, default=DEFAULT_AERO_GEOMETRY)
    parser.add_argument(
        "--mode",
        choices=["apparent_wind_prescribed", "wind_tether_state_prescribed"],
        default="apparent_wind_prescribed",
    )
    parser.add_argument("--case-id", default=None)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--continue-from-nearest", action="store_true")
    parser.add_argument("--max-cases", type=int, default=None)
    parser.add_argument("--max-iter", type=int, default=None)
    parser.add_argument("--include-gravity", action="store_true")
    parser.add_argument(
        "--include-bridle-drag",
        action="store_true",
        default=True,
        help="Deprecated: Ch. 9 validation always includes bridle drag.",
    )
    parser.add_argument("--include-tether-drag", action="store_true")
    parser.add_argument(
        "--include-kcu-drag",
        action="store_true",
        default=True,
        help="Deprecated: Ch. 9 validation always includes KCU drag.",
    )
    parser.add_argument(
        "--force-source-for-comparison",
        choices=["wing_only", "aero_total", "reaction_total"],
        default="aero_total",
    )
    return parser


def _row_value(row, names, default=np.nan):
    return get_first(row, names, default)


def _case_id(row, idx):
    value = _row_value(row, ["case_id", "id", "case"], None)
    if value is None or pd.isna(value):
        return f"case_{idx:04d}"
    return str(value)


def _target_depower(row):
    direct = _row_value(
        row,
        ["target_depower_tape_length_m", "depower_tape_length_m", "ell_d_m"],
        np.nan,
    )
    if pd.notna(direct):
        return float(direct)
    udp = _row_value(row, ["u_dp_ch9", "u_dp", "depower"], np.nan)
    if pd.notna(udp):
        return udp_to_depower_tape_length_m(udp)
    return None


def _config_overrides(row, args):
    overrides = {
        "is_with_gravity": bool(args.include_gravity),
        "is_with_aero_bridle": True,
        "is_with_aero_kcu": True,
        "is_with_aero_tether": bool(args.include_tether_drag),
    }
    if args.mode == "wind_tether_state_prescribed":
        mapping = {
            "distance_radial": ["tether_length_m", "distance_radial"],
            "speed_radial": ["reel_speed_ms", "speed_radial"],
            "angle_elevation_deg": ["elevation_deg", "angle_elevation_deg"],
            "angle_azimuth_deg": ["azimuth_deg", "angle_azimuth_deg"],
            "angle_course_deg": ["course_deg", "angle_course_deg"],
            "wind_speed_wind_ref": [
                "wind_speed_horizontal_ms",
                "wind_speed_wind_ref",
                "wind_speed_ms",
            ],
        }
        for key, names in mapping.items():
            value = _row_value(row, names, np.nan)
            if pd.notna(value):
                overrides[key] = float(value)
    return overrides


def _failure_row(case_id, case_dir, args, row, exc):
    target_depower = _target_depower(row)
    return {
        "case_id": case_id,
        "campaign": _row_value(row, ["campaign"], ""),
        "year": _row_value(row, ["year"], np.nan),
        "mode": args.mode,
        "converged": False,
        "n_iter": 0,
        "results_dir": str(case_dir),
        "requested_V_a_ms": finite_or_nan(
            _row_value(row, ["V_a_ms", "Va_ms", "va"], np.nan)
        ),
        "solved_V_a_ms": np.nan,
        "requested_u_dp_ch9": finite_or_nan(
            _row_value(row, ["u_dp_ch9", "u_dp"], np.nan)
        ),
        "target_depower_tape_length_m": finite_or_nan(target_depower),
        "target_depower_tape_u_dp": finite_or_nan(
            depower_tape_length_m_to_udp(target_depower)
            if target_depower is not None
            else np.nan
        ),
        "final_depower_tape_length_m": np.nan,
        "requested_tether_length_m": finite_or_nan(
            _row_value(row, ["tether_length_m"], np.nan)
        ),
        "requested_reel_speed_ms": finite_or_nan(
            _row_value(row, ["reel_speed_ms"], np.nan)
        ),
        "requested_elevation_deg": finite_or_nan(
            _row_value(row, ["elevation_deg"], np.nan)
        ),
        "requested_azimuth_deg": finite_or_nan(
            _row_value(row, ["azimuth_deg"], np.nan)
        ),
        "requested_course_deg": finite_or_nan(_row_value(row, ["course_deg"], np.nan)),
        "measured_tether_force_N": finite_or_nan(
            _row_value(
                row,
                [
                    "tether_force_kite_N",
                    "ground_tether_force_N",
                    "tether_force_N",
                    "force_N",
                ],
                np.nan,
            )
        ),
        "measured_CL_ekf": finite_or_nan(
            _row_value(row, ["C_L_kite", "CL_kite_ekf", "C_L", "CL", "CL_ekf"], np.nan)
        ),
        "measured_CD_ekf": finite_or_nan(
            _row_value(row, ["C_D_kite", "CD_kite_ekf", "C_D", "CD", "CD_ekf"], np.nan)
        ),
        "measured_L_over_D_ekf": finite_or_nan(
            _row_value(
                row,
                [
                    "L_over_D_kite_ekf",
                    "L_over_D_ekf",
                    "L_over_D",
                    "CL_over_CD",
                ],
                np.nan,
            )
        ),
        "measured_CL_wing_ekf": finite_or_nan(
            _row_value(row, ["C_L_wing", "CL_wing_ekf", "CL_ekf"], np.nan)
        ),
        "measured_CD_wing_ekf": finite_or_nan(
            _row_value(row, ["C_D_wing", "CD_wing_ekf", "CD_ekf"], np.nan)
        ),
        "sim_tether_or_reaction_force_N": np.nan,
        "sim_wing_force_N": np.nan,
        "sim_CL_wing": np.nan,
        "sim_CD_wing": np.nan,
        "sim_L_over_D_wing": np.nan,
        "sim_CL_total": np.nan,
        "sim_CD_total": np.nan,
        "sim_L_over_D_total": np.nan,
        "sim_CL_kite": np.nan,
        "sim_CD_kite": np.nan,
        "sim_L_over_D_kite": np.nan,
        "projected_span_m": np.nan,
        "projected_area_m2": np.nan,
        "midspan_pitch_deg": np.nan,
        "mean_twist_deg": np.nan,
        "max_abs_twist_deg": np.nan,
        "error": str(exc),
    }


def _summary_row(case_id, case_dir, args, ekf_row, meta, tracking):
    requested_va = finite_or_nan(_row_value(ekf_row, ["V_a_ms", "Va_ms", "va"], np.nan))
    udp = finite_or_nan(_row_value(ekf_row, ["u_dp_ch9", "u_dp"], np.nan))
    target_depower = _target_depower(ekf_row)
    base = summary_from_case(
        case_id,
        case_dir,
        meta,
        tracking,
        extra={
            "mode": args.mode,
            "requested_V_a_ms": requested_va,
            "requested_u_dp_ch9": udp,
            "target_depower_tape_length_m": target_depower,
            "target_depower_tape_u_dp": (
                depower_tape_length_m_to_udp(target_depower)
                if target_depower is not None
                else np.nan
            ),
        },
    )
    base.update(
        {
            "campaign": _row_value(ekf_row, ["campaign"], ""),
            "year": _row_value(ekf_row, ["year"], np.nan),
            "requested_tether_length_m": finite_or_nan(
                _row_value(ekf_row, ["tether_length_m"], np.nan)
            ),
            "requested_reel_speed_ms": finite_or_nan(
                _row_value(ekf_row, ["reel_speed_ms"], np.nan)
            ),
            "requested_elevation_deg": finite_or_nan(
                _row_value(ekf_row, ["elevation_deg"], np.nan)
            ),
            "requested_azimuth_deg": finite_or_nan(
                _row_value(ekf_row, ["azimuth_deg"], np.nan)
            ),
            "requested_course_deg": finite_or_nan(
                _row_value(ekf_row, ["course_deg"], np.nan)
            ),
            "measured_tether_force_N": finite_or_nan(
                _row_value(
                    ekf_row,
                    [
                        "tether_force_kite_N",
                        "ground_tether_force_N",
                        "tether_force_N",
                        "force_N",
                    ],
                    np.nan,
                )
            ),
            "measured_CL_ekf": finite_or_nan(
                _row_value(
                    ekf_row,
                    ["C_L_kite", "CL_kite_ekf", "C_L", "CL", "CL_ekf"],
                    np.nan,
                )
            ),
            "measured_CD_ekf": finite_or_nan(
                _row_value(
                    ekf_row,
                    ["C_D_kite", "CD_kite_ekf", "C_D", "CD", "CD_ekf"],
                    np.nan,
                )
            ),
            "measured_L_over_D_ekf": finite_or_nan(
                _row_value(
                    ekf_row,
                    [
                        "L_over_D_kite_ekf",
                        "L_over_D_ekf",
                        "L_over_D",
                        "CL_over_CD",
                    ],
                    np.nan,
                )
            ),
            "measured_CL_wing_ekf": finite_or_nan(
                _row_value(ekf_row, ["C_L_wing", "CL_wing_ekf", "CL_ekf"], np.nan)
            ),
            "measured_CD_wing_ekf": finite_or_nan(
                _row_value(ekf_row, ["C_D_wing", "CD_wing_ekf", "CD_ekf"], np.nan)
            ),
            "force_source_for_comparison": args.force_source_for_comparison,
            "include_kcu_drag_requested": bool(args.include_kcu_drag),
        }
    )
    return base


def main():
    args = build_parser().parse_args()
    if not args.cases_csv.exists():
        raise SystemExit(f"EKF case table not found yet: {args.cases_csv}")
    cases = pd.read_csv(args.cases_csv)
    if args.case_id is not None:
        cases = cases[
            cases.apply(lambda row: _case_id(row, row.name) == args.case_id, axis=1)
        ]
    if args.max_cases is not None:
        cases = cases.head(args.max_cases)

    args.output_root.mkdir(parents=True, exist_ok=True)
    rows = []
    for count, (_, ekf_row) in enumerate(cases.iterrows(), start=1):
        case_id = _case_id(ekf_row, count)
        case_dir = args.output_root / str(case_id)
        if (case_dir / "sim_output.h5").exists() and not args.force:
            try:
                from ch9_analysis_utils import load_case

                meta, tracking, _ = load_case(case_dir)
                row = _summary_row(case_id, case_dir, args, ekf_row, meta, tracking)
                rows.append(row)
                continue
            except Exception:
                pass

        try:
            requested_va = finite_or_nan(
                _row_value(ekf_row, ["V_a_ms", "Va_ms", "va"], np.nan)
            )
            if args.mode == "apparent_wind_prescribed" and not np.isfinite(
                requested_va
            ):
                raise ValueError(f"Case {case_id} has no V_a_ms column/value.")
            meta, tracking, _ = run_coupled_case(
                case_dir=case_dir,
                config_path=args.config,
                struc_geometry_path=args.struc_geometry,
                aero_geometry_path=args.aero_geometry,
                solver_mode=(
                    "level_1" if args.mode == "apparent_wind_prescribed" else "qsm"
                ),
                requested_va_ms=(
                    requested_va if args.mode == "apparent_wind_prescribed" else None
                ),
                target_depower_tape_length_m=_target_depower(ekf_row),
                steering_tape_final_extension_m=0.0,
                config_overrides=_config_overrides(ekf_row, args),
                max_iter=args.max_iter,
            )
            row = _summary_row(case_id, case_dir, args, ekf_row, meta, tracking)
            rows.append(row)
            write_json(case_dir / "case_summary.json", row)
        except Exception as exc:
            row = _failure_row(case_id, case_dir, args, ekf_row, exc)
            rows.append(row)
            write_json(case_dir / "case_summary.json", row)
        print(f"[{count}/{len(cases)}] {case_id}: converged={rows[-1]['converged']}")

    summary_csv = args.output_root / "ch9_3_2_askite_case_summary.csv"
    write_csv(summary_csv, rows)
    write_markdown_table(args.output_root / "ch9_3_2_askite_case_summary.md", rows)
    write_json(
        args.output_root / "ch9_3_2_askite_run_manifest.json",
        {
            "cases_csv": str(args.cases_csv),
            "mode": args.mode,
            "n_cases": len(cases),
            "config": str(args.config),
            "struc_geometry": str(args.struc_geometry),
            "aero_geometry": str(args.aero_geometry),
            "include_gravity": args.include_gravity,
            "include_bridle_drag": args.include_bridle_drag,
            "include_tether_drag": args.include_tether_drag,
            "include_kcu_drag_requested": args.include_kcu_drag,
            "force_source_for_comparison": args.force_source_for_comparison,
            "summary_csv": str(summary_csv),
            "continue_from_nearest": args.continue_from_nearest,
        },
    )


if __name__ == "__main__":
    main()
