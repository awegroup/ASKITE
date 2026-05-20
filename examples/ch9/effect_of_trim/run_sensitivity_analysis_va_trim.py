import argparse
from pathlib import Path
import sys

import numpy as np

CH9_DIR = Path(__file__).resolve().parents[1]
if str(CH9_DIR) not in sys.path:
    sys.path.insert(0, str(CH9_DIR))

from ch9_analysis_utils import (
    DEFAULT_AERO_GEOMETRY,
    DEFAULT_CONFIG,
    DEFAULT_STRUC_GEOMETRY,
    PROJECT_DIR,
    parse_float_list,
    depower_tape_length_m_to_udp,
    run_coupled_case,
    summary_from_case,
    udp_to_depower_tape_length_m,
    write_csv,
    write_json,
    write_markdown_table,
)

# V3 kite velocities during 2019-2025 flight are 10m to 25m/s
# so 17.5 is a good middle way, for the sweep
DEFAULT_VA_VALUES = "17.5"
# DEFAULT_UDP_VALUES = "0.22,0.25,0.28,0.31,0.33,0.36, 0.39"
DEFAULT_UDP_VALUES = "0.18,0.25,0.3,0.35,0.41"


def build_parser():
    parser = argparse.ArgumentParser(
        description="Run controlled Section 9.3.1 V_a and depower sensitivity cases."
    )
    parser.add_argument("--sweep", choices=["va", "depower", "both"], default="both")
    parser.add_argument(
        "--output-root",
        type=Path,
        default=PROJECT_DIR / "results" / "ch9" / "effect_of_trim",
    )
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--struc-geometry", type=Path, default=DEFAULT_STRUC_GEOMETRY)
    parser.add_argument("--aero-geometry", type=Path, default=DEFAULT_AERO_GEOMETRY)
    parser.add_argument("--solver-mode", choices=["level_1", "qsm"], default="level_1")
    parser.add_argument("--va-values", default=DEFAULT_VA_VALUES)
    parser.add_argument("--udp-values", default=DEFAULT_UDP_VALUES)
    parser.add_argument("--baseline-udp", type=float, default=0.35)
    parser.add_argument("--baseline-va", type=float, default=14.0)
    parser.add_argument("--continue-from-previous", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--max-cases", type=int, default=None)
    parser.add_argument("--max-iter", type=int, default=None)
    return parser


def _case_folder(sweep_type, va, udp):
    va_tag = f"va_{int(round(float(va) * 10)):04d}"
    udp_tag = f"udp_{int(round(float(udp) * 1000)):04d}"
    return f"{sweep_type}_{va_tag}_{udp_tag}"


def _failure_row(case_id, case_dir, sweep_type, va, udp, exc):
    target_depower_tape_length_m = udp_to_depower_tape_length_m(udp)
    return {
        "case_id": case_id,
        "sweep_type": sweep_type,
        "requested_va_ms": float(va),
        "solved_va_ms": np.nan,
        "requested_u_dp": float(udp),
        "target_depower_tape_length_m": target_depower_tape_length_m,
        "target_depower_tape_u_dp": depower_tape_length_m_to_udp(
            target_depower_tape_length_m
        ),
        "final_depower_tape_length_m": np.nan,
        "power_tape_final_extension_m": np.nan,
        "converged": False,
        "n_iter": 0,
        "final_residual_norm": np.nan,
        "C_L_wing": np.nan,
        "C_D_wing": np.nan,
        "glide_ratio_wing": np.nan,
        "resultant_force_N": np.nan,
        "tether_or_reaction_force_N": np.nan,
        "projected_span_m": np.nan,
        "projected_area_m2": np.nan,
        "midspan_pitch_deg": np.nan,
        "mean_twist_deg": np.nan,
        "max_abs_twist_deg": np.nan,
        "results_dir": str(case_dir),
        "error": str(exc),
    }


def _summary_row_for_sensitivity(base_row, sweep_type, va, udp):
    target_depower_tape_length_m = udp_to_depower_tape_length_m(udp)
    return {
        "case_id": base_row["case_id"],
        "sweep_type": sweep_type,
        "requested_va_ms": float(va),
        "solved_va_ms": base_row["solved_V_a_ms"],
        "requested_u_dp": float(udp),
        "initial_depower_tape_length_m": base_row.get(
            "initial_depower_tape_length_m", np.nan
        ),
        "initial_u_dp": base_row.get("initial_u_dp", np.nan),
        "target_depower_tape_length_m": target_depower_tape_length_m,
        "target_depower_tape_u_dp": depower_tape_length_m_to_udp(
            target_depower_tape_length_m
        ),
        "final_depower_tape_length_m": base_row["final_depower_tape_length_m"],
        "power_tape_final_extension_m": base_row["power_tape_final_extension_m"],
        "converged": base_row["converged"],
        "n_iter": base_row["n_iter"],
        "final_residual_norm": base_row["final_residual_norm"],
        "C_L_wing": base_row["sim_CL_wing"],
        "C_D_wing": base_row["sim_CD_wing"],
        "glide_ratio_wing": base_row["sim_L_over_D_wing"],
        "resultant_force_N": base_row["sim_aero_total_force_N"],
        "tether_or_reaction_force_N": base_row["sim_tether_or_reaction_force_N"],
        "projected_span_m": base_row["projected_span_m"],
        "projected_area_m2": base_row["projected_area_m2"],
        "midspan_pitch_deg": base_row["midspan_pitch_deg"],
        "mean_twist_deg": base_row["mean_twist_deg"],
        "max_abs_twist_deg": base_row["max_abs_twist_deg"],
        "results_dir": base_row["results_dir"],
        "error": "",
    }


def main():
    args = build_parser().parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)
    va_values = parse_float_list(args.va_values)
    udp_values = parse_float_list(args.udp_values)

    requested_cases = []
    if args.sweep in ("va", "both"):
        for va in va_values:
            requested_cases.append(("va", va, args.baseline_udp))
    if args.sweep in ("depower", "both"):
        for udp in udp_values:
            requested_cases.append(("depower", args.baseline_va, udp))
    if args.max_cases is not None:
        requested_cases = requested_cases[: args.max_cases]

    rows = []
    for idx, (sweep_type, va, udp) in enumerate(requested_cases, start=1):
        target_depower_tape_length_m = udp_to_depower_tape_length_m(udp)
        case_id = _case_folder(sweep_type, va, udp)
        case_dir = args.output_root / case_id
        if (case_dir / "sim_output.h5").exists() and not args.force:
            try:
                from ch9_analysis_utils import load_case

                meta, tracking, _ = load_case(case_dir)
                base = summary_from_case(
                    case_id,
                    case_dir,
                    meta,
                    tracking,
                    extra={
                        "mode": args.solver_mode,
                        "requested_V_a_ms": va,
                        "requested_u_dp_ch9": udp,
                        "initial_length_power_tape_m": np.nan,
                        "initial_u_dp": np.nan,
                        "target_depower_tape_length_m": target_depower_tape_length_m,
                    },
                )
                rows.append(_summary_row_for_sensitivity(base, sweep_type, va, udp))
                continue
            except Exception:
                pass

        try:
            meta, tracking, case_inputs = run_coupled_case(
                case_dir=case_dir,
                config_path=args.config,
                struc_geometry_path=args.struc_geometry,
                aero_geometry_path=args.aero_geometry,
                solver_mode=args.solver_mode,
                requested_va_ms=va,
                target_depower_tape_length_m=target_depower_tape_length_m,
                steering_tape_final_extension_m=0.0,
                max_iter=args.max_iter,
            )
            base = summary_from_case(
                case_id,
                case_dir,
                meta,
                tracking,
                extra={
                    "mode": args.solver_mode,
                    "requested_V_a_ms": va,
                    "requested_u_dp_ch9": udp,
                    "initial_length_power_tape_m": (
                        case_inputs.get("initial_length_power_tape_m", np.nan)
                        if isinstance(case_inputs, dict)
                        else np.nan
                    ),
                    "initial_u_dp": (
                        case_inputs.get("initial_u_dp", np.nan)
                        if isinstance(case_inputs, dict)
                        else np.nan
                    ),
                    "target_depower_tape_length_m": target_depower_tape_length_m,
                },
            )
            rows.append(_summary_row_for_sensitivity(base, sweep_type, va, udp))
            write_json(case_dir / "case_summary.json", rows[-1])
        except Exception as exc:
            row = _failure_row(case_id, case_dir, sweep_type, va, udp, exc)
            rows.append(row)
            write_json(case_dir / "case_summary.json", row)
        print(
            f"[{idx}/{len(requested_cases)}] {case_id}: converged={rows[-1]['converged']}"
        )

    csv_path = args.output_root / "sensitivity_va_trim_summary.csv"
    write_csv(csv_path, rows)
    write_markdown_table(args.output_root / "sensitivity_va_trim_summary.md", rows)
    write_json(
        args.output_root / "run_manifest.json",
        {
            "sweep": args.sweep,
            "solver_mode": args.solver_mode,
            "n_cases": len(requested_cases),
            "config": str(args.config),
            "struc_geometry": str(args.struc_geometry),
            "aero_geometry": str(args.aero_geometry),
            "continue_from_previous": args.continue_from_previous,
            "summary_csv": str(csv_path),
        },
    )


if __name__ == "__main__":
    main()
