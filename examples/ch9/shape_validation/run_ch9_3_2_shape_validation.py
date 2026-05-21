"""Run the single Chapter 9 shape-validation ASKITE case.

Run from the ASKITE repo root:

python examples/ch9/shape_validation/run_ch9_3_2_shape_validation.py \
  --max-iter 750
"""

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
    depower_tape_length_m_to_udp,
    reset_processed_data_dir,
    run_coupled_case,
    summary_from_case,
    udp_to_depower_tape_length_m,
    write_csv,
    write_json,
    write_markdown_table,
)

DEFAULT_VA_MS = 16.75
DEFAULT_UDP = 0.4151


def build_parser():
    parser = argparse.ArgumentParser(
        description="Run the single ASKITE case used for Ch. 9 shape validation."
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=PROJECT_DIR / "results" / "ch9" / "shape_validation",
    )
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--struc-geometry", type=Path, default=DEFAULT_STRUC_GEOMETRY)
    parser.add_argument("--aero-geometry", type=Path, default=DEFAULT_AERO_GEOMETRY)
    parser.add_argument("--solver-mode", choices=["level_1", "qsm"], default="level_1")
    parser.add_argument("--va", type=float, default=DEFAULT_VA_MS)
    parser.add_argument("--udp", type=float, default=DEFAULT_UDP)
    parser.add_argument("--max-iter", type=int, default=None)
    return parser


def _case_folder(va, udp):
    va_tag = f"va_{int(round(float(va) * 100)):04d}"
    udp_tag = f"udp_{int(round(float(udp) * 10000)):05d}"
    return f"shape_validation_{va_tag}_{udp_tag}"


def _failure_row(case_id, case_dir, args, exc):
    target_depower_tape_length_m = udp_to_depower_tape_length_m(args.udp)
    return {
        "case_id": case_id,
        "mode": args.solver_mode,
        "converged": False,
        "n_iter": 0,
        "results_dir": str(case_dir),
        "requested_V_a_ms": float(args.va),
        "requested_u_dp_ch9": float(args.udp),
        "target_depower_tape_length_m": target_depower_tape_length_m,
        "target_depower_tape_u_dp": depower_tape_length_m_to_udp(
            target_depower_tape_length_m
        ),
        "solved_V_a_ms": np.nan,
        "final_depower_tape_length_m": np.nan,
        "final_residual_norm": np.nan,
        "error": str(exc),
    }


def main():
    args = build_parser().parse_args()
    processed_root = reset_processed_data_dir(args.output_root)
    case_id = _case_folder(args.va, args.udp)
    case_dir = processed_root / case_id
    target_depower_tape_length_m = udp_to_depower_tape_length_m(args.udp)

    try:
        meta, tracking, case_inputs = run_coupled_case(
            case_dir=case_dir,
            config_path=args.config,
            struc_geometry_path=args.struc_geometry,
            aero_geometry_path=args.aero_geometry,
            solver_mode=args.solver_mode,
            requested_va_ms=args.va,
            target_depower_tape_length_m=target_depower_tape_length_m,
            steering_tape_final_extension_m=0.0,
            max_iter=args.max_iter,
        )
        row = summary_from_case(
            case_id,
            case_dir,
            meta,
            tracking,
            extra={
                "mode": args.solver_mode,
                "requested_V_a_ms": args.va,
                "requested_u_dp_ch9": args.udp,
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
                "target_depower_tape_u_dp": depower_tape_length_m_to_udp(
                    target_depower_tape_length_m
                ),
            },
        )
        write_json(case_dir / "case_summary.json", row)
        print(f"{case_id}: converged={row['converged']}")
    except Exception as exc:
        row = _failure_row(case_id, case_dir, args, exc)
        write_json(case_dir / "case_summary.json", row)
        write_csv(processed_root / "shape_validation_summary.csv", [row])
        write_markdown_table(processed_root / "shape_validation_summary.md", [row])
        raise

    summary_csv = processed_root / "shape_validation_summary.csv"
    write_csv(summary_csv, [row])
    write_markdown_table(processed_root / "shape_validation_summary.md", [row])
    write_json(
        processed_root / "shape_validation_run_manifest.json",
        {
            "case_id": case_id,
            "n_cases": 1,
            "output_root": str(args.output_root),
            "processed_data_dir": str(processed_root),
            "case_dir": str(case_dir),
            "summary_csv": str(summary_csv),
            "solver_mode": args.solver_mode,
            "requested_V_a_ms": float(args.va),
            "requested_u_dp_ch9": float(args.udp),
            "target_depower_tape_length_m": target_depower_tape_length_m,
            "config": str(args.config),
            "struc_geometry": str(args.struc_geometry),
            "aero_geometry": str(args.aero_geometry),
        },
    )


if __name__ == "__main__":
    main()
