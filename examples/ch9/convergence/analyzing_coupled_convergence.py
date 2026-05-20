import argparse
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

from kitesim.utils import load_sim_output, load_yaml

CH9_DIR = Path(__file__).resolve().parents[1]
if str(CH9_DIR) not in sys.path:
    sys.path.insert(0, str(CH9_DIR))

from ch9_analysis_utils import (
    DEFAULT_AERO_GEOMETRY,
    DEFAULT_CONFIG,
    DEFAULT_STRUC_GEOMETRY,
    PROJECT_DIR,
    parse_formats,
    run_coupled_case,
    save_figure,
    udp_to_depower_tape_length_m,
    write_csv,
    write_json,
    write_markdown_table,
    write_yaml,
)

MISSING_COEFF_MESSAGE = (
    "This result folder does not contain per-iteration aerodynamic coefficients.\n"
    "Re-run the coupled solver after adding the Section 9.3.1 tracking fields."
)


def build_parser():
    parser = argparse.ArgumentParser(
        description="Analyze Section 9.3.1 coupled VWT convergence from sim_output.h5."
    )
    parser.add_argument("--case-dir", required=True, type=Path)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_DIR / "results" / "ch9" / "convergence",
    )
    parser.add_argument("--case-label", default="ASKITE coupled VWT")
    parser.add_argument("--s-ref", type=float, default=None)
    parser.add_argument("--skip-initial-zero", action="store_true")
    parser.add_argument("--format", default="pdf")
    parser.add_argument("--stagnation-window", type=int, default=5)
    parser.add_argument(
        "--rerun",
        action="store_true",
        help="Run the coupled solver again and analyze that new result.",
    )
    parser.add_argument(
        "--rerun-case-dir",
        type=Path,
        default=None,
        help="Directory for the new rerun case. If omitted, a suffix is added to --case-dir.",
    )
    parser.add_argument(
        "--tol",
        type=float,
        default=1.0,
        help="Residual convergence tolerance used for rerun (aero_structural_solver.tol).",
    )
    parser.add_argument(
        "--rerun-max-iter",
        type=int,
        default=None,
        help="Optional max iterations for rerun. If omitted, uses config default.",
    )
    parser.add_argument(
        "--rerun-udp",
        type=float,
        default=None,
        help=(
            "Optional UDP value for rerun. This rewrites depower_tape l0 directly in "
            "a copied structural YAML and uses that file for rerun."
        ),
    )
    return parser


def _load_case_inputs(case_dir):
    path = Path(case_dir) / "case_inputs.json"
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def _existing_or_default(case_dir, name, default_path):
    candidate = Path(case_dir) / name
    return candidate if candidate.exists() else Path(default_path)


def _path_from_case_inputs(case_inputs, key, fallback):
    raw = case_inputs.get(key, None)
    if raw is None:
        return Path(fallback)
    candidate = Path(raw)
    return candidate if candidate.exists() else Path(fallback)


def _build_rerun_case_dir(case_dir, tol, rerun_udp=None):
    tol_tag = str(float(tol)).replace(".", "p")
    if rerun_udp is None:
        return Path(f"{case_dir}_rerun_tol_{tol_tag}")
    udp_tag = f"{int(round(float(rerun_udp) * 1000)):04d}"
    return Path(f"{case_dir}_rerun_tol_{tol_tag}_udp_{udp_tag}")


def _rewrite_depower_tape_l0_for_udp(src_struc_geometry_path, dst_case_dir, udp):
    target_len_m = float(udp_to_depower_tape_length_m(float(udp)))
    data = load_yaml(Path(src_struc_geometry_path))

    if "bridle_elements" not in data:
        raise ValueError("Missing 'bridle_elements' in structural geometry YAML.")

    bridle_elements = data["bridle_elements"]
    headers = list(bridle_elements.get("headers", []))
    rows = bridle_elements.get("data", [])

    if "name" not in headers or "l0" not in headers:
        raise ValueError("Expected 'name' and 'l0' headers in bridle_elements.")

    idx_name = headers.index("name")
    idx_l0 = headers.index("l0")

    found = False
    for row in rows:
        if len(row) <= max(idx_name, idx_l0):
            continue
        if str(row[idx_name]) == "depower_tape":
            row[idx_l0] = target_len_m
            found = True
            break

    if not found:
        raise ValueError("Could not find depower_tape entry in bridle_elements.data.")

    dst_case_dir = Path(dst_case_dir)
    dst_case_dir.mkdir(parents=True, exist_ok=True)
    udp_tag = f"{int(round(float(udp) * 1000)):04d}"
    dst_path = dst_case_dir / f"struc_geometry_udp_{udp_tag}.yaml"
    write_yaml(dst_path, data)
    return dst_path


def relative_change(values, window):
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size <= 1:
        return np.nan
    window = min(max(1, int(window)), finite.size - 1)
    denom = max(abs(finite[-1]), 1e-12)
    return float(abs(finite[-1] - finite[-1 - window]) / denom)


def main():
    args = build_parser().parse_args()
    case_dir = args.case_dir

    if args.rerun:
        case_inputs = _load_case_inputs(case_dir)
        rerun_case_dir = (
            args.rerun_case_dir
            if args.rerun_case_dir is not None
            else _build_rerun_case_dir(case_dir, args.tol, args.rerun_udp)
        )

        config_path = _path_from_case_inputs(
            case_inputs,
            "config_path",
            _existing_or_default(case_dir, "config.yaml", DEFAULT_CONFIG),
        )
        struc_geometry_path = _path_from_case_inputs(
            case_inputs,
            "struc_geometry_path",
            _existing_or_default(
                case_dir, "struc_geometry.yaml", DEFAULT_STRUC_GEOMETRY
            ),
        )
        aero_geometry_path = _path_from_case_inputs(
            case_inputs,
            "aero_geometry_path",
            _existing_or_default(case_dir, "aero_geometry.yaml", DEFAULT_AERO_GEOMETRY),
        )

        target_depower_tape_length_m = case_inputs.get(
            "target_depower_tape_length_m", None
        )
        power_tape_final_extension_m = case_inputs.get(
            "power_tape_final_extension_m", None
        )

        if args.rerun_udp is not None:
            struc_geometry_path = _rewrite_depower_tape_l0_for_udp(
                src_struc_geometry_path=struc_geometry_path,
                dst_case_dir=rerun_case_dir,
                udp=args.rerun_udp,
            )
            # When UDP is controlled by direct YAML l0 rewrite, do not also
            # request actuation through target/final extension inputs.
            target_depower_tape_length_m = None
            power_tape_final_extension_m = None

        run_coupled_case(
            case_dir=rerun_case_dir,
            config_path=config_path,
            struc_geometry_path=struc_geometry_path,
            aero_geometry_path=aero_geometry_path,
            solver_mode=case_inputs.get("solver_mode", "level_1"),
            requested_va_ms=case_inputs.get("requested_va_ms", None),
            target_depower_tape_length_m=target_depower_tape_length_m,
            power_tape_final_extension_m=power_tape_final_extension_m,
            steering_tape_final_extension_m=float(
                case_inputs.get("steering_tape_final_extension_m", 0.0)
            ),
            config_overrides={"aero_structural_solver.tol": float(args.tol)},
            max_iter=args.rerun_max_iter,
        )
        case_dir = rerun_case_dir

    h5_path = case_dir / "sim_output.h5"
    meta, tracking = load_sim_output(h5_path)

    required = [
        "C_L_wing",
        "C_D_wing",
        "glide_ratio_wing",
        "residual_norm",
        "positions",
    ]
    if any(name not in tracking for name in required):
        raise SystemExit(MISSING_COEFF_MESSAGE)
    if not np.isfinite(np.asarray(tracking["C_L_wing"], dtype=float)).any():
        raise SystemExit(MISSING_COEFF_MESSAGE)

    n_iter = int(meta.get("n_iter", len(tracking["C_L_wing"])))
    n_iter = min(n_iter, len(tracking["C_L_wing"]))
    iteration = np.arange(n_iter)
    cl = np.asarray(tracking["C_L_wing"][:n_iter], dtype=float)
    cd = np.asarray(tracking["C_D_wing"][:n_iter], dtype=float)
    glide = np.asarray(tracking["glide_ratio_wing"][:n_iter], dtype=float)
    residual = np.asarray(tracking["residual_norm"][:n_iter], dtype=float)
    geom_rel = np.asarray(
        tracking.get("geometry_update_rel_norm", np.full(n_iter, np.nan))[:n_iter],
        dtype=float,
    )
    aero_update = np.asarray(
        tracking.get("aero_force_update_norm", np.full(n_iter, np.nan))[:n_iter],
        dtype=float,
    )

    if args.skip_initial_zero and n_iter > 1:
        start = 1
        iteration = iteration[start:]
        cl = cl[start:]
        cd = cd[start:]
        glide = glide[start:]
        residual = residual[start:]
        geom_rel = geom_rel[start:]
        aero_update = aero_update[start:]

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.2), sharex=True)
    panels = [
        (axes[0, 0], cl, "$C_L$"),
        (axes[0, 1], cd, "$C_D$"),
        (axes[1, 0], glide, "$C_L/C_D$"),
    ]
    for ax, values, label in panels:
        ax.plot(iteration, values, "o-", color="0.15", linewidth=1.0, markersize=3)
        finite = np.isfinite(values)
        if finite.any():
            ax.plot(iteration[finite][-1], values[finite][-1], "o", color="0.0")
        ax.set_ylabel(label)
        ax.grid(True, color="0.88", linewidth=0.6)

    ax = axes[1, 1]
    ax.semilogy(
        iteration, residual, "o-", label="Residual", linewidth=1.0, markersize=3
    )
    if np.isfinite(geom_rel).any():
        ax.semilogy(
            iteration,
            geom_rel,
            "s-",
            label="Geometry update",
            linewidth=1.0,
            markersize=3,
        )
    if np.isfinite(aero_update).any():
        ax.semilogy(
            iteration,
            aero_update,
            "^-",
            label="Aero-force update",
            linewidth=1.0,
            markersize=3,
        )
    ax.set_ylabel("Norm")
    ax.grid(True, color="0.88", linewidth=0.6)
    ax.legend(frameon=False, fontsize=8)

    for ax in axes[1, :]:
        ax.set_xlabel("Coupling iteration")
    fig.suptitle(args.case_label)
    fig.tight_layout()
    save_figure(
        fig,
        args.output_dir,
        "fig_9_3_1_4_coupled_vwt_convergence",
        parse_formats(args.format),
    )

    finite_idx = np.where(np.isfinite(cl) & np.isfinite(cd) & np.isfinite(glide))[0]
    final_local_idx = int(finite_idx[-1]) if finite_idx.size else len(cl) - 1
    summary = {
        "case_label": args.case_label,
        "case_dir": str(case_dir),
        "source_h5": str(h5_path),
        "rerun_enabled": bool(args.rerun),
        "rerun_tol": float(args.tol) if args.rerun else np.nan,
        "V_a": float(meta.get("final_V_a", meta.get("va", np.nan))),
        "depower_tape_length_m": float(meta.get("final_depower_tape_length_m", np.nan)),
        "n_iter": int(meta.get("n_iter", n_iter)),
        "converged": bool(meta.get("converged", False)),
        "final_residual_norm": float(residual[final_local_idx]),
        "final_C_L": float(cl[final_local_idx]),
        "final_C_D": float(cd[final_local_idx]),
        "final_L_over_D": float(glide[final_local_idx]),
        "relative_CL_change_last_N": relative_change(cl, args.stagnation_window),
        "relative_CD_change_last_N": relative_change(cd, args.stagnation_window),
        "relative_L_over_D_change_last_N": relative_change(
            glide, args.stagnation_window
        ),
        "S_ref_m2": float(
            args.s_ref if args.s_ref is not None else meta.get("S_ref_m2", np.nan)
        ),
        "S_ref_source": str(meta.get("S_ref_source", "")),
    }
    rows = [summary]
    write_csv(args.output_dir / "coupled_convergence_summary.csv", rows)
    write_markdown_table(args.output_dir / "coupled_convergence_summary.md", rows)
    write_json(args.output_dir / "coupled_convergence_summary.json", summary)


if __name__ == "__main__":
    main()
