from __future__ import annotations

import argparse
import logging
import os
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# Keep matplotlib, imported indirectly by aero2struc_level_1, from trying to
# write its cache into the user's home directory when this script runs headless.
os.environ.setdefault("MPLCONFIGDIR", "/tmp/askite_matplotlib")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

# Allow running this script directly from the repository root without requiring
# an editable install or PYTHONPATH setup.
PROJECT_DIR = Path(__file__).resolve().parents[2]
SRC_DIR = PROJECT_DIR / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from kitesim import aero2struc_level_1, aero2struc_level_1_old_pre_20_05_2026
from kitesim.utils import load_yaml


FORCE_REL_TOL = 1e-10
MOMENT_NORM_REL_TOL = 4e-3
MOMENT_COMPONENT_REL_TOL = 5e-3


@dataclass(frozen=True)
class Panel:
    aerodynamic_center: np.ndarray


@dataclass(frozen=True)
class CaseResult:
    implementation: str
    case_name: str
    force_rel_norm: float
    max_force_component_rel: float
    moment_rel_norm: float
    max_moment_component_rel: float
    passed: bool


def _resolve_level_1_geometry_path(project_dir: Path) -> Path:
    candidate_paths = [
        project_dir
        / "data"
        / "TUDELFT_V3_KITE"
        / "struc_geometry_level_1_manual_JULIA.yaml",
        project_dir / "data" / "TUDELFT_V3_KITE" / "struc_geometry_level_1_manual.yaml",
    ]
    for path in candidate_paths:
        if path.exists():
            return path
    raise FileNotFoundError(
        "Could not find a level-1 geometry YAML. Checked:\n"
        + "\n".join(str(p) for p in candidate_paths)
    )


def _load_wing_nodes(struc_geometry: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build an indexed wing-node array without parsing the bridle subsystem."""
    wing_rows = struc_geometry["wing_particles"]["data"]
    max_node_idx = max(int(row[0]) for row in wing_rows)
    nodes = np.zeros((max_node_idx + 1, 3), dtype=float)
    nodes[0] = np.asarray(struc_geometry.get("bridle_point_node", [0.0, 0.0, 0.0]))

    le_indices = []
    te_indices = []
    for node_idx, x, y, z in wing_rows:
        node_idx = int(node_idx)
        nodes[node_idx] = [float(x), float(y), float(z)]
        if node_idx % 2:
            le_indices.append(node_idx)
        else:
            te_indices.append(node_idx)

    return nodes, np.asarray(le_indices, dtype=int), np.asarray(te_indices, dtype=int)


def _rotate_about_global_y(point: np.ndarray, pivot: np.ndarray, angle_deg: float):
    """Rotate so positive angle raises a positive-x trailing edge in z."""
    theta = np.deg2rad(angle_deg)
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)
    rel = point - pivot
    rotated = np.array(
        [
            cos_t * rel[0] - sin_t * rel[2],
            rel[1],
            sin_t * rel[0] + cos_t * rel[2],
        ],
        dtype=float,
    )
    return pivot + rotated


def _apply_asymmetric_airfoil_twist(
    nodes: np.ndarray,
    le_indices: np.ndarray,
    te_indices: np.ndarray,
    left_deg: float,
    right_deg: float,
) -> np.ndarray:
    """
    Rotate each LE/TE rib pair about its mid-chord.

    Positive-y ribs are treated as the left half of the wing. Positive angles
    rotate the local trailing edge upward.
    """
    twisted = nodes.copy()
    for le_idx, te_idx in zip(le_indices, te_indices):
        pivot = 0.5 * (twisted[le_idx] + twisted[te_idx])
        angle_deg = left_deg if pivot[1] >= 0.0 else right_deg
        twisted[le_idx] = _rotate_about_global_y(twisted[le_idx], pivot, angle_deg)
        twisted[te_idx] = _rotate_about_global_y(twisted[te_idx], pivot, angle_deg)
    return twisted


def _build_panel_centers(
    nodes: np.ndarray,
    le_indices: np.ndarray,
    te_indices: np.ndarray,
    panels_per_struct_section: int,
    chord_fraction: float = 0.25,
) -> np.ndarray:
    """Create aerodynamic centers from subdivided structural LE/TE sections."""
    le_order = np.argsort(nodes[le_indices, 1])
    le_sorted = le_indices[le_order]
    te_sorted = te_indices[le_order]

    panel_centers = []
    for section_idx in range(len(le_sorted) - 1):
        le_lo = nodes[le_sorted[section_idx]]
        le_hi = nodes[le_sorted[section_idx + 1]]
        te_lo = nodes[te_sorted[section_idx]]
        te_hi = nodes[te_sorted[section_idx + 1]]

        for panel_idx in range(panels_per_struct_section):
            eta = (panel_idx + 0.5) / panels_per_struct_section
            le_mid = (1.0 - eta) * le_lo + eta * le_hi
            te_mid = (1.0 - eta) * te_lo + eta * te_hi
            panel_centers.append(le_mid + chord_fraction * (te_mid - le_mid))

    return np.asarray(panel_centers, dtype=float)


def _synthetic_panel_forces(panel_centers: np.ndarray) -> np.ndarray:
    """
    Deterministic VSM-like forces with nonzero total moments in all components.

    The exact magnitudes are not important for the mapping test. They are chosen
    to keep all moment-component denominators well away from zero.
    """
    y = panel_centers[:, 1]
    y_norm = y / np.max(np.abs(y))
    drag = 55.0 * (1.0 + 0.25 * y_norm)
    side = 200.0 * (0.6 + 0.3 * y_norm)
    lift = 750.0 * (1.0 - 0.1 * np.abs(y_norm)) * (1.0 + 0.06 * y_norm)
    return np.column_stack([drag, side, lift])


def _mapping_metrics(
    panel_centers: np.ndarray,
    panel_forces: np.ndarray,
    structural_nodes: np.ndarray,
    mapped_forces: np.ndarray,
    ref_point: np.ndarray,
) -> dict:
    force_aero = np.sum(panel_forces, axis=0)
    force_struct = np.sum(mapped_forces, axis=0)
    force_error = force_struct - force_aero

    moment_aero = np.sum(np.cross(panel_centers - ref_point, panel_forces), axis=0)
    moment_struct = np.sum(
        np.cross(structural_nodes - ref_point, mapped_forces),
        axis=0,
    )
    moment_error = moment_struct - moment_aero

    return {
        "force_rel_norm": _relative_norm(force_error, force_aero),
        "max_force_component_rel": float(
            np.max(_relative_components(force_error, force_aero))
        ),
        "moment_rel_norm": _relative_norm(moment_error, moment_aero),
        "max_moment_component_rel": float(
            np.max(_relative_components(moment_error, moment_aero))
        ),
    }


def _relative_norm(error: np.ndarray, reference: np.ndarray) -> float:
    reference_norm = np.linalg.norm(reference)
    if reference_norm < 1e-12:
        return 0.0
    return float(np.linalg.norm(error) / reference_norm)


def _relative_components(error: np.ndarray, reference: np.ndarray) -> np.ndarray:
    return np.abs(error) / np.maximum(np.abs(reference), 1e-12)


def _run_mapping(
    implementation_name: str,
    module,
    case_name: str,
    nodes: np.ndarray,
    le_indices: np.ndarray,
    te_indices: np.ndarray,
    config_aero2struc: dict,
    panels_per_struct_section: int,
    ref_point: np.ndarray,
) -> CaseResult:
    panel_centers = _build_panel_centers(
        nodes,
        le_indices,
        te_indices,
        panels_per_struct_section=panels_per_struct_section,
    )
    panels = [Panel(center) for center in panel_centers]
    panel_forces = _synthetic_panel_forces(panel_centers)
    mapping = module.initialize_mapping(panels, nodes, le_indices, te_indices)

    # Intentionally call the current implementation without moment_aero_panel so
    # it exercises the same force-only bilinear path as the archived script.
    mapped_forces = module.main(
        config_aero2struc["coupling_method"],
        panel_forces,
        nodes,
        panel_centers,
        mapping,
        False,
        config_aero2struc,
    )
    metrics = _mapping_metrics(
        panel_centers=panel_centers,
        panel_forces=panel_forces,
        structural_nodes=nodes,
        mapped_forces=mapped_forces,
        ref_point=ref_point,
    )
    passed = (
        metrics["force_rel_norm"] <= FORCE_REL_TOL
        and metrics["max_force_component_rel"] <= FORCE_REL_TOL
        and metrics["moment_rel_norm"] <= MOMENT_NORM_REL_TOL
        and metrics["max_moment_component_rel"] <= MOMENT_COMPONENT_REL_TOL
    )
    return CaseResult(
        implementation=implementation_name,
        case_name=case_name,
        force_rel_norm=metrics["force_rel_norm"],
        max_force_component_rel=metrics["max_force_component_rel"],
        moment_rel_norm=metrics["moment_rel_norm"],
        max_moment_component_rel=metrics["max_moment_component_rel"],
        passed=passed,
    )


def _format_percent(value: float) -> str:
    return f"{100.0 * value:.6f}"


def _print_results(results: list[CaseResult], panels_per_struct_section: int) -> None:
    print("Aero-to-structure force/moment mapping verification")
    print(f"  panels per structural section: {panels_per_struct_section}")
    print(f"  force relative tolerance: {FORCE_REL_TOL:.1e}")
    print(f"  moment norm tolerance [%]: {100.0 * MOMENT_NORM_REL_TOL:.3f}")
    print(f"  moment component tolerance [%]: {100.0 * MOMENT_COMPONENT_REL_TOL:.3f}")
    print("")
    header = (
        "implementation",
        "case",
        "force norm [%]",
        "force comp max [%]",
        "moment norm [%]",
        "moment comp max [%]",
        "pass",
    )
    widths = (40, 32, 16, 20, 17, 21, 6)
    print(
        " ".join(
            text.ljust(width) for text, width in zip(header, widths)
        )
    )
    print("-" * (sum(widths) + len(widths) - 1))
    for result in results:
        row = (
            result.implementation,
            result.case_name,
            _format_percent(result.force_rel_norm),
            _format_percent(result.max_force_component_rel),
            _format_percent(result.moment_rel_norm),
            _format_percent(result.max_moment_component_rel),
            "yes" if result.passed else "no",
        )
        print(
            " ".join(
                text.ljust(width) for text, width in zip(row, widths)
            )
        )


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Verify force conservation and moment-error bounds for the "
            "level-1 aero-to-structure mapping scripts."
        )
    )
    parser.add_argument(
        "--panels-per-struct-section",
        type=int,
        default=None,
        help=(
            "Aerodynamic panels per structural section. Defaults to the "
            "TUDELFT_V3_KITE config value."
        ),
    )
    args = parser.parse_args()

    logging.getLogger().setLevel(logging.ERROR)

    config = load_yaml(PROJECT_DIR / "data" / "TUDELFT_V3_KITE" / "config.yaml")
    config_aero2struc = config["aero2struc"]
    panels_per_struct_section = (
        args.panels_per_struct_section
        if args.panels_per_struct_section is not None
        else int(config["aerodynamic"]["n_aero_panels_per_struc_section"])
    )
    if panels_per_struct_section <= 0:
        raise ValueError("--panels-per-struct-section must be positive.")

    geometry_path = _resolve_level_1_geometry_path(PROJECT_DIR)
    struc_geometry = load_yaml(geometry_path)
    nodes, le_indices, te_indices = _load_wing_nodes(struc_geometry)
    ref_point = np.asarray(config["aerodynamic"]["reference_point"], dtype=float)

    cases = [
        ("undeformed", nodes),
        (
            "left_up_5deg_right_down_10deg",
            _apply_asymmetric_airfoil_twist(
                nodes,
                le_indices,
                te_indices,
                left_deg=5.0,
                right_deg=-10.0,
            ),
        ),
        (
            "left_down_5deg_right_up_10deg",
            _apply_asymmetric_airfoil_twist(
                nodes,
                le_indices,
                te_indices,
                left_deg=-5.0,
                right_deg=10.0,
            ),
        ),
    ]
    implementations = [
        ("aero2struc_level_1", aero2struc_level_1),
        (
            "aero2struc_level_1_old_pre_20_05_2026",
            aero2struc_level_1_old_pre_20_05_2026,
        ),
    ]

    results = []
    for implementation_name, module in implementations:
        for case_name, case_nodes in cases:
            results.append(
                _run_mapping(
                    implementation_name=implementation_name,
                    module=module,
                    case_name=case_name,
                    nodes=case_nodes,
                    le_indices=le_indices,
                    te_indices=te_indices,
                    config_aero2struc=config_aero2struc,
                    panels_per_struct_section=panels_per_struct_section,
                    ref_point=ref_point,
                )
            )

    _print_results(results, panels_per_struct_section)
    return 0 if all(result.passed for result in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
