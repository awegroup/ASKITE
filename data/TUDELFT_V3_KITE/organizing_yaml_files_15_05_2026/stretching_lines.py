"""
Stretch bridle lines while holding the wing/canopy nodes fixed.

This is a small example script for the TU Delft V3 full/reduced geometry files.
It uses local two-stage bridle-relaxation helpers for both PSM and FEM.

The script saves a 2x2 plot with initial and stretched states for both models.
"""

from __future__ import annotations

import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterable

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import numpy as np
import yaml

SCRIPT_DIR = Path(__file__).resolve().parent
ASKITE_SRC = Path("/home/jellepoland/ownCloud/phd/code/ASKITE/src")
KITE_FEM_SRC = Path("/home/jellepoland/ownCloud/phd/code/kite_fem/src")
PSS_SRC = Path("/home/jellepoland/ownCloud/phd/code/Particle_System_Simulator/src")

for path in (ASKITE_SRC, KITE_FEM_SRC, PSS_SRC):
    if path.exists() and str(path) not in sys.path:
        sys.path.insert(0, str(path))

from kitesim import read_struc_geometry_yaml_level_2
from kite_fem.FEMStructure import FEM_structure
from kite_fem.Functions import fix_nodes
from PSS.particleSystem import ParticleSystem
from PSS.particleSystem.SpringDamper import SpringDamperType

RESULTS_DIR = SCRIPT_DIR / "results"

PSM_GEOMETRY_PATH = (
    SCRIPT_DIR / "results" / "struc_geometry_PSM_reduced_stretched_mapped_to_FEM.yaml"
)
FEM_GEOMETRY_PATH = SCRIPT_DIR / "struc_geometry_FEM_full.yaml"
OUTPUT_PNG = RESULTS_DIR / "stretching_lines.png"
OUTPUT_PDF = RESULTS_DIR / "stretching_lines.pdf"
OUTPUT_SVG = RESULTS_DIR / "stretching_lines.svg"
OUTPUT_PDF_AZ_180 = RESULTS_DIR / "stretching_lines_az_-180.pdf"
OUTPUT_SVG_AZ_180 = RESULTS_DIR / "stretching_lines_az_-180.svg"
OUTPUT_PDF_AZ_90 = RESULTS_DIR / "stretching_lines_az_-90.pdf"
OUTPUT_SVG_AZ_90 = RESULTS_DIR / "stretching_lines_az_-90.svg"
OUTPUT_REPORT_MD = RESULTS_DIR / "stretching_lines_report.md"

# Show an interactive final-only PSM/FEM comparison after writing outputs.
SHOW_FINAL_3D_OUTPUT = True

# PSM integration / linear-solve settings.
PSS_DT = 0.0005
PSS_T_STEPS = 10000
PSS_ABS_TOL = 1e-50
PSS_REL_TOL = 1e-5
PSS_MAX_ITER = 500

# PSM kinetic-damping convergence settings.
CONVERGENCE_TOL = 1e-6
SETTLE_WINDOW = 10
SYMMETRY_TOL = 1e-5

# Bridle pre-relaxation settings. Negative z pulls the bridle origin down,
# and are shared by both PSM and FEM two-stage relaxation.
BRIDLE_RELAX_PULL_FORCE = -50
BRIDLE_RELAX_SETTLE_FORCE = -25  # -50
BRIDLE_RELAX_ENABLE_SETTLE_STAGE = True
RELAX_MAX_STEPS_PER_STAGE = None

# FEM nonlinear-solver settings during two-stage bridle relaxation.
FEM_RELAX_MAX_ITERATIONS = 2000
FEM_RELAX_TOLERANCE = 0.01

BuildPSystem = Callable[
    [list[np.ndarray], list[float], list[list[Any]], list[float], dict, list[int]], Any
]


@dataclass
class PSMRelaxationResult:
    psystem: Any
    coords_initial: np.ndarray
    coords_relaxed: np.ndarray
    fixed_indices: list[int]
    stage_histories: dict[str, np.ndarray]
    stage_converged: dict[str, bool]
    stage_iterations: dict[str, int | None]


def two_stage_force_schedule(
    pull_force: float,
    settle_force: float,
    *,
    include_settle_stage: bool = True,
) -> list[tuple[str, float]]:
    schedule = [("pull", float(pull_force))]
    if include_settle_stage:
        schedule.append(("settle", float(settle_force)))
    return schedule


def _coerce_solver_converged(solve_result: Any) -> bool | None:
    if isinstance(solve_result, tuple) and len(solve_result) > 0:
        return bool(solve_result[0])
    if isinstance(solve_result, (bool, np.bool_)):
        return bool(solve_result)
    return None


def _yes_no(value: bool | None) -> str:
    if value is None:
        return "unknown"
    return "yes" if value else "no"


def _iter_text(value: int | None) -> str:
    return "unknown" if value is None else str(int(value))


def load_yaml(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as stream:
        return yaml.safe_load(stream)


def node_ids_from_table(table: dict) -> list[int]:
    return [int(row[0]) for row in table["data"]]


def infer_wing_node_indices(struc_geometry: dict) -> list[int]:
    return node_ids_from_table(struc_geometry["wing_particles"])


def infer_fixed_node_indices(struc_geometry: dict) -> list[int]:
    fixed_indices = [int(idx) for idx in struc_geometry.get("fixed_point_indices", [])]
    if not fixed_indices:
        raise ValueError(
            "Geometry YAML must define at least one fixed_point_indices entry."
        )
    return fixed_indices


def parse_psm_geometry(struc_geometry: dict):
    """Parse a reduced structural geometry YAML into PSS inputs."""
    nodes = [np.array(struc_geometry["bridle_point_node"], dtype=float)]
    masses = [float(struc_geometry.get("kcu_mass", 0.0))]

    def pss_idx(yaml_node_id):
        return int(yaml_node_id)

    for row in struc_geometry["wing_particles"]["data"]:
        node_id = int(row[0])
        nodes.append(np.array([float(row[1]), float(row[2]), float(row[3])]))
        masses.append(0.0)
        assert pss_idx(node_id) == len(nodes) - 1, (
            f"Wing particle id {node_id} does not match PSS index {len(nodes) - 1}. "
            "Particle IDs must be sequential and must account for KCU node 0."
        )

    wing_elem_dict = {
        row[0]: dict(zip(struc_geometry["wing_elements"]["headers"][1:], row[1:]))
        for row in struc_geometry["wing_elements"]["data"]
    }

    connectivity = []
    l0_arr = []
    k_arr = []
    c_arr = []
    linktype_arr = []

    for conn_name, ci, cj in struc_geometry["wing_connections"]["data"]:
        ci, cj = pss_idx(ci), pss_idx(cj)
        e = wing_elem_dict[conn_name]
        l0 = float(e["l0"])
        k = float(e["k"])
        c = float(e["c"])
        m_elem = float(e["m"])
        masses[ci] += m_elem / 2.0
        masses[cj] += m_elem / 2.0
        linktype = str(e["linktype"]).lower()
        connectivity.append([ci, cj, k, c, SpringDamperType(linktype)])
        l0_arr.append(l0)
        k_arr.append(k)
        c_arr.append(c)
        linktype_arr.append(linktype)

    for row in struc_geometry["bridle_particles"]["data"]:
        node_id = int(row[0])
        nodes.append(np.array([float(row[1]), float(row[2]), float(row[3])]))
        masses.append(0.0)
        assert pss_idx(node_id) == len(nodes) - 1, (
            f"Bridle particle id {node_id} does not match PSS index "
            f"{len(nodes) - 1}."
        )

    bridle_elem_dict = {
        row[0]: dict(zip(struc_geometry["bridle_elements"]["headers"][1:], row[1:]))
        for row in struc_geometry["bridle_elements"]["data"]
    }

    material_props = struc_geometry
    pulley_mass = float(struc_geometry.get("pulley_mass", 0.1))
    pulley_dict = {}
    conn_idx = len(connectivity)

    for conn_data in struc_geometry["bridle_connections"]["data"]:
        conn_name = conn_data[0]
        ci = pss_idx(conn_data[1])
        cj = pss_idx(conn_data[2])
        e = bridle_elem_dict[conn_name]
        l0 = float(e["l0"])
        mat = material_props[str(e["material"])]
        d = float(e["d"])
        area = np.pi * (d / 2.0) ** 2
        k_line = float(mat["youngs_modulus"]) * area / l0
        c_line = float(mat["damping_per_stiffness"]) * k_line
        m_line = float(mat["density"]) * area * l0
        masses[ci] += m_line / 2.0
        masses[cj] += m_line / 2.0
        linktype = str(e["linktype"]).lower()

        if len(conn_data) == 4:
            ck = pss_idx(conn_data[3])
            masses[cj] += pulley_mass

            len_ci_cj = np.linalg.norm(nodes[ci] - nodes[cj])
            len_cj_ck = np.linalg.norm(nodes[cj] - nodes[ck])
            total_len = len_ci_cj + len_cj_ck
            l0_ci_cj = (len_ci_cj / total_len) * l0
            l0_cj_ck = (len_cj_ck / total_len) * l0

            connectivity.append([ci, cj, k_line, c_line, SpringDamperType("pulley")])
            l0_arr.append(l0_ci_cj)
            k_arr.append(k_line)
            c_arr.append(c_line)
            linktype_arr.append("pulley")
            pulley_dict[str(conn_idx)] = np.array(
                [cj, ck, l0_cj_ck, l0_ci_cj, ci], dtype=float
            )
            conn_idx += 1

            connectivity.append([cj, ck, k_line, c_line, SpringDamperType("pulley")])
            l0_arr.append(l0_cj_ck)
            k_arr.append(k_line)
            c_arr.append(c_line)
            linktype_arr.append("pulley")
            pulley_dict[str(conn_idx)] = np.array(
                [cj, ci, l0_ci_cj, l0_cj_ck, ci], dtype=float
            )
            conn_idx += 1
        else:
            connectivity.append([ci, cj, k_line, c_line, SpringDamperType(linktype)])
            l0_arr.append(l0)
            k_arr.append(k_line)
            c_arr.append(c_line)
            linktype_arr.append(linktype)
            conn_idx += 1

    fixed_indices = [pss_idx(idx) for idx in infer_fixed_node_indices(struc_geometry)]

    return (
        nodes,
        masses,
        connectivity,
        l0_arr,
        k_arr,
        c_arr,
        linktype_arr,
        pulley_dict,
        fixed_indices,
    )


def build_psystem(
    nodes,
    masses,
    connectivity,
    l0_arr,
    pulley_dict,
    fixed_indices,
    dt=PSS_DT,
    t_steps=PSS_T_STEPS,
    abs_tol=PSS_ABS_TOL,
    rel_tol=PSS_REL_TOL,
    max_iter=PSS_MAX_ITER,
):
    """Instantiate a ParticleSystem and override rest lengths with YAML values."""
    vel_ini = np.zeros((len(nodes), 3))
    initial_conditions = []
    for i, (pos, vel, m) in enumerate(zip(nodes, vel_ini, masses)):
        initial_conditions.append([pos, vel, float(m), i in fixed_indices])

    pss_params = {
        "dt": dt,
        "t_steps": t_steps,
        "abs_tol": abs_tol,
        "rel_tol": rel_tol,
        "max_iter": max_iter,
        "pulley_other_line_pair": {k: v[:3] for k, v in pulley_dict.items()},
    }

    psystem = ParticleSystem(connectivity, initial_conditions, pss_params)

    for idx, current_l0 in enumerate(psystem.extract_rest_length):
        if str(idx) in pulley_dict:
            _cj, _ck, _l0_cj_ck, l0_ci_cj, _ci = pulley_dict[str(idx)]
            target_l0 = l0_ci_cj
        else:
            target_l0 = l0_arr[idx]
        psystem.update_rest_length(idx, target_l0 - current_l0)

    return psystem


def build_symmetry_mapping(coords, tol=SYMMETRY_TOL):
    pos_idx = [i for i, p in enumerate(coords) if p[1] > tol]
    neg_idx = [i for i, p in enumerate(coords) if p[1] < -tol]
    center_idx = [i for i, p in enumerate(coords) if abs(p[1]) <= tol]

    pairs = []
    for pi in pos_idx:
        mirrored = np.array([coords[pi][0], -coords[pi][1], coords[pi][2]])
        for ni in neg_idx:
            if np.allclose(coords[ni], mirrored, atol=tol):
                pairs.append((pi, ni))
                break

    pairs = np.array(pairs, dtype=int) if pairs else np.empty((0, 2), dtype=int)
    print(
        f"[PSM] symmetry mapping: {len(pairs)} mirror pairs, "
        f"{len(center_idx)} centre-plane nodes"
    )
    return {"pairs": pairs, "center_indices": center_idx}


def force_symmetry_psm(psystem, sym_map) -> None:
    particles = psystem._ParticleSystem__particles

    for src, mir in sym_map["pairs"]:
        x_src = particles[src].x.copy()
        particles[mir].update_pos(np.array([x_src[0], -x_src[1], x_src[2]]))

        v_src = particles[src].v.copy()
        particles[mir].update_vel(np.array([v_src[0], -v_src[1], v_src[2]]))

    for ci in sym_map["center_indices"]:
        x_c = particles[ci].x.copy()
        x_c[1] = 0.0
        particles[ci].update_pos(x_c)

        v_c = particles[ci].v.copy()
        v_c[1] = 0.0
        particles[ci].update_vel(v_c)


def _as_index_list(indices: Iterable[int]) -> list[int]:
    return [int(idx) for idx in indices]


def set_particle_fixed_states(
    psystem: Any,
    fixed_indices: Iterable[int],
    *,
    line_constraint_indices: Iterable[int] = (),
    line_constraint: Iterable[float] = (0.0, 0.0, 1.0),
) -> None:
    fixed = set(_as_index_list(fixed_indices))
    line_fixed = set(_as_index_list(line_constraint_indices))
    line_constraint = [float(value) for value in line_constraint]

    for idx, particle in enumerate(psystem.particles):
        if idx in line_fixed:
            particle.set_fixed(True, line_constraint, "line")
        elif idx in fixed:
            particle.set_fixed(True)
        else:
            particle.set_fixed(False)


def kinetic_damping_solve(
    psystem: Any,
    f_ext: np.ndarray | None = None,
    *,
    tol: float = 1e-6,
    settle_window: int = 10,
    max_steps: int | None = None,
    symmetry_map: dict | None = None,
    force_symmetry: Callable[[Any, dict], None] | None = None,
    debug_callback: Callable[[np.ndarray, int, float], None] | None = None,
) -> tuple[bool, np.ndarray]:
    if f_ext is None:
        f_ext = np.zeros(psystem.n * 3)
    else:
        f_ext = np.asarray(f_ext, dtype=float)

    t_steps = int(psystem.params["t_steps"] if max_steps is None else max_steps)
    e_kin_history: list[float] = []
    converged = False

    for step_i in range(t_steps):
        psystem.kin_damp_sim(f_ext)

        if symmetry_map is not None and force_symmetry is not None:
            force_symmetry(psystem, symmetry_map)

        x, v = psystem.x_v_current
        e_kin = float(np.linalg.norm(v * v))
        e_kin_history.append(e_kin)

        if debug_callback is not None:
            coords_now = psystem.x_v_current_3D[0].copy()
            debug_callback(coords_now, step_i, e_kin)

        if step_i >= settle_window and np.mean(e_kin_history[-settle_window:]) <= tol:
            converged = True
            break

        if np.isnan(x).any() or np.isnan(e_kin):
            break

    return converged, np.asarray(e_kin_history)


def make_nodal_force_vector(
    n_nodes: int,
    node_indices: Iterable[int],
    force: Iterable[float],
) -> np.ndarray:
    f_ext = np.zeros(int(n_nodes) * 3)
    force = np.asarray(force, dtype=float)
    for node_idx in node_indices:
        start = int(node_idx) * 3
        f_ext[start : start + 3] += force
    return f_ext


def _uses_cleaned_yaml_indexing(nodes: list[np.ndarray], psystem: Any) -> bool:
    return psystem.n == len(nodes) - 1


def _to_psystem_indices(
    indices: Iterable[int], nodes: list[np.ndarray], psystem: Any
) -> list[int]:
    indices = _as_index_list(indices)
    if _uses_cleaned_yaml_indexing(nodes, psystem):
        return [idx - 1 for idx in indices]
    return indices


def _restore_yaml_coordinates(
    coords_initial: np.ndarray,
    coords_psystem: np.ndarray,
) -> np.ndarray:
    if len(coords_psystem) == len(coords_initial) - 1:
        coords_full = coords_initial.copy()
        coords_full[1:] = coords_psystem
        return coords_full
    return coords_psystem.copy()


def relax_bridles_psm(
    *,
    nodes: list[np.ndarray],
    masses: list[float],
    connectivity: list[list[Any]],
    l0_arr: list[float],
    pulley_dict: dict,
    fixed_indices: Iterable[int],
    wing_node_indices: Iterable[int],
    bridle_origin_indices: Iterable[int],
    build_psystem: BuildPSystem,
    force_symmetry: Callable[[Any, dict], None] | None = None,
    build_symmetry_mapping: Callable[[np.ndarray], dict] | None = None,
    pull_force: float = 100.0,
    settle_force: float = 1.0,
    tol: float = 1e-6,
    settle_window: int = 10,
    max_steps_per_stage: int | None = None,
    keep_origin_position: bool = True,
    apply_origin_z_offset: bool = False,
    include_settle_stage: bool = True,
    debug_callback: Callable[[str, np.ndarray, int, float], None] | None = None,
) -> PSMRelaxationResult:
    fixed_indices = _as_index_list(fixed_indices)
    wing_node_indices = _as_index_list(wing_node_indices)
    bridle_origin_indices = _as_index_list(bridle_origin_indices)

    coords_initial = np.asarray(nodes, dtype=float)

    relaxation_fixed = sorted(set(fixed_indices).union(wing_node_indices))
    for idx in bridle_origin_indices:
        if idx in relaxation_fixed:
            relaxation_fixed.remove(idx)

    psystem = build_psystem(
        [np.asarray(node, dtype=float).copy() for node in nodes],
        list(masses),
        connectivity,
        l0_arr,
        pulley_dict,
        relaxation_fixed,
    )
    relaxation_fixed_ps = _to_psystem_indices(relaxation_fixed, nodes, psystem)
    bridle_origin_indices_ps = _to_psystem_indices(
        bridle_origin_indices, nodes, psystem
    )
    set_particle_fixed_states(
        psystem,
        relaxation_fixed_ps,
        line_constraint_indices=bridle_origin_indices_ps,
        line_constraint=(0.0, 0.0, 1.0),
    )

    symmetry_map = None
    if build_symmetry_mapping is not None:
        symmetry_map = build_symmetry_mapping(psystem.x_v_current_3D[0].copy())

    stage_histories: dict[str, np.ndarray] = {}
    stage_converged: dict[str, bool] = {}
    stage_iterations: dict[str, int | None] = {}

    for label, force_z in two_stage_force_schedule(
        pull_force,
        settle_force,
        include_settle_stage=include_settle_stage,
    ):
        f_ext = make_nodal_force_vector(
            psystem.n, bridle_origin_indices_ps, (0.0, 0.0, force_z)
        )
        stage_debug_callback = None
        if debug_callback is not None:
            stage_debug_callback = (
                lambda coords, step_i, e_kin, stage=label: debug_callback(
                    stage,
                    _restore_yaml_coordinates(coords_initial, coords),
                    step_i,
                    e_kin,
                )
            )
        converged, history = kinetic_damping_solve(
            psystem,
            f_ext,
            tol=tol,
            settle_window=settle_window,
            max_steps=max_steps_per_stage,
            symmetry_map=symmetry_map,
            force_symmetry=force_symmetry,
            debug_callback=stage_debug_callback,
        )
        stage_histories[label] = history
        stage_converged[label] = converged
        stage_iterations[label] = (len(history) - 1) if len(history) > 0 else None

    coords_relaxed = _restore_yaml_coordinates(
        coords_initial, psystem.x_v_current_3D[0].copy()
    )
    if bridle_origin_indices:
        initial_origin = coords_initial[bridle_origin_indices].mean(axis=0)
        relaxed_origin = coords_relaxed[bridle_origin_indices].mean(axis=0)
        z_offset = relaxed_origin[2] - initial_origin[2]
        if keep_origin_position:
            coords_relaxed += initial_origin - relaxed_origin
        if apply_origin_z_offset:
            # Match the FEM post-processing convention.
            coords_relaxed[:, 2] += z_offset

    final_psystem = build_psystem(
        [coord.copy() for coord in coords_relaxed],
        list(masses),
        connectivity,
        l0_arr,
        pulley_dict,
        fixed_indices,
    )

    return PSMRelaxationResult(
        psystem=final_psystem,
        coords_initial=coords_initial,
        coords_relaxed=coords_relaxed,
        fixed_indices=fixed_indices,
        stage_histories=stage_histories,
        stage_converged=stage_converged,
        stage_iterations=stage_iterations,
    )


def relax_bridles_fem_two_stage(
    kite: FEM_structure,
    wing_node_indices: Iterable[int],
    bridle_origin_indices: Iterable[int],
    *,
    pull_force: float,
    settle_force: float,
    include_settle_stage: bool = True,
    max_iterations: int = FEM_RELAX_MAX_ITERATIONS,
    tolerance: float = FEM_RELAX_TOLERANCE,
) -> tuple[FEM_structure, dict[str, bool | None], dict[str, int | None]]:
    wing_nodes = _as_index_list(wing_node_indices)
    origin_nodes = _as_index_list(bridle_origin_indices)

    kite = fix_nodes(kite, wing_nodes)
    initial_conditions = kite.initial_conditions
    spring_matrix = kite.spring_matrix
    pulley_matrix = kite.pulley_matrix
    beam_matrix = kite.beam_matrix

    fe = np.zeros(kite.N)
    for node_id in origin_nodes:
        # Free z DOF on origin nodes so an external z-force can move them.
        kite.bc[6 * node_id + 2] = True
        kite.fixed[6 * node_id + 2] = False

    stage_converged: dict[str, bool | None] = {}
    stage_iterations: dict[str, int | None] = {}
    for label, force_z in two_stage_force_schedule(
        pull_force,
        settle_force,
        include_settle_stage=include_settle_stage,
    ):
        for node_id in origin_nodes:
            fe[6 * node_id + 2] = force_z
        solve_result = kite.solve(
            fe,
            max_iterations=max_iterations,
            tolerance=tolerance,
            print_info=False,
        )
        stage_converged[label] = _coerce_solver_converged(solve_result)
        if getattr(kite, "iteration_history", None):
            stage_iterations[label] = int(kite.iteration_history[-1])
        else:
            stage_iterations[label] = None

    initcoords = np.reshape(kite.coords_init, (-1, 3))
    newcoords = np.reshape(kite.coords_current, (-1, 3))

    z_offset = newcoords[origin_nodes[0]][2] - initcoords[origin_nodes[0]][2]

    initial_conditions_new = []
    for node_id, (pos, vel, mass, fixed) in enumerate(initial_conditions):
        posnew = newcoords[node_id].copy()
        posnew[2] += z_offset
        initial_conditions_new.append([posnew, vel, mass, fixed])

    return (
        FEM_structure(
            initial_conditions_new, spring_matrix, pulley_matrix, beam_matrix
        ),
        stage_converged,
        stage_iterations,
    )


def build_raw_fem_kite(struc_geometry: dict):
    """Create a FEM kite without running the automatic ASKITE bridle relaxation."""
    (
        struc_nodes,
        m_arr,
        _struc_node_le_indices,
        _struc_node_te_indices,
        _power_tape_index,
        _steering_tape_indices,
        _pulley_node_indices,
        canopy_sections,
        strut_sections,
        _simplified_bridle_points,
        kite_connectivity_arr,
        _bridle_connectivity_arr,
        _bridle_diameter_arr,
        l0_arr,
        k_arr,
        c_arr,
        linktype_arr,
        _pulley_line_indices,
        pulley_line_to_other_node_pair_dict,
    ) = read_struc_geometry_yaml_level_2.main(struc_geometry)

    initial_conditions = []
    vel_ini = np.zeros((len(struc_nodes), 3))
    fixed_set = {int(i) for i in struc_geometry.get("fixed_point_indices", [])}
    for node_id, pos in enumerate(struc_nodes):
        initial_conditions.append(
            [
                np.asarray(pos, dtype=float),
                vel_ini[node_id],
                float(m_arr[node_id]),
                node_id in fixed_set,
            ]
        )

    spring_matrix = []
    pulley_matrix = []
    beam_matrix = []
    seen_pulley_triplets = set()

    for idx, (cicj, k, c, l0, linktype) in enumerate(
        zip(kite_connectivity_arr, k_arr, c_arr, l0_arr, linktype_arr)
    ):
        ci, cj = int(cicj[0]), int(cicj[1])
        lt = str(linktype).lower()

        if lt == "pulley":
            map_val = pulley_line_to_other_node_pair_dict.get(str(idx))
            if map_val is None:
                spring_matrix.append([ci, cj, float(k), float(c), float(l0), lt])
                continue

            try:
                ci_map, cj_map, ck = int(map_val[3]), int(map_val[4]), int(map_val[5])
            except Exception:
                cj_map, ck = int(map_val[0]), int(map_val[1])
                ci_map = ci

            triplet = (ci_map, cj_map, ck)
            if triplet in seen_pulley_triplets:
                continue
            seen_pulley_triplets.add(triplet)

            if ci_map != cj_map:
                pulley_matrix.append(
                    [ci_map, cj_map, ck, float(k), float(c), float(l0)]
                )
        elif lt == "inflatable_beam":
            beam_matrix.append([ci, cj, float(k), float(c), float(l0)])
        else:
            spring_matrix.append([ci, cj, float(k), float(c), float(l0), lt])

    kite = FEM_structure(
        initial_conditions=initial_conditions,
        spring_matrix=spring_matrix,
        pulley_matrix=pulley_matrix,
        beam_matrix=beam_matrix,
    )
    wing_nodes = sorted(
        {int(node) for section in canopy_sections + strut_sections for node in section}
    )
    return kite, wing_nodes


def origin_location_change(
    coords_initial: np.ndarray,
    coords_final: np.ndarray,
    origin_nodes: Iterable[int],
) -> tuple[np.ndarray, np.ndarray]:
    origin_indices = _as_index_list(origin_nodes)
    initial_origin = np.asarray(coords_initial[origin_indices], dtype=float).mean(
        axis=0
    )
    final_origin = np.asarray(coords_final[origin_indices], dtype=float).mean(axis=0)
    return initial_origin, final_origin


def build_bridle_line_catalog(struc_geometry: dict) -> list[dict[str, Any]]:
    if "bridle_elements" in struc_geometry:
        line_table = struc_geometry["bridle_elements"]
        rest_length_key = "l0"
    else:
        line_table = struc_geometry["bridle_lines"]
        rest_length_key = "rest_length"

    line_props = {
        row[0]: dict(zip(line_table["headers"][1:], row[1:]))
        for row in line_table["data"]
    }

    name_counts: dict[str, int] = {}
    line_catalog = []
    for row in struc_geometry["bridle_connections"]["data"]:
        name = str(row[0])
        name_counts[name] = name_counts.get(name, 0) + 1
        node_ids = [int(value) for value in row[1:]]
        props = line_props[name]
        line_catalog.append(
            {
                "name": name,
                "occurrence": name_counts[name],
                "node_ids": node_ids,
                "rest_length": float(props[rest_length_key]),
                "is_pulley": len(node_ids) == 3,
            }
        )

    return line_catalog


def bridle_line_lengths(
    coords: np.ndarray,
    node_ids: Iterable[int],
) -> tuple[float, list[float]]:
    node_ids = _as_index_list(node_ids)
    if len(node_ids) == 2:
        segment_length = float(
            np.linalg.norm(coords[node_ids[0]] - coords[node_ids[1]])
        )
        return segment_length, [segment_length]
    if len(node_ids) == 3:
        branch_1 = float(np.linalg.norm(coords[node_ids[0]] - coords[node_ids[1]]))
        branch_2 = float(np.linalg.norm(coords[node_ids[1]] - coords[node_ids[2]]))
        return branch_1 + branch_2, [branch_1, branch_2]
    raise ValueError(f"Unsupported bridle connection with nodes={node_ids}")


def build_bridle_line_stretch_rows(
    struc_geometry: dict,
    coords_initial: np.ndarray,
    coords_final: np.ndarray,
) -> list[dict[str, Any]]:
    rows = []
    for line_data in build_bridle_line_catalog(struc_geometry):
        name = line_data["name"]
        occurrence = line_data["occurrence"]
        node_ids = line_data["node_ids"]
        rest_length = line_data["rest_length"]

        initial_length, initial_segments = bridle_line_lengths(coords_initial, node_ids)
        final_length, final_segments = bridle_line_lengths(coords_final, node_ids)

        node_label = "-".join(str(node_id) for node_id in node_ids)
        rows.append(
            {
                "name": name,
                "occurrence": occurrence,
                "type": "pulley" if line_data["is_pulley"] else "line",
                "nodes": node_label,
                "rest_length": rest_length,
                "initial_length": initial_length,
                "final_length": final_length,
                "stretch_initial": initial_length - rest_length,
                "stretch_final": final_length - rest_length,
                "stretch_initial_pct": (
                    100.0 * (initial_length - rest_length) / rest_length
                    if rest_length > 0.0
                    else np.nan
                ),
                "stretch_final_pct": (
                    100.0 * (final_length - rest_length) / rest_length
                    if rest_length > 0.0
                    else np.nan
                ),
                "delta_length": final_length - initial_length,
                "initial_segments": initial_segments,
                "final_segments": final_segments,
            }
        )
    return rows


def _format_vec(vec: np.ndarray) -> str:
    return f"[{vec[0]:.6f}, {vec[1]:.6f}, {vec[2]:.6f}]"


def _format_segments(segments: list[float]) -> str:
    if len(segments) == 1:
        return f"{segments[0]:.4f}"
    return f"({segments[0]:.4f}, {segments[1]:.4f})"


def summarize_final_stretch(rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    min_row = min(rows, key=lambda row: row["stretch_final"])
    max_row = max(rows, key=lambda row: row["stretch_final"])
    return {"min": min_row, "max": max_row}


def summarize_final_stretch_pct(
    rows: list[dict[str, Any]],
) -> dict[str, dict[str, Any]]:
    finite_rows = [row for row in rows if np.isfinite(row["stretch_final_pct"])]
    min_row = min(finite_rows, key=lambda row: row["stretch_final_pct"])
    max_row = max(finite_rows, key=lambda row: row["stretch_final_pct"])
    return {"min": min_row, "max": max_row}


def _line_ref(row: dict[str, Any]) -> str:
    return f"{row['name']}#{row['occurrence']} ({row['nodes']})"


def build_markdown_table(rows: list[dict[str, Any]]) -> str:
    lines = [
        "| Line | Type | Nodes | l0 [m] | Init [m] | Final [m] | Stretch init [m] | Stretch final [m] | Stretch final [%] | Delta [m] | Init branches [m] | Final branches [m] |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |",
    ]
    for row in rows:
        lines.append(
            "| "
            f"{row['name']}#{row['occurrence']} | "
            f"{row['type']} | "
            f"{row['nodes']} | "
            f"{row['rest_length']:.4f} | "
            f"{row['initial_length']:.4f} | "
            f"{row['final_length']:.4f} | "
            f"{row['stretch_initial']:+.4f} | "
            f"{row['stretch_final']:+.4f} | "
            f"{row['stretch_final_pct']:+.2f} | "
            f"{row['delta_length']:+.4f} | "
            f"{_format_segments(row['initial_segments'])} | "
            f"{_format_segments(row['final_segments'])} |"
        )
    return "\n".join(lines)


def compare_models_by_line_percent(
    psm_rows: list[dict[str, Any]],
    fem_rows: list[dict[str, Any]],
) -> dict[str, Any]:
    psm_map = {(row["name"], row["occurrence"]): row for row in psm_rows}
    fem_map = {(row["name"], row["occurrence"]): row for row in fem_rows}
    shared_keys = sorted(set(psm_map).intersection(fem_map))

    comparisons = []
    for key in shared_keys:
        psm_row = psm_map[key]
        fem_row = fem_map[key]
        if not (
            np.isfinite(psm_row["stretch_final_pct"])
            and np.isfinite(fem_row["stretch_final_pct"])
        ):
            continue
        diff_pct = psm_row["stretch_final_pct"] - fem_row["stretch_final_pct"]
        comparisons.append(
            {
                "key": key,
                "line_label": f"{key[0]}#{key[1]}",
                "psm_pct": psm_row["stretch_final_pct"],
                "fem_pct": fem_row["stretch_final_pct"],
                "diff_pct": diff_pct,
                "abs_diff_pct": abs(diff_pct),
            }
        )

    if not comparisons:
        return {"count": 0, "comparisons": []}

    max_abs = max(comparisons, key=lambda row: row["abs_diff_pct"])
    mean_abs = float(np.mean([row["abs_diff_pct"] for row in comparisons]))
    return {
        "count": len(comparisons),
        "comparisons": comparisons,
        "max_abs": max_abs,
        "mean_abs": mean_abs,
    }


def analyze_cross_model_wing_node_proximity(
    coords_psm: np.ndarray,
    wing_nodes_psm: Iterable[int],
    coords_fem: np.ndarray,
    wing_nodes_fem: Iterable[int],
    *,
    y_abs_limit: float = 1.0,
    duplicate_tol: float = 1e-9,
    close_tol: float = 0.25,
) -> dict[str, Any]:
    selected_psm = [
        int(node)
        for node in wing_nodes_psm
        if 0 <= int(node) < len(coords_psm)
        and abs(float(coords_psm[int(node)][1])) < y_abs_limit
    ]
    selected_fem = [
        int(node)
        for node in wing_nodes_fem
        if 0 <= int(node) < len(coords_fem)
        and abs(float(coords_fem[int(node)][1])) < y_abs_limit
    ]
    selected_psm = sorted(set(selected_psm))
    selected_fem = sorted(set(selected_fem))

    duplicates = []
    close_pairs = []
    for ni in selected_psm:
        for nj in selected_fem:
            dist = float(np.linalg.norm(coords_psm[ni] - coords_fem[nj]))
            if dist <= duplicate_tol:
                duplicates.append((ni, nj, dist))
            elif dist <= close_tol:
                close_pairs.append((ni, nj, dist))

    close_pairs.sort(key=lambda item: item[2])
    return {
        "selected_nodes_psm": selected_psm,
        "selected_nodes_fem": selected_fem,
        "duplicates": duplicates,
        "close_pairs": close_pairs,
    }


def format_cross_model_wing_proximity_section(proximity: dict[str, Any]) -> list[str]:
    lines = ["### Cross-Model Wing Node Proximity (PSM vs FEM, |y| < 1 m)"]
    lines.append(f"Selected PSM wing nodes: {len(proximity['selected_nodes_psm'])}")
    lines.append(f"Selected FEM wing nodes: {len(proximity['selected_nodes_fem'])}")

    if proximity["duplicates"]:
        lines.append("Exact cross-model coordinate duplicates:")
        for ni, nj, _dist in proximity["duplicates"]:
            lines.append(f"- PSM {ni} and FEM {nj}")
    else:
        lines.append("Exact cross-model coordinate duplicates: none")

    if proximity["close_pairs"]:
        lines.append("Close cross-model pairs (distance <= 0.25 m):")
        for ni, nj, dist in proximity["close_pairs"][:10]:
            lines.append(f"- PSM {ni} and FEM {nj}: {dist:.4f} m")
    else:
        lines.append("Close cross-model pairs (distance <= 0.25 m): none")

    lines.append("")
    return lines


def write_markdown_report(
    psm_result: dict[str, Any], fem_result: dict[str, Any]
) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    psm_origin_initial, psm_origin_final = origin_location_change(
        psm_result["coords_initial"],
        psm_result["coords_final"],
        psm_result["origin_nodes"],
    )
    fem_origin_initial, fem_origin_final = origin_location_change(
        fem_result["coords_initial"],
        fem_result["coords_final"],
        fem_result["origin_nodes"],
    )
    origin_final_delta = fem_origin_final - psm_origin_final

    line_comparison = compare_models_by_line_percent(
        psm_result["line_rows"], fem_result["line_rows"]
    )

    cross_model_proximity = analyze_cross_model_wing_node_proximity(
        psm_result["coords_initial"],
        psm_result["wing_nodes"],
        fem_result["coords_initial"],
        fem_result["wing_nodes"],
    )

    sections = [
        "# Stretching Lines Report",
        "",
        "## Model Comparison Summary",
        "",
        "Solver convergence (pull/settle): "
        f"PSM {_yes_no(psm_result['stage_converged'].get('pull'))}/{_yes_no(psm_result['stage_converged'].get('settle'))}, "
        f"FEM {_yes_no(fem_result['stage_converged'].get('pull'))}/{_yes_no(fem_result['stage_converged'].get('settle'))}",
        "Solver stage iterations (pull/settle): "
        f"PSM {_iter_text(psm_result['stage_iterations'].get('pull'))}/{_iter_text(psm_result['stage_iterations'].get('settle'))}, "
        f"FEM {_iter_text(fem_result['stage_iterations'].get('pull'))}/{_iter_text(fem_result['stage_iterations'].get('settle'))}",
        "",
        f"Final origin PSM: {_format_vec(psm_origin_final)}  ",
        f"Final origin FEM: {_format_vec(fem_origin_final)}  ",
        f"Final origin delta (FEM - PSM): {_format_vec(origin_final_delta)}  ",
        f"Final origin delta norm: {np.linalg.norm(origin_final_delta):.6f} m",
        "",
    ]

    if line_comparison["count"] > 0:
        max_abs = line_comparison["max_abs"]
        sections.extend(
            [
                f"Shared line instances compared by final stretch %: {line_comparison['count']}  ",
                f"Mean absolute difference in final stretch %: {line_comparison['mean_abs']:.2f} %-points  ",
                "Largest absolute difference in final stretch %: "
                f"{max_abs['line_label']} = {max_abs['abs_diff_pct']:.2f} %-points "
                f"(PSM {max_abs['psm_pct']:+.2f}%, FEM {max_abs['fem_pct']:+.2f}%)",
                "",
            ]
        )
    else:
        sections.extend(
            [
                "Shared line instances compared by final stretch %: none",
                "",
            ]
        )

    sections.extend(format_cross_model_wing_proximity_section(cross_model_proximity))

    for label, result in (("PSM", psm_result), ("FEM", fem_result)):
        initial_origin, final_origin = origin_location_change(
            result["coords_initial"], result["coords_final"], result["origin_nodes"]
        )
        stretch_summary = summarize_final_stretch(result["line_rows"])
        stretch_summary_pct = summarize_final_stretch_pct(result["line_rows"])

        sections.extend(
            [
                f"## {label}",
                "",
                f"Origin initial: {_format_vec(initial_origin)}  ",
                f"Origin final: {_format_vec(final_origin)}  ",
                "Solver converged (pull/settle): "
                f"{_yes_no(result['stage_converged'].get('pull'))}/{_yes_no(result['stage_converged'].get('settle'))}  ",
                "Solver stage iterations (pull/settle): "
                f"{_iter_text(result['stage_iterations'].get('pull'))}/{_iter_text(result['stage_iterations'].get('settle'))}  ",
                f"Min final stretch: {_line_ref(stretch_summary['min'])} = {stretch_summary['min']['stretch_final']:+.4f} m  ",
                f"Max final stretch: {_line_ref(stretch_summary['max'])} = {stretch_summary['max']['stretch_final']:+.4f} m  ",
                f"Min final stretch %: {_line_ref(stretch_summary_pct['min'])} = {stretch_summary_pct['min']['stretch_final_pct']:+.2f}%  ",
                f"Max final stretch %: {_line_ref(stretch_summary_pct['max'])} = {stretch_summary_pct['max']['stretch_final_pct']:+.2f}%",
                "",
                build_markdown_table(result["line_rows"]),
                "",
            ]
        )

    OUTPUT_REPORT_MD.write_text("\n".join(sections), encoding="utf-8")
    print(f"[stretching_lines] wrote {OUTPUT_REPORT_MD}")


def run_psm_stretch():
    struc_geometry = load_yaml(PSM_GEOMETRY_PATH)
    (
        nodes,
        masses,
        connectivity,
        l0_arr,
        _k_arr,
        _c_arr,
        _linktype_arr,
        pulley_dict,
        fixed_indices,
    ) = parse_psm_geometry(struc_geometry)

    wing_nodes = infer_wing_node_indices(struc_geometry)
    result = relax_bridles_psm(
        nodes=nodes,
        masses=masses,
        connectivity=connectivity,
        l0_arr=l0_arr,
        pulley_dict=pulley_dict,
        fixed_indices=fixed_indices,
        wing_node_indices=wing_nodes,
        bridle_origin_indices=fixed_indices,
        build_psystem=build_psystem,
        force_symmetry=force_symmetry_psm,
        build_symmetry_mapping=build_symmetry_mapping,
        pull_force=BRIDLE_RELAX_PULL_FORCE,
        settle_force=BRIDLE_RELAX_SETTLE_FORCE,
        include_settle_stage=BRIDLE_RELAX_ENABLE_SETTLE_STAGE,
        tol=CONVERGENCE_TOL,
        settle_window=SETTLE_WINDOW,
        max_steps_per_stage=RELAX_MAX_STEPS_PER_STAGE,
        keep_origin_position=False,
        apply_origin_z_offset=True,
    )

    print(
        "[PSM] bridle stretch: "
        f"pull converged={result.stage_converged.get('pull')}, "
        f"settle converged={result.stage_converged.get('settle')}"
    )
    print(
        "[PSM] bridle stretch iterations: "
        f"pull={_iter_text(result.stage_iterations.get('pull'))}, "
        f"settle={_iter_text(result.stage_iterations.get('settle'))}"
    )
    return {
        "coords_initial": result.coords_initial,
        "coords_final": result.coords_relaxed,
        "connectivity": [(int(row[0]), int(row[1])) for row in connectivity],
        "wing_nodes": wing_nodes,
        "origin_nodes": fixed_indices,
        "stage_converged": result.stage_converged,
        "stage_iterations": result.stage_iterations,
        "line_rows": build_bridle_line_stretch_rows(
            struc_geometry, result.coords_initial, result.coords_relaxed
        ),
    }


def run_fem_stretch():
    struc_geometry = load_yaml(FEM_GEOMETRY_PATH)
    kite, wing_nodes = build_raw_fem_kite(struc_geometry)
    fixed_nodes = infer_fixed_node_indices(struc_geometry)
    origin_nodes = fixed_nodes

    kite = fix_nodes(kite, fixed_nodes)
    coords_initial = kite.coords_current.reshape((-1, 3)).copy()

    kite_stretched, stage_converged, stage_iterations = relax_bridles_fem_two_stage(
        kite,
        wing_nodes,
        origin_nodes,
        pull_force=BRIDLE_RELAX_PULL_FORCE,
        settle_force=BRIDLE_RELAX_SETTLE_FORCE,
        include_settle_stage=BRIDLE_RELAX_ENABLE_SETTLE_STAGE,
    )
    print(
        "[FEM] bridle stretch: "
        f"pull converged={_yes_no(stage_converged.get('pull'))}, "
        f"settle converged={_yes_no(stage_converged.get('settle'))}"
    )
    print(
        "[FEM] bridle stretch iterations: "
        f"pull={_iter_text(stage_iterations.get('pull'))}, "
        f"settle={_iter_text(stage_iterations.get('settle'))}"
    )
    kite_stretched = fix_nodes(kite_stretched, fixed_nodes)
    coords_final = kite_stretched.coords_current.reshape((-1, 3)).copy()

    return {
        "coords_initial": coords_initial,
        "coords_final": coords_final,
        "structure": kite_stretched,
        "wing_nodes": wing_nodes,
        "origin_nodes": origin_nodes,
        "stage_converged": stage_converged,
        "stage_iterations": stage_iterations,
        "line_rows": build_bridle_line_stretch_rows(
            struc_geometry, coords_initial, coords_final
        ),
    }


def _split_line_ending(line: str) -> tuple[str, str]:
    if line.endswith("\r\n"):
        return line[:-2], "\r\n"
    if line.endswith("\n"):
        return line[:-1], "\n"
    return line, ""


def _format_yaml_float(value: float) -> str:
    value = float(value)
    if abs(value) < 5e-13:
        value = 0.0
    text = f"{value:.8f}".rstrip("0").rstrip(".")
    if text in {"", "-0"}:
        return "0"
    return text


def _format_yaml_position(values: Iterable[float]) -> str:
    return ", ".join(_format_yaml_float(value) for value in values)


def coords_shifted_to_force_origin(
    coords: np.ndarray,
    force_node_indices: Iterable[int],
) -> np.ndarray:
    coords_shifted = np.asarray(coords, dtype=float).copy()
    force_nodes = _as_index_list(force_node_indices)
    if not force_nodes:
        raise ValueError("At least one force node is required to shift stretched YAML.")

    out_of_bounds = [
        idx for idx in force_nodes if idx < 0 or idx >= len(coords_shifted)
    ]
    if out_of_bounds:
        raise ValueError(
            f"Force node indices out of bounds for stretched YAML: {out_of_bounds}"
        )

    force_origin_z = float(coords_shifted[force_nodes, 2].mean())
    coords_shifted[:, 2] -= force_origin_z
    return coords_shifted


def save_stretched_geometry_yaml(
    input_path: Path,
    coords_final: np.ndarray,
    force_node_indices: Iterable[int],
    *,
    output_dir: Path = RESULTS_DIR,
) -> Path:
    """Copy a geometry YAML and replace only nodal coordinates."""
    coords_output = coords_shifted_to_force_origin(coords_final, force_node_indices)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{input_path.stem}_stretched{input_path.suffix}"

    bridle_point_re = re.compile(r"^(\s*bridle_point_node:\s*)\[[^\]]*\](.*)$")
    node_row_re = re.compile(r"^(\s*-\s*)\[\s*([+-]?\d+)\s*,[^\]]*\](.*)$")
    top_level_key_re = re.compile(r"^([A-Za-z_][A-Za-z0-9_]*)\s*:")
    node_tables = {"wing_particles", "bridle_particles"}

    active_node_table: str | None = None
    updated_bridle_point = False
    updated_node_ids: set[int] = set()
    output_lines: list[str] = []

    for line in input_path.read_text(encoding="utf-8").splitlines(keepends=True):
        body, ending = _split_line_ending(line)

        top_level_key = top_level_key_re.match(body)
        if top_level_key is not None:
            key = top_level_key.group(1)
            active_node_table = key if key in node_tables else None

        bridle_point_match = bridle_point_re.match(body)
        if bridle_point_match is not None:
            if len(coords_output) == 0:
                raise ValueError("Cannot update bridle_point_node without node 0.")
            body = (
                f"{bridle_point_match.group(1)}"
                f"[{_format_yaml_position(coords_output[0])}]"
                f"{bridle_point_match.group(2)}"
            )
            updated_bridle_point = True
        elif active_node_table is not None:
            node_match = node_row_re.match(body)
            if node_match is not None:
                node_id = int(node_match.group(2))
                if node_id < 0 or node_id >= len(coords_output):
                    raise ValueError(
                        f"Node {node_id} in {input_path} is out of bounds for "
                        f"{len(coords_output)} stretched coordinates."
                    )
                body = (
                    f"{node_match.group(1)}"
                    f"[{node_id}, {_format_yaml_position(coords_output[node_id])}]"
                    f"{node_match.group(3)}"
                )
                updated_node_ids.add(node_id)

        output_lines.append(body + ending)

    if not updated_bridle_point:
        raise ValueError(f"No bridle_point_node entry found in {input_path}.")
    if not updated_node_ids:
        raise ValueError(f"No node coordinate rows updated in {input_path}.")

    output_path.write_text("".join(output_lines), encoding="utf-8")
    print(f"[stretching_lines] wrote {output_path}")
    return output_path


def save_stretched_geometry_yamls(
    psm_result: dict[str, Any],
    fem_result: dict[str, Any],
) -> None:
    save_stretched_geometry_yaml(
        PSM_GEOMETRY_PATH,
        psm_result["coords_final"],
        psm_result["origin_nodes"],
    )
    save_stretched_geometry_yaml(
        FEM_GEOMETRY_PATH,
        fem_result["coords_final"],
        fem_result["origin_nodes"],
    )


def fem_segments(structure: FEM_structure) -> list[tuple[int, int, str]]:
    segments = []
    for beam in structure.beam_elements:
        segments.append((int(beam.beam.n1), int(beam.beam.n2), "beam"))
    for spring in structure.spring_elements:
        segments.append(
            (int(spring.spring.n1), int(spring.spring.n2), spring.springtype)
        )
    return segments


def set_equal_3d_limits(ax, coords: np.ndarray) -> None:
    coords = np.asarray(coords, dtype=float)
    finite = np.isfinite(coords).all(axis=1)
    pts = coords[finite]
    if len(pts) == 0:
        return
    mins = pts.min(axis=0)
    maxs = pts.max(axis=0)
    center = 0.5 * (mins + maxs)
    half_range = 0.55 * float(np.max(maxs - mins))
    if half_range == 0:
        half_range = 1.0
    ax.set_xlim(center[0] - half_range, center[0] + half_range)
    ax.set_ylim(center[1] - half_range, center[1] + half_range)
    ax.set_zlim(center[2] - half_range, center[2] + half_range)
    ax.set_box_aspect((1.0, 1.0, 1.0))


def clean_axis(ax) -> None:
    ax.grid(True, alpha=0.35)
    ax.set_axis_on()
    ax.view_init(elev=0, azim=-180)


def plot_segments(ax, coords, segments, wing_nodes: set[int]) -> None:
    for n1, n2, kind in segments:
        if n1 < 0 or n2 < 0 or n1 >= len(coords) or n2 >= len(coords):
            continue
        p1 = coords[n1]
        p2 = coords[n2]
        if not np.isfinite([*p1, *p2]).all():
            continue

        if kind == "beam":
            color, lw, zorder = "black", 1.4, 4
        elif n1 in wing_nodes and n2 in wing_nodes:
            color, lw, zorder = "#0072B2", 0.8, 3
        elif kind == "pulley":
            color, lw, zorder = "darkgreen", 0.8, 5
        else:
            color, lw, zorder = "0.45", 0.7, 2

        ax.plot(
            [p1[0], p2[0]],
            [p1[1], p2[1]],
            [p1[2], p2[2]],
            color=color,
            lw=lw,
            zorder=zorder,
        )


def plot_model(ax, title, coords, segments, wing_nodes, origin_nodes) -> None:
    wing_nodes = set(int(idx) for idx in wing_nodes)
    origin_nodes = [int(idx) for idx in origin_nodes]

    plot_segments(ax, coords, segments, wing_nodes)

    wing = [idx for idx in wing_nodes if 0 <= idx < len(coords)]
    if wing:
        pts = coords[wing]
        ax.scatter(
            pts[:, 0],
            pts[:, 1],
            pts[:, 2],
            s=4,
            color="black",
            depthshade=False,
            zorder=8,
        )

    origin = [
        idx
        for idx in origin_nodes
        if 0 <= idx < len(coords) and np.isfinite(coords[idx]).all()
    ]
    if origin:
        pts = coords[origin]
        ax.scatter(
            pts[:, 0],
            pts[:, 1],
            pts[:, 2],
            s=24,
            color="#D55E00",
            depthshade=False,
            zorder=12,
        )

    set_equal_3d_limits(ax, coords)
    clean_axis(ax)
    ax.set_title(title, fontsize=10)


def psm_segments_for_plot(psm_result) -> list[tuple[int, int, str]]:
    psm_wing_nodes = set(psm_result["wing_nodes"])
    return [
        (
            n1,
            n2,
            "wing" if n1 in psm_wing_nodes and n2 in psm_wing_nodes else "bridle",
        )
        for n1, n2 in psm_result["connectivity"]
    ]


def show_final_3d_output(psm_result, fem_result) -> None:
    psm_segments = psm_segments_for_plot(psm_result)
    fem_plot_segments = fem_segments(fem_result["structure"])

    fig = plt.figure(figsize=(10.0, 4.8))
    axes = [
        fig.add_subplot(1, 2, 1, projection="3d"),
        fig.add_subplot(1, 2, 2, projection="3d"),
    ]

    plot_model(
        axes[0],
        "PSM final",
        psm_result["coords_final"],
        psm_segments,
        psm_result["wing_nodes"],
        psm_result["origin_nodes"],
    )
    plot_model(
        axes[1],
        "FEM final",
        fem_result["coords_final"],
        fem_plot_segments,
        fem_result["wing_nodes"],
        fem_result["origin_nodes"],
    )

    fig.tight_layout(pad=0.8)
    plt.show()
    plt.close(fig)


def save_plot(psm_result, fem_result, show: bool = False) -> None:
    psm_segments = psm_segments_for_plot(psm_result)
    fem_plot_segments = fem_segments(fem_result["structure"])

    fig = plt.figure(figsize=(9.5, 7.0))
    axes = [
        fig.add_subplot(2, 2, 1, projection="3d"),
        fig.add_subplot(2, 2, 2, projection="3d"),
        fig.add_subplot(2, 2, 3, projection="3d"),
        fig.add_subplot(2, 2, 4, projection="3d"),
    ]

    plot_model(
        axes[0],
        "PSM before stretching",
        psm_result["coords_initial"],
        psm_segments,
        psm_result["wing_nodes"],
        psm_result["origin_nodes"],
    )
    plot_model(
        axes[1],
        "PSM after stretching",
        psm_result["coords_final"],
        psm_segments,
        psm_result["wing_nodes"],
        psm_result["origin_nodes"],
    )
    plot_model(
        axes[2],
        "FEM before stretching",
        fem_result["coords_initial"],
        fem_plot_segments,
        fem_result["wing_nodes"],
        fem_result["origin_nodes"],
    )
    plot_model(
        axes[3],
        "FEM after stretching",
        fem_result["coords_final"],
        fem_plot_segments,
        fem_result["wing_nodes"],
        fem_result["origin_nodes"],
    )

    fig.tight_layout(pad=0.6)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    view_outputs = [
        (-180, OUTPUT_PDF_AZ_180, OUTPUT_SVG_AZ_180),
        (-90, OUTPUT_PDF_AZ_90, OUTPUT_SVG_AZ_90),
    ]

    for azim, pdf_path, svg_path in view_outputs:
        for ax in axes:
            ax.view_init(elev=0, azim=azim)
        fig.savefig(pdf_path, bbox_inches="tight")
        fig.savefig(svg_path, bbox_inches="tight")
        print(f"[stretching_lines] wrote {pdf_path}")
        print(f"[stretching_lines] wrote {svg_path}")

    if show:
        plt.show()
    plt.close(fig)


def main(show: bool = False) -> None:
    psm_result = run_psm_stretch()
    fem_result = run_fem_stretch()
    write_markdown_report(psm_result, fem_result)
    save_stretched_geometry_yamls(psm_result, fem_result)
    save_plot(psm_result, fem_result)
    if show or SHOW_FINAL_3D_OUTPUT:
        show_final_3d_output(psm_result, fem_result)


if __name__ == "__main__":
    main(show="--show" in {arg.lower() for arg in sys.argv[1:]})
