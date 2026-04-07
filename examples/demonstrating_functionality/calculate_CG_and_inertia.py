from pathlib import Path
import sys

import numpy as np

# Allow running this script directly from the repository root without
# requiring an editable install or PYTHONPATH setup.
PROJECT_DIR = Path(__file__).resolve().parents[2]
SRC_DIR = PROJECT_DIR / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from kitesim import read_struc_geometry_yaml_level_1
from kitesim.utils import calculate_cg, calculate_inertia, load_yaml


def _print_inertia_tensor(label, tensor):
    """Print all 9 inertia tensor components with a label."""
    print(f"{label} [kg*m^2]")
    print(f"  Ixx: {tensor[0, 0]:.6f}")
    print(f"  Ixy: {tensor[0, 1]:.6f}")
    print(f"  Ixz: {tensor[0, 2]:.6f}")
    print(f"  Iyx: {tensor[1, 0]:.6f}")
    print(f"  Iyy: {tensor[1, 1]:.6f}")
    print(f"  Iyz: {tensor[1, 2]:.6f}")
    print(f"  Izx: {tensor[2, 0]:.6f}")
    print(f"  Izy: {tensor[2, 1]:.6f}")
    print(f"  Izz: {tensor[2, 2]:.6f}")


def _build_wing_structural_mass_arr(struc_geometry, n_nodes):
    """
    Build nodal masses from wing structure only (wing_elements + wing_connections).

    This excludes all bridle-related mass contributions.
    """
    wing_mass_arr = np.zeros(n_nodes, dtype=float)

    wing_elements_headers = {
        str(name): idx for idx, name in enumerate(struc_geometry["wing_elements"]["headers"])
    }
    wing_connections_headers = {
        str(name): idx
        for idx, name in enumerate(struc_geometry["wing_connections"]["headers"])
    }

    we_name_idx = wing_elements_headers.get("name", 0)
    we_mass_idx = wing_elements_headers.get("m", 4)
    wc_name_idx = wing_connections_headers.get("name", 0)
    wc_ci_idx = wing_connections_headers.get("ci", 1)
    wc_cj_idx = wing_connections_headers.get("cj", 2)

    # Same name-based lookup convention as ASKITE level-1 geometry parser.
    wing_element_mass_by_name = {
        str(row[we_name_idx]): float(row[we_mass_idx])
        for row in struc_geometry["wing_elements"]["data"]
    }

    for row in struc_geometry["wing_connections"]["data"]:
        conn_name = str(row[wc_name_idx])
        ci = int(row[wc_ci_idx])
        cj = int(row[wc_cj_idx])
        if conn_name not in wing_element_mass_by_name:
            raise KeyError(
                f"Connection '{conn_name}' not found in wing_elements mass table."
            )
        m_element = wing_element_mass_by_name[conn_name]
        wing_mass_arr[ci] += 0.5 * m_element
        wing_mass_arr[cj] += 0.5 * m_element

    return wing_mass_arr


def _resolve_level_1_geometry_path(project_dir):
    """Resolve TUDELFT_V3_KITE level-1 geometry file path."""
    candidate_paths = [
        # project_dir
        # / "data"
        # / "TUDELFT_V3_KITE"
        # / "struc_geometry_level_1_manual_JULIA.yaml",
        project_dir
        / "data"
        / "TUDELFT_V3_KITE"
        / "struc_geometry_level_1_manual.yaml",
    ]
    for path in candidate_paths:
        if path.exists():
            return path
    raise FileNotFoundError(
        "Could not find a level-1 geometry YAML. Checked:\n"
        + "\n".join(str(p) for p in candidate_paths)
    )


def main():
    struc_geometry_path = _resolve_level_1_geometry_path(PROJECT_DIR)
    struc_geometry = load_yaml(struc_geometry_path)

    geometry_data = read_struc_geometry_yaml_level_1.main(struc_geometry)
    struc_nodes = np.asarray(geometry_data[0], dtype=float)
    m_arr = np.asarray(geometry_data[1], dtype=float)
    point_mass_nodes = [(struc_nodes[i], m_arr[i]) for i in range(len(struc_nodes))]
    wing_node_indices = sorted(
        {int(ci) for _, ci, _ in struc_geometry["wing_connections"]["data"]}
        | {int(cj) for _, _, cj in struc_geometry["wing_connections"]["data"]}
    )
    wing_structural_mass_arr = _build_wing_structural_mass_arr(
        struc_geometry, len(struc_nodes)
    )
    wing_point_mass_nodes = [
        (struc_nodes[i], wing_structural_mass_arr[i]) for i in wing_node_indices
    ]
    wing_nodes_arr = np.asarray([struc_nodes[i] for i in wing_node_indices], dtype=float)
    wing_masses_arr = np.asarray(
        [wing_structural_mass_arr[i] for i in wing_node_indices], dtype=float
    )

    cg = calculate_cg(struc_nodes=struc_nodes, m_arr=m_arr)
    wing_cg = calculate_cg(struc_nodes=wing_nodes_arr, m_arr=wing_masses_arr)
    inertia_origin = calculate_inertia(point_mass_nodes, desired_point=(0.0, 0.0, 0.0))
    inertia_cg = calculate_inertia(point_mass_nodes, desired_point=cg)
    wing_only_inertia_cg = calculate_inertia(
        wing_point_mass_nodes, desired_point=wing_cg
    )

    print("TUDELFT_V3_KITE level-1 mass properties")
    print(f"  geometry file: {struc_geometry_path}")
    print(f"  n_nodes: {len(struc_nodes)}")
    print(f"  total mass [kg]: {np.sum(m_arr):.6f}")
    print("")
    print("Center of Gravity [m]")
    print(f"  x_cg: {cg[0]:.6f}")
    print(f"  y_cg: {cg[1]:.6f}")
    print(f"  z_cg: {cg[2]:.6f}")
    print("")
    _print_inertia_tensor("Inertia Tensor about Origin", inertia_origin)
    print("")
    _print_inertia_tensor("Inertia Tensor about CG", inertia_cg)
    print("")
    print(
        "Wing-only subset (structural only): "
        f"{len(wing_node_indices)} nodes, mass [kg]: {np.sum([m for _, m in wing_point_mass_nodes]):.6f}"
    )
    print("Wing-only CG [m]")
    print(f"  x_cg_wing: {wing_cg[0]:.6f}")
    print(f"  y_cg_wing: {wing_cg[1]:.6f}")
    print(f"  z_cg_wing: {wing_cg[2]:.6f}")
    _print_inertia_tensor(
        "Inertia Tensor about Wing-only CG (wing nodes only)",
        wing_only_inertia_cg,
    )


if __name__ == "__main__":
    main()
