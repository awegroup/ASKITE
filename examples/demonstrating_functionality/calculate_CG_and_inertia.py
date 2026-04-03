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
from kitesim.utils import calculate_cg, calculate_moments_of_inertia, load_yaml


def _resolve_level_1_geometry_path(project_dir):
    """Resolve TUDELFT_V3_KITE level-1 geometry file path."""
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


def main():
    struc_geometry_path = _resolve_level_1_geometry_path(PROJECT_DIR)
    struc_geometry = load_yaml(struc_geometry_path)

    geometry_data = read_struc_geometry_yaml_level_1.main(struc_geometry)
    struc_nodes = np.asarray(geometry_data[0], dtype=float)
    m_arr = np.asarray(geometry_data[1], dtype=float)

    cg = calculate_cg(struc_nodes=struc_nodes, m_arr=m_arr)
    inertia_origin = calculate_moments_of_inertia(struc_nodes=struc_nodes, m_arr=m_arr)
    inertia_cg = calculate_moments_of_inertia(
        struc_nodes=struc_nodes,
        m_arr=m_arr,
        point=cg,
    )

    Ixx_o, Iyy_o, Izz_o = inertia_origin
    Ixx_cg, Iyy_cg, Izz_cg = inertia_cg

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
    print("Moments of Inertia about Origin [kg*m^2]")
    print(f"  Ixx: {Ixx_o:.6f}")
    print(f"  Iyy: {Iyy_o:.6f}")
    print(f"  Izz: {Izz_o:.6f}")
    print("")
    print("Moments of Inertia about CG [kg*m^2]")
    print(f"  Ixx: {Ixx_cg:.6f}")
    print(f"  Iyy: {Iyy_cg:.6f}")
    print(f"  Izz: {Izz_cg:.6f}")


if __name__ == "__main__":
    main()
