import numpy as np
import logging
import os
from pathlib import Path
from VSM.WingGeometry import Wing
from VSM.WingAerodynamics import WingAerodynamics
from VSM.Solver import Solver
from VSM.plotting import plot_polars, plot_distribution, plot_geometry

# Find the root directory of the repository
root_dir = os.path.abspath(os.path.dirname(__file__))
while not os.path.isfile(os.path.join(root_dir, ".gitignore")):
    root_dir = os.path.abspath(os.path.join(root_dir, ".."))
    if root_dir == "/":
        raise FileNotFoundError("Could not find the root directory of the repository.")
save_results_folder = Path(root_dir) / "examples" / "VSM" / "results"

# Defining settings
n_panels = 18
spanwise_panel_distribution = "split_provided"

### rib_list_from_CAD_LE_TE_and_surfplan_d_tube_camber_10ribs
csv_file_path = (
    Path(root_dir)
    / "data"
    / "V3_25"
    / "rib_list_from_CAD_LE_TE_and_surfplan_d_tube_camber_10ribs.csv"
)
(
    LE_x_array,
    LE_y_array,
    LE_z_array,
    TE_x_array,
    TE_y_array,
    TE_z_array,
    d_tube_array,
    camber_array,
) = np.loadtxt(csv_file_path, delimiter=",", skiprows=1, unpack=True)
rib_list_from_CAD_LE_TE_and_surfplan_d_tube_camber_10ribs = []
for i in range(len(LE_x_array)):
    LE = np.array([LE_x_array[i], LE_y_array[i], LE_z_array[i]])
    TE = np.array([TE_x_array[i], TE_y_array[i], TE_z_array[i]])
    print(f"i, {i}, LE: {LE}, TE: {TE}")
    rib_list_from_CAD_LE_TE_and_surfplan_d_tube_camber_10ribs.append(
        [LE, TE, ["lei_airfoil_breukels", [d_tube_array[i], camber_array[i]]]]
    )
CAD_wing = Wing(n_panels, spanwise_panel_distribution)

for i, CAD_rib_i in enumerate(
    rib_list_from_CAD_LE_TE_and_surfplan_d_tube_camber_10ribs
):
    CAD_wing.add_section(CAD_rib_i[0], CAD_rib_i[1], CAD_rib_i[2])
wing_aero_CAD_10ribs = WingAerodynamics([CAD_wing])

### Setting va for each wing
aoa_rad = np.deg2rad(10)
side_slip = 0
yaw_rate = 0
Umag = 12
va = (
    np.array(
        [
            np.cos(aoa_rad) * np.cos(side_slip),
            np.sin(side_slip),
            np.sin(aoa_rad),
        ]
    )
    * Umag
)
va = np.array([12, 0, 2])
wing_aero_CAD_10ribs.va = va, yaw_rate


# printing the wing geometry
for panel in wing_aero_CAD_10ribs.panels:
    print(f"panel ac: {panel.aerodynamic_center}")
    print(f"panel cp: {panel.control_point}")

### plotting shapes
plot_geometry(
    wing_aero_CAD_10ribs,
    title="wing_aero_CAD_10ribs",
    data_type=".pdf",
    save_path=Path(save_results_folder) / "geometry",
    is_save=True,
    is_show=True,
    view_elevation=15,
    view_azimuth=-120,
)


## Running the solver
solver = Solver(aerodynamic_model_type="VSM", is_with_artificial_damping=True)
results = solver.solve(wing_aero_CAD_10ribs)


def is_symmetric_1d(array, tol=1e-8):
    return np.allclose(array, array[::-1], atol=tol)


print(f"\n VSM is symmetric: {is_symmetric_1d(results['gamma_distribution'])}")
