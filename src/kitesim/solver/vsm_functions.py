from pathlib import Path
import yaml
import numpy as np
import pandas as pd
from VSM.BodyAerodynamics import BodyAerodynamics
from VSM.WingGeometry import Wing
from VSM.Solver import Solver
from VSM.interactive import interactive_plot
from kitesim.utils import load_yaml


def initialize_vsm(
    config_kite, n_panels: int, spanwise_panel_distribution: str
) -> BodyAerodynamics:
    """
    Load kite configuration from a YAML file and initialize the VSM BodyAerodynamics
    object with one Wing instance constructed from the 'airfoils' data.

    Args:
        config_path (Path): Path to the YAML configuration file.
        n_panels (int): Number of panels for the wing.
        spanwise_panel_distribution (str): Type of spanwise distribution.

    Returns:
        BodyAerodynamics: Initialized body aerodynamic model.
        Solver: Initialized solver object.
    """

    # 1) Extract airfoil table
    af = config_kite["airfoils"]
    headers = af["headers"]  # e.g. ["LE_x", "LE_y", ...]
    data = af["data"]  # list of lists

    # 2) Build DataFrame for convenience
    df = pd.DataFrame(data, columns=headers)

    # 3) Create Wing instance
    wing = Wing(
        n_panels=n_panels, spanwise_panel_distribution=spanwise_panel_distribution
    )

    # 4) Add each airfoil section to the wing
    for _, row in df.iterrows():
        LE = np.array([row["LE_x"], row["LE_y"], row["LE_z"]])
        TE = np.array([row["TE_x"], row["TE_y"], row["TE_z"]])
        airfoil_data = ["lei_airfoil_breukels", [row["d_tube"], row["y_camber"]]]
        wing.add_section(LE, TE, airfoil_data)

    return BodyAerodynamics([wing]), Solver()


def plot_vsm_geometry(body_aero):
    # angle_of_attack = np.rad2deg(np.arctan(vel_app[2] / vel_app[0]))
    # interactive_plot(
    #     body_aero,
    #     np.linalg.norm(vel_app),
    #     angle_of_attack=angle_of_attack,
    #     side_slip=0,
    #     yaw_rate=0,
    #     title="Interactive plot",
    #     is_with_aerodynamic_details=False,
    #     save_path=None,
    #     is_save=False,
    #     filename="wing_geometry",
    #     is_show=True,
    # )
    from VSM.plotting import plot_geometry

    plot_geometry(
        body_aero,
        title="VSM Geometry",
        data_type=None,
        save_path=None,
        is_save=False,
        is_show=True,
        view_elevation=15,
        view_azimuth=-120,
    )


# def initialize_vsm(
#     geometry_csv_path,
#     polar_data_dir,
#     n_panels,
#     spanwise_panel_distribution="uniform",
#     is_half_wing=True,
#     is_with_corrected_polar=True,
# ):
#     """
#     Initialize VSM simulation components.

#     Returns:
#         (body_aero, solver): Instantiated aerodynamic body and solver objects.
#     """
#     body_aero = BodyAerodynamics.from_file(
#         geometry_csv_path,
#         n_panels=n_panels,
#         spanwise_panel_distribution=spanwise_panel_distribution,
#         is_with_corrected_polar=is_with_corrected_polar,
#         polar_data_dir=polar_data_dir,
#         is_half_wing=is_half_wing,
#     )

#     solver = Solver()
#     return body_aero, solver


def run_vsm_package(
    body_aero,
    solver,
    le_arr,
    te_arr,
    va_vector,
    aero_input_type="reuse_initial_polar_data",
    initial_polar_data=None,
    yaw_rate=0.0,
):
    """
    Run the VSM simulation with updated geometry and velocity.

    Returns:
        force_distribution: Nx3 array of aerodynamic force vectors.
    """
    # redefine the points
    body_aero.update_from_points(
        le_arr,
        te_arr,
        aero_input_type=aero_input_type,
        initial_polar_data=initial_polar_data,
    )
    # set again where velocity vector is coming from
    body_aero.va = (va_vector, yaw_rate)
    # solve the problem
    results = solver.solve(body_aero)
    return np.array(results["F_distribution"]), body_aero, results
