import numpy as np
import pandas as pd
from VSM.BodyAerodynamics import BodyAerodynamics
from VSM.WingGeometry import Wing
from VSM.Solver import Solver
from VSM.plotting import plot_geometry


def initialize_vsm(
    config_kite, n_panels: int, spanwise_panel_distribution: str
) -> BodyAerodynamics:
    """
    Load kite configuration and initialize the VSM BodyAerodynamics object with one Wing instance.

    Args:
        config_kite (dict): Kite configuration dictionary.
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
    """
    Plot the VSM geometry using the provided aerodynamic body.

    Args:
        body_aero (BodyAerodynamics): Aerodynamic body object.

    Returns:
        None. Displays a 3D plot.
    """
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


def run_vsm_package(
    body_aero,
    solver,
    le_arr,
    te_arr,
    va_vector,
    aero_input_type="reuse_initial_polar_data",
    initial_polar_data=None,
    yaw_rate=0.0,
    is_with_plot=False,
):
    """
    Run the VSM simulation with updated geometry and velocity.

    Args:
        body_aero (BodyAerodynamics): Aerodynamic body object.
        solver (Solver): VSM solver object.
        le_arr (np.ndarray): Leading edge points (n,3).
        te_arr (np.ndarray): Trailing edge points (n,3).
        va_vector (np.ndarray): Apparent wind vector (3,).
        aero_input_type (str): Type of aerodynamic input.
        initial_polar_data (list or None): Initial polar data for panels.
        yaw_rate (float): Yaw rate for the simulation.
        is_with_plot (bool): If True, plot the geometry.

    Returns:
        np.ndarray: Aerodynamic force distribution (n_panels,3).
        BodyAerodynamics: Updated aerodynamic body.
        dict: Results dictionary from the solver.
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
    if is_with_plot:
        plot_vsm_geometry(body_aero)

    # solve the problem
    results = solver.solve(body_aero)
    return np.array(results["F_distribution"]), body_aero, results
