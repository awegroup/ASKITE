import numpy as np
import pandas as pd
from pathlib import Path
from VSM.core.BodyAerodynamics import BodyAerodynamics
from VSM.core.WingGeometry import Wing
from VSM.core.Solver import Solver
from VSM.plot_geometry_matplotlib import plot_geometry


def initialize(
    kite_name,
    PROJECT_DIR,
    config,
    n_panels_aero: int,
) -> BodyAerodynamics:
    """
    Load kite configuration and initialize the VSM BodyAerodynamics object with one Wing instance.

    Args:
        geometry_dict (dict): Kite configuration dictionary.
        n_panels (int): Number of panels for the wing.
        spanwise_panel_distribution (str): Type of spanwise distribution.

    Returns:
        BodyAerodynamics: Initialized body aerodynamic model.
        Solver: Initialized solver object.
    """
    body_aero = BodyAerodynamics.instantiate(
        n_panels=int(n_panels_aero),
        file_path=(Path(PROJECT_DIR) / "data" / f"{kite_name}" / "aero_geometry.yaml"),
        spanwise_panel_distribution=config["aerodynamic"][
            "spanwise_panel_distribution"
        ],
        # is_with_bridles=False,
    )

    vsm_solver = Solver(
        max_iterations=config["aerodynamic"]["max_iterations"],
        allowed_error=config["aerodynamic"]["allowed_error"],
        relaxation_factor=config["aerodynamic"]["relaxation_factor"],
        reference_point=config["aerodynamic"]["reference_point"],
        mu=config["mu"],
        rho=config["rho"],
    )

    vel_app = np.array(config["vel_wind"]) - np.array(config["vel_kite"])
    body_aero.va = (vel_app, 0)
    wing = body_aero.wings[0]
    new_sections = wing.refine_aerodynamic_mesh()
    initial_polar_data = []
    for new_section in new_sections:
        initial_polar_data.append(new_section.polar_data)

    return body_aero, vsm_solver, vel_app, initial_polar_data


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
