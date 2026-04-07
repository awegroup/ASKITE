import numpy as np
import pandas as pd
from pathlib import Path
from VSM.core.BodyAerodynamics import BodyAerodynamics
from VSM.core.WingGeometry import Wing
from VSM.core.Solver import Solver
from VSM.plot_geometry_matplotlib import plot_geometry
from VSM.quasi_steady_state import solve_quasi_steady_state


# Bounds and defaults (aoa, sideslip, course_rate_body)
kite_speed_bounds = (2.0, 80.0)  # m/s
pitch_bounds = (-5, 5)  # deg
yaw_bounds = (-6, 6)  # deg
course_rate_bounds = (
    -3,
    3,
)  # rad/s, small course rate allowed for numerical reasons; not a physical course rate
roll_bounds = (
    -5,
    5,
)  # deg, small roll allowed for numerical reasons; not a physical roll

DEFAULT_GUESS_QS = np.array(
    [40.0, 0.0, 0.0, 0.0, 0.0]
)  # [kite_speed, roll, pitch, yaw, course_rate_body


def initialize(
    aero_geometry_path,
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
        file_path=aero_geometry_path,
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
    body_aero.va = vel_app
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
    system_model,
    center_of_gravity,
    le_arr,
    te_arr,
    # va_vector,
    aero_input_type="reuse_initial_polar_data",
    initial_polar_data=None,
    reference_point=[0.0, 0.0, 0.0],
    is_with_plot=False,
    current_guess=None,
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
    # The VSM va setter accepts keyword arguments but properties don't support that in Python
    # So we call the underlying setter method directly using the descriptor protocol
    # type(body_aero).va.fset(body_aero, va_vector)

    bounds_lower = np.array(
        [
            kite_speed_bounds[0],
            roll_bounds[0],
            pitch_bounds[0],
            yaw_bounds[0],
            course_rate_bounds[0],
        ]
    )
    bounds_upper = np.array(
        [
            kite_speed_bounds[1],
            roll_bounds[1],
            pitch_bounds[1],
            yaw_bounds[1],
            course_rate_bounds[1],
        ]
    )
    if current_guess is None:
        current_guess = DEFAULT_GUESS_QS
    # solve the problem
    results, body_aero = solve_quasi_steady_state(
        body_aero=body_aero,
        center_of_gravity=center_of_gravity,
        reference_point=reference_point,
        system_model=system_model,
        x_guess=current_guess,
        solver=solver,
        bounds_lower=bounds_lower,
        bounds_upper=bounds_upper,
    )
    if is_with_plot:
        plot_vsm_geometry(body_aero)
    return np.array(results["F_distribution"]), body_aero, results
