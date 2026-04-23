import numpy as np
import pandas as pd
from pathlib import Path
import copy
from VSM.core.BodyAerodynamics import BodyAerodynamics
from VSM.core.WingGeometry import Wing
from VSM.core.Solver import Solver
from VSM.plot_geometry_matplotlib import plot_geometry
from VSM.quasi_steady_state import (
    solve_quasi_steady_state,
    DEFAULT_TRANSFORMATION_C_FROM_VSM,
)


# Bounds and defaults (aoa, sideslip, course_rate_body)
kite_speed_bounds = (1.0, 50.0)  # m/s
pitch_bounds = (-4, 4)  # deg
yaw_bounds = (-4, 4)  # deg
course_rate_bounds = (
    -3,
    3,
)  # rad/s, small course rate allowed for numerical reasons; not a physical course rate
roll_bounds = (
    -4,
    4,
)  # deg, small roll allowed for numerical reasons; not a physical roll

DEFAULT_GUESS_QS = np.array(
    [30.0, 0.0, 0.0, 0.0, -0.0]
)  # [kite_speed, roll, pitch, yaw, course_rate_body]


def _run_vsm_direct_fallback(body_aero, solver, system_model, current_guess):
    """
    Fallback used when quasi-steady trim fails.

    This path bypasses trim optimization and runs a direct aerodynamic solve on the
    current body geometry. The returned dictionary mirrors the keys expected by the
    coupled solver pipeline.
    """
    res = solver.solve(body_aero)

    # Convert AWETrim vectors into the same course-frame convention used by QSM outputs.
    trans = np.asarray(DEFAULT_TRANSFORMATION_C_FROM_VSM, dtype=float)
    inertial_force = -float(system_model.mass_wing) * np.asarray(
        trans @ np.asarray(system_model.acceleration_course_body, dtype=float),
        dtype=float,
    ).reshape(3)
    gravity_force = np.asarray(
        trans @ np.asarray(system_model.force_gravity, dtype=float),
        dtype=float,
    ).reshape(3)

    guess = np.asarray(current_guess, dtype=float).reshape(5)
    fallback_opt_x = np.array(DEFAULT_GUESS_QS, dtype=float)
    print(
        f"Falling back to direct VSM solve with guess {guess} (quasi-steady optimization failed)."
    )
    print(
        f"Direct VSM results: {res.get('F_distribution')}, cmx: {res.get('cmx')}, cmy: {res.get('cmy')}, cmz: {res.get('cmz')}, side_slip_deg: {res.get('side_slip_deg')}, side_slip_course_deg: {res.get('side_slip_course_deg')}"
    )
    return {
        "opt_x": fallback_opt_x,
        "success": False,
        "inertial_force": inertial_force,
        "gravity_force": gravity_force,
        "panel_cp_locations": res.get("panel_cp_locations"),
        "F_distribution": res.get("F_distribution"),
        "alpha_at_ac": res.get("alpha_at_ac"),
    }


def initialize(
    aero_geometry_path,
    config,
    n_panels_aero: int,
    bridle_path=None,
) -> BodyAerodynamics:
    """
    Initialize aerodynamic model and VSM solver.

    Args:
        aero_geometry_path: Path to aerodynamic geometry file.
        config (dict): Main ASKITE configuration dictionary.
        n_panels_aero (int): Number of aerodynamic panels.
        bridle_path: Optional structural geometry path used by VSM to build bridle lines.

    Returns:
        tuple: (body_aero, vsm_solver, vel_app, initial_polar_data)
    """
    body_aero = BodyAerodynamics.instantiate(
        n_panels=int(n_panels_aero),
        file_path=aero_geometry_path,
        spanwise_panel_distribution=config["aerodynamic"][
            "spanwise_panel_distribution"
        ],
        bridle_path=bridle_path,
    )

    vsm_solver = Solver(
        max_iterations=config["aerodynamic"]["max_iterations"],
        allowed_error=config["aerodynamic"]["allowed_error"],
        relaxation_factor=config["aerodynamic"]["relaxation_factor"],
        reference_point=config["aerodynamic"]["reference_point"],
        mu=config["mu"],
        rho=config["rho"],
    )

    # For QSM, wind speed comes from system model configuration (wind_speed_wind_ref).
    # Kite velocity is computed by the optimizer, so we initialize with wind direction.
    wind_speed_ref = float(config.get("wind_speed_wind_ref", 6.0))
    vel_app = np.array([wind_speed_ref, 0.0, 0.0])
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
    include_gravity=False,
    is_with_plot=False,
    current_guess=None,
):
    """
    Run quasi-steady aerodynamic solve for the current structural geometry.

    Args:
        body_aero (BodyAerodynamics): Aerodynamic body object.
        solver (Solver): VSM solver object.
        system_model: AWETrim system model used by quasi-steady trim.
        center_of_gravity (np.ndarray): Current center of gravity in solver frame.
        le_arr (np.ndarray): Leading edge points (n,3).
        te_arr (np.ndarray): Trailing edge points (n,3).
        aero_input_type (str): Type of aerodynamic input.
        initial_polar_data (list or None): Initial polar data for panels.
        reference_point (list[float]): Reference point for moments and rotations.
        include_gravity (bool): Include gravity in quasi-steady force/moment balance.
        is_with_plot (bool): If True, plot the geometry.
        current_guess (np.ndarray or None): Initial guess for quasi-steady optimizer.

    Returns:
        tuple: (F_distribution, body_aero, results)
    """
    # Update aerodynamic mesh from the latest structural leading/trailing-edge points.
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

    # Primary path: quasi-steady trim solve.
    try:
        # Preserve the pre-trim body state in case we need direct-solve fallback.
        body_fallback = copy.deepcopy(body_aero)
        results, body_aero = solve_quasi_steady_state(
            body_aero=body_aero,
            center_of_gravity=center_of_gravity,
            reference_point=reference_point,
            system_model=system_model,
            x_guess=current_guess,
            solver=solver,
            bounds_lower=bounds_lower,
            bounds_upper=bounds_upper,
            include_gravity=include_gravity,
        )
        if not results.get("success", False):
            print(
                "Quasi-steady optimization did not converge to a valid trim state. "
                "Falling back to direct VSM solver.solve(body_aero)."
            )
            results = _run_vsm_direct_fallback(
                body_aero=body_fallback,
                solver=solver,
                system_model=system_model,
                current_guess=current_guess,
            )
    except ValueError as exc:
        # Typical case: non-finite residual in initial optimizer point.
        print(
            f"QSM failed ({type(exc).__name__}: {exc}). "
            "Falling back to direct VSM solver.solve(body_aero)."
        )
        results = _run_vsm_direct_fallback(
            body_aero=body_aero,
            solver=solver,
            system_model=system_model,
            current_guess=current_guess,
        )
    if is_with_plot:
        plot_vsm_geometry(body_aero)
    return np.array(results["F_distribution"]), body_aero, results
