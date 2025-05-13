# VSM Initialisation and Execution Functions

from VSM.BodyAerodynamics import BodyAerodynamics
from VSM.Solver import Solver


def initialize_vsm(
    geometry_csv_path,
    polar_data_dir,
    n_panels,
    spanwise_panel_distribution="uniform",
    is_half_wing=True,
    is_with_corrected_polar=True,
):
    """
    Initialize VSM simulation components.

    Returns:
        (body_aero, solver): Instantiated aerodynamic body and solver objects.
    """
    body_aero = BodyAerodynamics.from_file(
        geometry_csv_path,
        n_panels=n_panels,
        spanwise_panel_distribution=spanwise_panel_distribution,
        is_with_corrected_polar=is_with_corrected_polar,
        polar_data_dir=polar_data_dir,
        is_half_wing=is_half_wing,
    )
    solver = Solver()
    return body_aero, solver


def run_vsm_package(
    body_aero,
    solver,
    le_arr,
    te_arr,
    va_vector,
    d_tube_arr=None,
    y_camber_arr=None,
    yaw_rate=0.0,
    aero_input_type="polar_data",
):
    """
    Run the VSM simulation with updated geometry and velocity.

    Returns:
        force_distribution: Nx3 array of aerodynamic force vectors.
    """
    body_aero.update_from_points(
        le_arr, te_arr, d_tube_arr, y_camber_arr, aero_input_type=aero_input_type
    )
    body_aero.va = (va_vector, yaw_rate)
    results = solver.solve(body_aero)
    return results["F_distribution"], body_aero
