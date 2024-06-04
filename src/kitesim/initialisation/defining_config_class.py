from dataclasses import dataclass, astuple
import numpy as np
from attr import attrs, define, frozen

# attrs has slots as a default, meaning there is only a limited #attributes -immutable
# can also make with with @define(slot=False), to allow you to add attributes dynamically in-code!
# dataclasses doesn't have slot as a default


# TODO: figure out how to do it like this, instead of writing boilerplate code
# # child-classes, inherent the from_yaml loader functionality.
# class BaseConfig:
#     def create_from_dict(self, config_data: dict):
#         """Find keys in the yaml file and assign them to the class attributes."""
#         for key, value in config_data.items():
#             if hasattr(self, key):
#                 setattr(self, key, value)
#             else:
#                 raise AttributeError(
#                     f"Attribute {key} not found in class {self.__class__.__name__})"
#                 )


# Child Classes
@frozen
class AeroConfig:
    """
    Configuration settings related to aerodynamics.

    Attributes:
        ring_geometry (str): The geometry of the ring.
        model (str): The aerodynamic model.
        n_iter (int): Number of iterations.
        error_limit (float): Error limit.
        relax_factor (float): Relaxation factor.
        n_splits (int): Number of splits.
        bridle_cd (float): Bridle coefficient of drag.
        bridle_cd_shear (float): Bridle coefficient of drag for shear.
        cd_multiplier (float): Coefficient of drag multiplier.
    """

    ring_geometry: str
    model: str
    n_iter: int
    error_limit: float
    # used to be: 1e-5
    relax_factor: float
    n_splits: int
    n_chordwise_aero_nodes: int
    cd_cylinder: float
    cd_shear_cylinder: float
    cd_multiplier: float


@frozen
class SolverConfig:
    """
    Configuration settings related to the solver.

    Attributes:
        method (str): The solver method ('root' or 'newton').
        tol (float): Tolerance in N (tol < calculated_total_force).
        max_fev (int): Maximum number of function evaluations.
        integration_method (str): Integration method.
        newton_rtol (float): Newton solver relative tolerance.
        newton_max_iter (int): Maximum number of Newton iterations.
        newton_disp (bool): Whether to display Newton solver output.
    """

    method: str  # method: either 'root' or 'newton'
    tol: float  # tolerance in N, tol < calculated_total_force

    # root settings
    max_fev: int
    integration_method: str

    # newton settings
    newton_rtol: float
    newton_max_iter: int
    newton_disp: bool

    # alex his psm-setting
    damping_constant: float
    is_with_visc_damping: bool
    dt: float
    n_time_steps: int
    abs_tol: float
    rel_tol: float
    max_iter: int


@frozen
class AeroStructuralConfig:
    """
    Configuration settings related to aerostructural analysis.

    Attributes:
        max_iter (int): Maximum number of iterations.
        it_check (int): Check interval for iterations.
        tol (float): Tolerance.
        crosswind_max_iter (int): Maximum number of iterations for crosswind flight.
        crosswind_tol (float): Tolerance for crosswind flight.
        crosswind_relax_factor (float): Relaxation factor for crosswind flight.
    """

    max_iter: int
    it_check: int
    tol: float

    # # TODO: implement it like this? (then also remove above)
    # ## symmetric
    # symmetric_max_iter: int
    # symmetric_tol_iter_diff: float
    # symmetric_tol_solver: float
    # symmetric_it_check: float

    ## crosswind flight
    crosswind_max_iter: int
    # crosswind_tol_iter_diff: float
    crosswind_tol: float
    # crosswind_it_check: float
    # crosswind_tol_fx_ratio_to_fz: float
    crosswind_relax_factor: float

    # ## circular flight
    # circular_max_iter: int
    # circular_tol_iter_diff: float
    # circular_tol_solver: float
    # circular_it_check: float
    # circular_tol_fx_ratio_to_fz: float


@frozen
class TetherConfig:
    """
    Configuration settings related to the tether.

    Attributes:
        diameter (float): Diameter of the tether.
        length (float): Length of the tether.
        density (float): Density of the tether material.
    """

    diameter: float
    length: float
    density: float


@frozen
class BridleConfig:
    diameter: float
    density: float
    bridle_point_index: int
    depower_tape_index: int
    left_steering_tape_index: int
    right_steering_tape_index: int


@frozen
class PulleyConfig:
    point_indices: np.ndarray
    mass: float
    number_of_pulleys_in_back_lines: int
    line_indices: np.ndarray
    line_pair_indices: np.ndarray
    ci: np.ndarray
    cj: np.ndarray
    other_line_pair: dict

    # def __array__(self):
    #     return np.array(astuple(self))

    # def __len__(self):
    #     return astuple(self).__len__()

    # def __getitem__(self, item):
    #     return astuple(self).__getitem__(item)


@frozen
class KCUConfig:
    drag_coefficient: float
    diameter: float
    index: int
    mass: float


@frozen
class ConnectivityConfig:
    bridle_ci: np.ndarray
    bridle_cj: np.ndarray
    plate_point_indices: np.ndarray
    wing_ci: np.ndarray
    wing_cj: np.ndarray
    te_line_indices: np.ndarray
    tube_line_indices: np.ndarray


@frozen
class AirfoilGeometry:
    tube_diameters: np.ndarray
    is_tube_diameter_dimensionless: bool
    canopy_max_heights: np.ndarray
    is_canopy_max_height_dimensionless: bool


@frozen
class StiffnessConfig:
    bridle: float
    tube: float
    trailing_edge: float
    canopy: float
    # rotational
    k_bend_strut: float
    k_bend_leading_edge: float


@frozen
class KiteConfig:
    """
    Configuration settings related to the kite.

    Attributes:
        surfplan_filename (str): Filename for the kite's surfplan.
        area_projected (float): Projected area of the kite.
        ref_chord (float): Reference chord length.
        wing_mass (float): Mass of the kite.
        stiffness_bridle (float): Stiffness of the bridle.
        stiffness_tube (float): Stiffness of the tube.
        stiffness_trailing_edge (float): Stiffness of the trailing edge.
        stiffness_canopy (float): Stiffness of the canopy.
        is_with_elongation_limit (bool): Whether elongation limit is enabled.
        elongation_limit (float): Elongation limit.
        is_with_compression_limit (bool): Whether compression limit is enabled.
        compression_limit (float): Compression limit.
        limit_stiffness_factor (float): Stiffness multiplier if over the limit.
    """

    points_ini: np.ndarray
    n_points: int
    surfplan_filename: str

    # geometry
    area_surface: float
    area_projected: float
    ref_chord: float
    span: float
    height: float
    wing_mass: float

    # structural settings
    is_with_elongation_limit: bool
    elongation_limit: float
    is_with_compression_limit: bool
    compression_limit: float
    # the factor by which stiffness is increased if over the limit
    limit_stiffness_factor: float

    # billowing angles
    billowing_angles: np.ndarray

    # number of structural segments
    n_segments: int

    # initial rest lengths
    wing_rest_lengths_initial: np.ndarray
    bridle_rest_lengths_initial: np.ndarray

    # mass distributed over all the nodes
    mass_points: np.ndarray

    ## Child classes
    bridle: BridleConfig
    pulley: PulleyConfig
    kcu: KCUConfig
    connectivity: ConnectivityConfig
    airfoil: AirfoilGeometry
    stiffness: StiffnessConfig

    def __array__(self):
        return np.array(astuple(self))

    def __len__(self):
        return astuple(self).__len__()

    def __getitem__(self, item):
        return astuple(self).__getitem__(item)


# Parent Class
@frozen
class Config:
    """
    Configuration settings for the entire simulation.

    Attributes:
        kite_name (str): Name of the kite.
        vel_wind (np.ndarray): Wind velocity vector.
        vel_kite (np.ndarray): Kite velocity vector.
        acc_kite (np.ndarray): Kite acceleration vector.
        u_p (float): Parameter U_P.
        depower_tape_extension_percentage (float): Depower tape extension percentage.
        is_print_mid_results (bool): Whether to print intermediate results.
        is_with_initial_plot (bool): Whether to plot initial geometry.
        is_billowing_on (bool): Whether billowing is enabled.
        is_with_gravity (bool): Whether gravity is enabled.
        is_with_aero_bridle (bool): Whether aero bridle is enabled.
        bridle_initial_compression_factor (float): Initial compression factor for bridle.
        grav_constant (float): Gravitational constant.
        rho (float): Density.
        mu (float): Coefficient of viscosity.

    Children Classses:
        aero (AeroConfig): Aero configuration settings.
        solver (SolverConfig): Solver configuration settings.
        aero_structural (AeroStructuralConfig): Aerostructural configuration settings.
        tether (TetherConfig): Tether configuration settings.
        kite (KiteConfig): Kite configuration settings.
    """

    # KITE NAME
    kite_name: str

    # CASE SETTINGS
    sim_name: str
    is_with_vk_optimization: bool
    is_circular_case: bool
    is_run_only_1_time_step: bool
    is_print_intermediate_results: bool
    is_with_gravity: bool
    is_with_velocity_initialization: bool

    # INFLOW CONDITIONS
    vel_wind: np.ndarray
    vel_kite: np.ndarray
    acc_kite: np.ndarray

    # ACTUATION
    ## depower
    u_p: float
    depower_tape_extension_percentage: float
    ## steering
    delta_ls_ini: float
    ## new in meters
    depower_tape_extension_step: float
    depower_tape_final_extension: float
    steering_tape_extension_step: float
    steering_tape_final_extension: float

    # OUTPUT SETTINGS
    is_with_printing: bool
    is_with_plotting: bool
    is_with_animation: bool
    is_print_mid_results: bool
    is_with_initial_plot: bool
    is_with_initial_point_velocity: bool
    is_with_plotly_plot: bool
    is_with_aero_geometry: bool

    # SIMULATION SETTINGS
    ## initialisation
    bridle_initial_compression_factor: float
    geometric_scaling_factor: float
    n_vel_initialisation_steps: int
    ## physics
    is_billowing_on: bool
    is_with_aero_bridle: bool
    is_with_aero_tether: bool
    coupling_method: str

    # CASES
    ## reel-in symmetric straight
    vel_app_initial: np.ndarray
    ## crosswind flight settings
    tol_fx_ratio_to_fz: float
    tol_vk_optimization: float
    vk_x_initial_guess_factor_of_vw: float
    ## circular flight settings
    is_with_varying_va: bool
    r_0_initial: 200.0

    # PHYSICAL CONSTANTS
    grav_constant: np.ndarray
    rho: float
    mu: float

    # CHILD CLASSES
    aero: AeroConfig
    solver: SolverConfig
    aero_structural: AeroStructuralConfig
    tether: TetherConfig
    kite: KiteConfig

    def __array__(self):
        return np.array(astuple(self))

    def __len__(self):
        return astuple(self).__len__()

    def __getitem__(self, item):
        return astuple(self).__getitem__(item)
