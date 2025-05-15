from pydantic import BaseModel
from typing import List, Dict, Union


class Airfoil(BaseModel):
    tube_diameters: List[float]
    is_tube_diameter_dimensionless: bool
    canopy_max_heights: List[float]
    is_canopy_max_height_dimensionless: bool


class Connectivity(BaseModel):
    bridle_ci: List[int]
    bridle_cj: List[int]
    plate_point_indices: List[List[int]]
    wing_ci: List[int]
    wing_cj: List[int]
    te_line_indices: List[int]
    tube_line_indices: List[int]


class Kcu(BaseModel):
    drag_coefficient: float
    diameter: float
    index: int
    mass: float


class Pulley(BaseModel):
    point_indices: List[int]
    mass: float
    number_of_pulleys_in_back_lines: List[int]
    line_indices: List[int]
    line_pair_indices: List[List[str]]
    ci: List[int]
    cj: List[int]
    other_line_pair: List[Union[str, List[float]]]


class Bridle(BaseModel):
    diameter: float
    density: float
    bridle_point_index: int
    depower_tape_index: int
    left_steering_tape_index: int
    right_steering_tape_index: int


class Kite(BaseModel):
    points_ini: List[List[float]]
    n_points: int
    surfplan_filename: str
    area_projected: float
    area_surface: float
    ref_chord: float
    span: float
    height: float
    wing_mass: float
    is_with_elongation_limit: bool
    elongation_limit: float
    is_with_compression_limit: bool
    compression_limit: float
    limit_stiffness_factor: float
    billowing_angles: List[float]
    n_segments: int
    wing_rest_lengths_initial: List[float]
    bridle_rest_lengths_initial: List[float]
    mass_points: List[float]
    bridle: Bridle
    pulley: Pulley
    kcu: Kcu
    connectivity: Connectivity
    airfoil: Airfoil
    stiffness: Dict[str, float]


class Tether(BaseModel):
    diameter: float
    length: float
    density: float


class AeroStructural(BaseModel):
    max_iter: int
    it_check: int
    tol: float
    crosswind_max_iter: int
    crosswind_tol: float
    crosswind_relax_factor: float


class Solver(BaseModel):
    method: str
    tol: float
    max_fev: int
    integration_method: str
    newton_rtol: float
    newton_max_iter: int
    newton_disp: bool
    damping_constant: float
    is_with_visc_damping: bool
    dt: float
    n_time_steps: int
    abs_tol: float
    rel_tol: float
    max_iter: int


class Aero(BaseModel):
    ring_geometry: str
    model: str
    n_iter: int
    error_limit: float
    relax_factor: float
    n_splits: int
    n_chordwise_aero_nodes: int
    cd_cylinder: float
    cd_shear_cylinder: float
    cd_multiplier: float


class ConfigData(BaseModel):
    kite_name: str
    sim_name: str
    u_p: float
    depower_tape_extension_percentage: float
    delta_ls_ini: float
    depower_tape_extension_step: float
    depower_tape_final_extension: float
    steering_tape_extension_step: float
    steering_tape_final_extension: float
    is_with_printing: bool
    is_print_mid_results: bool
    is_with_initial_plot: bool
    is_with_initial_point_velocity: bool
    is_with_plotly_plot: bool
    is_with_aero_geometry: bool
    is_with_plotting: bool
    plot_format: str
    plot_elev: List[float]
    plot_azim: List[float]
    is_with_animation: bool
    animation_fps: int
    animation_dpi: int
    animation_bitrate: int
    animation_elev: float
    animation_azim: float
    is_run_only_1_time_step: bool
    bridle_initial_compression_factor: float
    geometric_scaling_factor: float
    is_with_velocity_initialization: bool
    vel_app_initial: List[float]
    n_vel_initialisation_steps: int
    is_billowing_on: bool
    is_with_aero_bridle: bool
    is_with_aero_tether: bool
    coupling_method: str
    grav_constant: List[float]
    rho: float
    mu: float
    aero: Aero
    solver: Solver
    aero_structural: AeroStructural
    tether: Tether
    tol_fx_ratio_to_fz: float
    tol_vk_optimization: float
    vk_x_initial_guess_factor_of_vw: float
    is_with_varying_va: bool
    r_0_initial: float
    is_with_vk_optimization: bool
    is_circular_case: bool
    is_print_intermediate_results: bool
    is_with_gravity: bool
    vel_wind: List[float]
    vel_kite: List[float]
    acc_kite: List[float]
    kite: Kite
