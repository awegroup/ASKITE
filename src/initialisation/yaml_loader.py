from dataclasses import dataclass, astuple
import numpy as np
import yaml
import os
import importlib

from attrs import define, frozen, field, Factory
from src.structural import structural_mesher
from src.initialisation.path_functions import load_module_from_path
from src.initialisation import mass_distribution
from scipy.spatial import ConvexHull

# attrs has slots as a default, meaning there is only a limited #attributes -immutable
# can also make with with @define(slot=False), to allow you to add attributes dynamically in-code!
# dataclasses doesn't have slot as a default


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
        is_from_filename (bool): Whether to load from a filename.
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
    is_with_save: bool
    is_print_mid_results: bool
    is_with_initial_plot: bool
    is_with_initial_point_velocity: bool
    is_with_plotly_plot: bool
    is_with_aero_geometry: bool
    output_path: str

    # SIMULATION SETTINGS
    ## initialisation
    is_from_filename: bool
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


# Loading the yaml config file
with open(("settings/config.yaml"), "r") as config_file:
    config_data = yaml.load(config_file, Loader=yaml.SafeLoader)

## Initialize child-classes first
# Initialize Aero configuration
aero_config = AeroConfig(
    ring_geometry=config_data["AeroConfig"]["ring_geometry"],
    model=config_data["AeroConfig"]["model"],
    n_iter=config_data["AeroConfig"]["n_iter"],
    error_limit=config_data["AeroConfig"]["error_limit"],
    relax_factor=config_data["AeroConfig"]["relax_factor"],
    n_splits=config_data["AeroConfig"]["n_splits"],
    n_chordwise_aero_nodes=config_data["AeroConfig"]["n_chordwise_aero_nodes"],
    cd_cylinder=config_data["AeroConfig"]["cd_cylinder"],
    cd_shear_cylinder=config_data["AeroConfig"]["cd_shear_cylinder"],
    cd_multiplier=config_data["AeroConfig"]["cd_multiplier"],
)

solver_config = SolverConfig(
    method=config_data["SolverConfig"]["method"],
    tol=config_data["SolverConfig"]["tol"],
    max_fev=config_data["SolverConfig"]["max_fev"],
    integration_method=config_data["SolverConfig"]["integration_method"],
    newton_rtol=config_data["SolverConfig"]["newton_rtol"],
    newton_max_iter=config_data["SolverConfig"]["newton_max_iter"],
    newton_disp=config_data["SolverConfig"]["newton_disp"],
    damping_constant=config_data["SolverConfig"]["damping_constant"],
    is_with_visc_damping=config_data["SolverConfig"]["is_with_visc_damping"],
    dt=config_data["SolverConfig"]["dt"],
    n_time_steps=config_data["SolverConfig"]["n_time_steps"],
    abs_tol=config_data["SolverConfig"]["abs_tol"],
    rel_tol=config_data["SolverConfig"]["rel_tol"],
    max_iter=config_data["SolverConfig"]["max_iter"],
)

aero_structural_config = AeroStructuralConfig(
    max_iter=config_data["AeroStructuralConfig"]["max_iter"],
    it_check=config_data["AeroStructuralConfig"]["it_check"],
    tol=config_data["AeroStructuralConfig"]["tol"],
    crosswind_max_iter=config_data["AeroStructuralConfig"]["crosswind_max_iter"],
    crosswind_tol=config_data["AeroStructuralConfig"]["crosswind_tol"],
    crosswind_relax_factor=config_data["AeroStructuralConfig"][
        "crosswind_relax_factor"
    ],
)

tether_config = TetherConfig(
    diameter=config_data["TetherConfig"]["diameter"],
    length=config_data["TetherConfig"]["length"],
    density=config_data["TetherConfig"]["density"],
)

## Initialize Kite-specific configuration

## Finding the Connectivity
# Reading out surfplan
# TODO: implement real surfplan-reader here surfplan_file = load_surfplan.fetch_surfplan_file(config.kite.surfplan_filename,sys.path) # reading out the surfplan file
surfplan_file = None

kite_name = config_data["kite_name"]

# Defining paths
folder_path_kite = f"data/input/{kite_name}"
folder_path_kite_data = f"{folder_path_kite}/processed_design_files"
folder_path_initialisation = f"src/initialisation/{kite_name}"
folder_path_output = f"data/output/{kite_name}"

extract_points_and_connectivity = load_module_from_path(
    kite_name, f"{folder_path_initialisation}/analyze_surfplan.py"
).extract_points_and_connectivity

(
    points_struc,
    bridle_ci,
    bridle_cj,
    plate_point_indices,
    wing_ci,
    wing_cj,
    te_line_indices,
    tube_line_indices,
) = extract_points_and_connectivity(folder_path_kite_data, surfplan_file)

# points
if config_data["is_from_filename"]:
    # if there is a file defined
    points_ini = np.load(
        f'{folder_path_output}/points/points_up_{int(config_data["u_p"]*100)}.npy'
    )
else:
    points_ini = np.copy(points_struc)  # if no initial file is defined, use surfplan
points_ini = np.array(points_ini) / config_data["geometric_scaling_factor"]

bridle_rest_lengths_initial = np.array(
    structural_mesher.calculate_edge_lengths(bridle_ci, bridle_cj, points_ini)
)
wing_rest_lengths_initial = np.array(
    structural_mesher.calculate_edge_lengths(wing_ci, wing_cj, points_ini)
)  # [m]

# Billowing (both TE and middle TE points to make it symmetrical), when config.is_billowing_on is True
if config_data["is_billowing_on"]:
    points_ini, wing_rest_lengths = structural_mesher.update_for_billowing(
        points_ini, wing_rest_lengths_initial, config_data["u_p"]
    )


# number of structural segments
n_segments = int(len(plate_point_indices))


### opening up the kite specific config file
with open((f"{folder_path_kite}/config_kite_{kite_name}.yaml"), "r") as config_file:
    config_data_kite = yaml.load(config_file, Loader=yaml.SafeLoader)


### Getting additional data
if kite_name == "V3_25":
    # TODO: Most variables here should be generated/imported from the surfplan file instead
    BILLOWING_ANGLES = (
        np.array(config_data_kite["billowing_angles"])
        / config_data["geometric_scaling_factor"]
    )
    TUBE_DIAMETERS = (
        np.array(config_data_kite["tube_diameters"])
        / config_data["geometric_scaling_factor"]
    )
    CANOPY_MAX_HEIGHTS = (
        np.array(config_data_kite["canopy_max_heights"])
        / config_data["geometric_scaling_factor"]
    )

    ## aditional data for bridles
    bridle_data = {}
    bridle_data["diameter"] = config_data_kite["bridle"]["diameter"]
    bridle_data["density"] = config_data_kite["bridle"]["density"]
    bridle_data["bridle_point_index"] = config_data_kite["bridle"]["bridle_point_index"]
    # element indices
    bridle_data["depower_tape_index"] = config_data_kite["bridle"]["depower_tape_index"]
    bridle_data["left_steering_tape_index"] = config_data_kite["bridle"][
        "left_steering_tape_index"
    ]
    bridle_data["right_steering_tape_index"] = config_data_kite["bridle"][
        "right_steering_tape_index"
    ]

    ## additional data for pulleys
    pulley_data = {}
    pulley_data["point_indices"] = config_data_kite["pulley"]["point_indices"]
    pulley_data["mass"] = config_data_kite["pulley"]["mass"]
    pulley_data["number_of_pulleys_in_back_lines"] = config_data_kite["pulley"][
        "number_of_pulleys_in_back_lines"
    ]

    ## kcu data
    kcu_data = {}
    kcu_data["index"] = config_data_kite["kcu"]["index"]
    kcu_data["extra"] = {}

    ## stiffness rotational ##TODO: can this be done cleaner?
    stiffness_bend_strut = 0
    stiffness_bend_leading_edge = 0

elif kite_name == "V9_60C":
    # TODO: Should be generated/imported directly from the surfplan file instead

    RIB_DB_WHOLE_STRUTS = np.load(
        f"{folder_path_kite_data}/rib_db_whole_model_struts.npy",
        allow_pickle=True,
    )
    extract_airfoil_geometry_data_from_ribs = load_module_from_path(
        kite_name, f"{folder_path_initialisation}/analyze_surfplan.py"
    ).extract_airfoil_geometry_data_from_ribs
    print(extract_airfoil_geometry_data_from_ribs)
    (
        TUBE_DIAMETERS,
        CANOPY_MAX_HEIGHTS,
        BILLOWING_ANGLES,
    ) = extract_airfoil_geometry_data_from_ribs(
        RIB_DB_WHOLE_STRUTS, folder_path_initialisation, config_data["kite_name"]
    )

    ## bridle data
    bridle_data = {}
    bridle_data["diameter"] = config_data_kite["bridle"]["diameter"]
    bridle_data["density"] = config_data_kite["bridle"]["density"]
    bridle_data["bridle_point_index"] = np.load(
        f"{folder_path_kite_data}/bridlepoint_index.npy", allow_pickle=True
    )
    # element indices
    bridle_data["depower_tape_index"] = config_data_kite["bridle"]["depower_tape_index"]
    bridle_data["left_steering_tape_index"] = config_data_kite["bridle"][
        "left_steering_tape_index"
    ]
    bridle_data["right_steering_tape_index"] = config_data_kite["bridle"][
        "right_steering_tape_index"
    ]

    ## pulley_data
    pulley_data = {}
    pulley_data["point_indices"] = np.load(
        f"{folder_path_kite_data}/pulley_point_indices.npy", allow_pickle=True
    )
    pulley_data["mass"] = config_data_kite["pulley"]["mass"]
    pulley_data["number_of_pulleys_in_back_lines"] = config_data_kite["pulley"][
        "number_of_pulleys_in_back_lines"
    ]

    ## kcu data
    kcu_data = {}
    kcu_data["extra"] = {
        "kcu_point_indices": np.load(
            f"{folder_path_kite_data}/kcu_point_indices.npy", allow_pickle=True
        ),
        "kcu_line_indices": np.load(
            f"{folder_path_kite_data}/kcu_line_indices.npy", allow_pickle=True
        ),
        "kcu_plate_indices": np.load(
            f"{folder_path_kite_data}/kcu_plate_indices.npy", allow_pickle=True
        ),
    }

    ## rotational stiffness
    stiffness_bend_strut = config_data_kite["stiffness_bend_strut"]
    stiffness_bend_leading_edge = config_data_kite["stiffness_bend_leading_edge"]


# get additional pulley_data
extract_pulley_connectivity = importlib.import_module(
    f"src.initialisation.pulley_connectivity"
).extract_pulley_connectivity
pulley_data = extract_pulley_connectivity(points_ini, bridle_ci, bridle_cj, pulley_data)


mass_points = mass_distribution.new_calculate_mass_distribution(
    points_ini,
    bridle_ci,
    bridle_cj,
    wing_ci,
    wing_cj,
    config_data_kite["wing_mass"],
    config_data_kite["bridle"]["density"],
    config_data_kite["bridle"]["diameter"],
    config_data_kite["kcu"]["mass"],
    config_data_kite["kcu"]["index"],
    pulley_data["point_indices"],
    config_data_kite["pulley"]["mass"],
)

## Calculating reference distance
ref_chord_calculated = max(points_ini[:, 0]) - min(points_ini[:, 0])


def calculate_projected_area(points):
    # Project points onto the x,y plane
    xy_points = points[:, :2]

    # Find the convex hull
    hull = ConvexHull(xy_points)
    hull_points = xy_points[hull.vertices]

    # Using the shoelace formula
    x = hull_points[:, 0]
    y = hull_points[:, 1]

    return 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))


wing_connectivity = np.hstack((wing_ci, wing_cj))
wing_connectivity = np.unique(wing_connectivity)
wing_nodes = np.array([points_ini[i] for i in wing_connectivity])
area_projected_calculated = calculate_projected_area(wing_nodes)
area_surface_calculated = (
    config_data_kite["area_surface"] / config_data["geometric_scaling_factor"]
)
span_calculated = max(wing_nodes[:, 1]) - min(wing_nodes[:, 1])
height_calculated = max(wing_nodes[:, 2]) - min(wing_nodes[:, 2])

kite_config = KiteConfig(
    points_ini=points_ini,
    n_points=int(len(points_ini)),
    surfplan_filename=config_data_kite["surfplan_filename"],
    area_projected=area_projected_calculated,
    area_surface=area_surface_calculated,
    ref_chord=ref_chord_calculated,
    span=span_calculated,
    height=height_calculated,
    wing_mass=config_data_kite["wing_mass"],
    is_with_elongation_limit=config_data_kite["is_with_elongation_limit"],
    elongation_limit=config_data_kite["elongation_limit"],
    is_with_compression_limit=config_data_kite["is_with_compression_limit"],
    compression_limit=config_data_kite["compression_limit"],
    limit_stiffness_factor=config_data_kite["limit_stiffness_factor"],
    billowing_angles=BILLOWING_ANGLES,
    n_segments=n_segments,
    wing_rest_lengths_initial=wing_rest_lengths_initial,
    bridle_rest_lengths_initial=bridle_rest_lengths_initial,
    mass_points=mass_points,
    ## Child classes
    bridle=BridleConfig(
        diameter=bridle_data["diameter"],
        density=bridle_data["density"],
        bridle_point_index=bridle_data["bridle_point_index"],
        depower_tape_index=bridle_data["depower_tape_index"],
        left_steering_tape_index=bridle_data["left_steering_tape_index"],
        right_steering_tape_index=bridle_data["right_steering_tape_index"],
    ),
    pulley=PulleyConfig(
        point_indices=np.array(pulley_data["point_indices"]),
        mass=np.array(pulley_data["mass"]),
        number_of_pulleys_in_back_lines=np.array(
            pulley_data["number_of_pulleys_in_back_lines"]
        ),
        line_indices=np.array(pulley_data["line_indices"]),
        line_pair_indices=np.array(pulley_data["line_pair_indices"]),
        ci=np.array(pulley_data["ci"]),
        cj=np.array(pulley_data["cj"]),
        other_line_pair=pulley_data["other_line_pair"],
    ),
    kcu=KCUConfig(
        drag_coefficient=config_data_kite["kcu"]["drag_coefficient"],
        diameter=config_data_kite["kcu"]["diameter"],
        index=config_data_kite["kcu"]["index"],
        mass=config_data_kite["kcu"]["mass"],
    ),
    connectivity=ConnectivityConfig(
        bridle_ci=bridle_ci,
        bridle_cj=bridle_cj,
        plate_point_indices=plate_point_indices,
        wing_ci=wing_ci,
        wing_cj=wing_cj,
        te_line_indices=te_line_indices,
        tube_line_indices=tube_line_indices,
    ),
    airfoil=AirfoilGeometry(
        tube_diameters=TUBE_DIAMETERS,
        is_tube_diameter_dimensionless=config_data_kite[
            "is_tube_diameter_dimensionless"
        ],
        canopy_max_heights=CANOPY_MAX_HEIGHTS,
        is_canopy_max_height_dimensionless=config_data_kite[
            "is_canopy_max_height_dimensionless"
        ],
    ),
    stiffness=StiffnessConfig(
        bridle=config_data_kite["stiffness_bridle"],
        tube=config_data_kite["stiffness_tube"],
        trailing_edge=config_data_kite["stiffness_trailing_edge"],
        canopy=config_data_kite["stiffness_canopy"],
        # rotational
        k_bend_strut=stiffness_bend_strut,
        k_bend_leading_edge=stiffness_bend_leading_edge,
    ),
)

# Initialize the overarching Config object
config = Config(
    # KITE NAME
    kite_name=config_data["kite_name"],
    # CASE SETTINGS
    sim_name=config_data["sim_name"],
    is_with_vk_optimization=config_data["is_with_vk_optimization"],
    is_circular_case=config_data["is_circular_case"],
    is_run_only_1_time_step=config_data["is_run_only_1_time_step"],
    is_print_intermediate_results=config_data["is_print_intermediate_results"],
    is_with_gravity=config_data["is_with_gravity"],
    is_with_velocity_initialization=config_data["is_with_velocity_initialization"],
    # INFLOW CONDITIONS
    vel_wind=np.array([float(i) for i in config_data["vel_wind"]]),
    vel_kite=np.array([float(i) for i in config_data["vel_kite"]]),
    acc_kite=np.array([float(i) for i in config_data["acc_kite"]]),
    # ACTUATION
    u_p=config_data["u_p"],
    depower_tape_extension_percentage=config_data["depower_tape_extension_percentage"],
    delta_ls_ini=config_data["delta_ls_ini"],
    depower_tape_extension_step=config_data["depower_tape_extension_step"],
    depower_tape_final_extension=config_data["depower_tape_final_extension"],
    steering_tape_extension_step=config_data["steering_tape_extension_step"],
    steering_tape_final_extension=config_data["steering_tape_final_extension"],
    # OUTPUT SETTINGS
    is_with_printing=config_data["is_with_printing"],
    is_with_plotting=config_data["is_with_plotting"],
    is_with_animation=config_data["is_with_animation"],
    is_with_save=config_data["is_with_save"],
    is_print_mid_results=config_data["is_print_mid_results"],
    is_with_initial_plot=config_data["is_with_initial_plot"],
    is_with_initial_point_velocity=config_data["is_with_initial_point_velocity"],
    is_with_plotly_plot=config_data["is_with_plotly_plot"],
    is_with_aero_geometry=config_data["is_with_aero_geometry"],
    output_path=folder_path_output,
    # SIMULATION SETTINGS
    ## initialisation
    is_from_filename=config_data["is_from_filename"],
    bridle_initial_compression_factor=config_data["bridle_initial_compression_factor"],
    geometric_scaling_factor=config_data["geometric_scaling_factor"],
    n_vel_initialisation_steps=config_data["n_vel_initialisation_steps"],
    ## physics
    is_billowing_on=config_data["is_billowing_on"],
    is_with_aero_bridle=config_data["is_with_aero_bridle"],
    is_with_aero_tether=config_data["is_with_aero_tether"],
    coupling_method=config_data["coupling_method"],
    # CASES
    ## reel-in symmetric straight
    vel_app_initial=np.array(config_data["vel_app_initial"]),
    ## crosswind flight settings
    tol_fx_ratio_to_fz=config_data["tol_fx_ratio_to_fz"],
    tol_vk_optimization=config_data["tol_vk_optimization"],
    vk_x_initial_guess_factor_of_vw=config_data["vk_x_initial_guess_factor_of_vw"],
    ## circular flight settings
    is_with_varying_va=config_data["is_with_varying_va"],
    r_0_initial=config_data["r_0_initial"],
    # PHYSICAL CONSTANTS
    grav_constant=np.array(config_data["grav_constant"]),
    rho=config_data["rho"],
    mu=config_data["mu"],
    # CHILD CLASSES
    aero=aero_config,
    solver=solver_config,
    aero_structural=aero_structural_config,
    tether=tether_config,
    kite=kite_config,
)

# @define
# class MutableVariables:
#     vel_app: np.ndarray = field(default=np.zeros(3))
#     bridle_rest_lengths: np.ndarray = field(default=np.zeros(len(bridle_ci)))
#     wing_rest_lengths: np.ndarray = field(default=np.zeros(len(wing_ci)))
#     mass_points: np.ndarray = field(default=np.zeros(len(points_ini)))
#     solver_max_error: float = field(default=1e2)

#     ## points
#     # ordered for aero
#     points_left_to_right: np.ndarray = field(default=np.zeros(len(points_ini)))
#     # initial
#     points_ini: np.ndarray = field(default=np.zeros(len(points_ini)))
#     # current points
#     points: np.ndarray = field(default=np.zeros(len(points_ini)))


# variables = MutableVariables()
