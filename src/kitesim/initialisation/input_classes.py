import numpy as np
from attrs import define, frozen, field, Factory
from typing import Any

## Defining the immutable input classes seperately
# Such that it uses less arguments to speed up the calculations.


@frozen
class InputVSM:
    model: str
    n_splits: int
    billowing_angles: np.ndarray
    ring_geometry: str
    n_segments: int
    tube_diameters: np.ndarray
    is_tube_diameter_dimensionless: bool
    canopy_max_heights: np.ndarray
    is_canopy_max_height_dimensionless: bool
    cd_multiplier: float
    n_iter: int
    error_limit: float
    relax_factor: float

    @classmethod
    def create(cls, config: Any) -> "InputVSM":
        return cls(
            model=config.aero.model,
            n_splits=config.aero.n_splits,
            billowing_angles=config.kite.billowing_angles,
            ring_geometry=config.aero.ring_geometry,
            n_segments=config.kite.n_segments,
            tube_diameters=config.kite.airfoil.tube_diameters,
            is_tube_diameter_dimensionless=config.kite.airfoil.is_tube_diameter_dimensionless,
            canopy_max_heights=config.kite.airfoil.canopy_max_heights,
            is_canopy_max_height_dimensionless=config.kite.airfoil.is_canopy_max_height_dimensionless,
            cd_multiplier=config.aero.cd_multiplier,
            n_iter=config.aero.n_iter,
            error_limit=config.aero.error_limit,
            relax_factor=config.aero.relax_factor,
        )


@frozen
class InputBridleAero:
    vel_wind: np.ndarray
    bridle_ci: np.ndarray
    bridle_cj: np.ndarray
    rho: float
    cd_cylinder: float
    cd_shear_cylinder: float
    diameter: float
    kcu_cd: float
    kcu_index: int
    kcu_diameter: float

    @classmethod
    def create(cls, config: Any) -> "InputBridleAero":
        return cls(
            vel_wind=config.vel_wind,
            bridle_ci=config.kite.connectivity.bridle_ci,
            bridle_cj=config.kite.connectivity.bridle_cj,
            rho=config.rho,
            cd_cylinder=config.aero.cd_cylinder,
            cd_shear_cylinder=config.aero.cd_shear_cylinder,
            diameter=config.kite.bridle.diameter,
            kcu_cd=config.kite.kcu.drag_coefficient,
            kcu_index=config.kite.kcu.index,
            kcu_diameter=config.kite.kcu.diameter,
        )


# TODO: IS NOT USED?

# @frozen
# class InputCalculateForceSpring:
#     bridle_ci: np.ndarray
#     bridle_cj: np.ndarray
#     wing_ci: np.ndarray
#     wing_cj: np.ndarray
#     te_line_indices: np.ndarray
#     tube_line_indices: np.ndarray
#     bridle_stiffness: float
#     tube_stiffness: float
#     te_stiffness: float
#     canopy_stiffness: float
#     pulley_line_indices: np.ndarray
#     pulley_line_pair_indices: np.ndarray
#     is_with_elongation_limit: bool
#     elongation_limit: float
#     is_with_compression_limit: bool
#     compression_limit: float
#     limit_stiffness_factor: float


# @frozen
# class InputCalculateTotalForce:
#     bridle_point_index: int
#     acc_kite: int
#     is_with_gravity: bool
#     mass_points: np.ndarray
#     INPUT_calculate_force_spring: InputCalculateForceSpring


# @frozen
# class InputStructuralSolver:
#     method: str
#     max_fev: int
#     newton_max_fev: int

#     INPUT_calculate_total_force: InputCalculateTotalForce


# input_structural_solver = InputStructuralSolver(
#     method=config.solver.method,
#     max_fev=config.solver.max_fev,
#     newton_max_fev=config.solver.newton_max_iter,
#     INPUT_calculate_total_force=InputCalculateTotalForce(
#         bridle_point_index=config.kite.bridle.bridle_point_index,
#         acc_kite=config.acc_kite,
#         is_with_gravity=config.is_with_gravity,
#         mass_points=config.kite.mass_points,
#         INPUT_calculate_force_spring=InputCalculateForceSpring(
#             bridle_ci=config.kite.connectivity.bridle_ci,
#             bridle_cj=config.kite.connectivity.bridle_cj,
#             wing_ci=config.kite.connectivity.wing_ci,
#             wing_cj=config.kite.connectivity.wing_cj,
#             te_line_indices=config.kite.connectivity.te_line_indices,
#             tube_line_indices=config.kite.connectivity.tube_line_indices,
#             bridle_stiffness=config.kite.stiffness.bridle,
#             tube_stiffness=config.kite.stiffness.tube,
#             te_stiffness=config.kite.stiffness.trailing_edge,
#             canopy_stiffness=config.kite.stiffness.canopy,
#             pulley_line_indices=config.kite.pulley.line_indices,
#             pulley_line_pair_indices=config.kite.pulley.line_pair_indices,
#             is_with_elongation_limit=config.kite.is_with_elongation_limit,
#             elongation_limit=config.kite.elongation_limit,
#             is_with_compression_limit=config.kite.is_with_compression_limit,
#             compression_limit=config.kite.compression_limit,
#             limit_stiffness_factor=config.kite.limit_stiffness_factor,
#         ),
#     ),
# )
