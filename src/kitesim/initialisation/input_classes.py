import numpy as np
from attrs import frozen
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
