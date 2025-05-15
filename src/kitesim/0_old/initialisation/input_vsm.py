from kitesim.initialisation import initialisation_utils


def create_input_VSM(config):
    """Create the input VSM class.

    Args:
        config: The configuration class.

    Returns:
        The input VSM class."""
    input_vsm_dict = {
        "model": config.aero.model,
        "n_splits": config.aero.n_splits,
        "billowing_angles": config.kite.billowing_angles,
        "ring_geometry": config.aero.ring_geometry,
        "n_segments": config.kite.n_segments,
        "tube_diameters": config.kite.airfoil.tube_diameters,
        "is_tube_diameter_dimensionless": config.kite.airfoil.is_tube_diameter_dimensionless,
        "canopy_max_heights": config.kite.airfoil.canopy_max_heights,
        "is_canopy_max_height_dimensionless": config.kite.airfoil.is_canopy_max_height_dimensionless,
        "cd_multiplier": config.aero.cd_multiplier,
        "n_iter": config.aero.n_iter,
        "error_limit": config.aero.error_limit,
        "relax_factor": config.aero.relax_factor,
    }
    return initialisation_utils.create_attr_class_from_dict("input_vsm", input_vsm_dict)
