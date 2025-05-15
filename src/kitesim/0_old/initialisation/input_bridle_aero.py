from kitesim.initialisation import initialisation_utils


def create_input_bridle_aero(config):
    """Create the input bridle aero class.

    Args:
        config: The configuration class.

    Returns:
        The input bridle aero class.
    """
    input_bridle_aero_dict = {
        "is_with_aero_bridle": config.is_with_aero_bridle,
        "vel_wind": config.vel_wind,
        "bridle_ci": config.kite.connectivity.bridle_ci,
        "bridle_cj": config.kite.connectivity.bridle_cj,
        "rho": config.rho,
        "cd_cylinder": config.aero.cd_cylinder,
        "cd_shear_cylinder": config.aero.cd_shear_cylinder,
        "diameter": config.kite.bridle.diameter,
        "kcu_cd": config.kite.kcu.drag_coefficient,
        "kcu_index": config.kite.kcu.index,
        "kcu_diameter": config.kite.kcu.diameter,
    }
    return initialisation_utils.create_attr_class_from_dict(
        "input_bridle_aero", input_bridle_aero_dict
    )
