from kitesim.initialisation import initialisation_utils


def create_tether_aero(config):
    """Create the input bridle aero class.

    Args:
        config: The configuration class.

    Returns:
        The input bridle aero class.
    """
    input_tether_aero_dict = {
        "is_with_aero_tether": config.is_with_aero_tether,
        "diameter": config.tether.diameter,
        "length": config.tether.length,
        "rho": config.rho,
        "cd_cylinder": config.aero.cd_cylinder,
    }
    return initialisation_utils.create_attr_class_from_dict(
        "input_tether_aero", input_tether_aero_dict
    )
