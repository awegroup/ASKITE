from src.initialisation.yaml_loader import config
from src.initialisation import (
    input_particleSystem,
    particles_with_rotational_resistance,
)
from src.post_processing import functions_plot
from src.particleSystem.ParticleSystem import ParticleSystem


def get_mutable_variables():

    ## Mutable variables
    # Initializing Mutable Variables
    points = config.kite.points_ini
    # defining vel_app (vector of vel_app_norm)
    vel_app = config.vel_wind - config.vel_kite

    # Should be the same for each kite
    (
        connectivity_matrix,
        wing_connectivity,
    ) = input_particleSystem.define_connectivity_matrix(config)
    params_dict = input_particleSystem.define_params(
        config, wing_connectivity, connectivity_matrix
    )
    initial_conditions = input_particleSystem.define_initial_conditions_kite(config)
    points_between_dict = (
        particles_with_rotational_resistance.extract_points_between_dict(config)
    )
    if config.is_with_initial_plot:
        functions_plot.plot_initial_geometry(config, points_between_dict)

    is_with_rotational_resistance = False
    if config.kite_name == "V9_60C":
        is_with_rotational_resistance = True

    if is_with_rotational_resistance:
        (
            leadingedge_rotational_resistance_dict,
            strut_rotational_resistance_dict,
        ) = particles_with_rotational_resistance.extract_rotational_resistances_dicts(
            points_between_dict, config
        )
        # first do the struts
        params_dict = particles_with_rotational_resistance.initialize_bending_spring(
            config.kite.stiffness.k_bend_strut,
            initial_conditions,
            params_dict,
            connectivity_matrix,
            strut_rotational_resistance_dict,
        )
        # secondly do the leading-edge
        params_dict = particles_with_rotational_resistance.initialize_bending_spring(
            config.kite.stiffness.k_bend_leading_edge,
            initial_conditions,
            params_dict,
            connectivity_matrix,
            leadingedge_rotational_resistance_dict,
        )
    # Should be the same for each kite
    psystem = ParticleSystem(connectivity_matrix, initial_conditions, params_dict)

    return (
        points,
        vel_app,
        params_dict,
        psystem,
    )
