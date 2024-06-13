def calculate_tether_drag(input_tether_aero, vel_app):
    ## TETHER DRAG (commented is pauls method, uncommented is the mentioned paper method)
    # unit_vector_bridle = np.array([0,0,1])
    # vel_bridle_abs = np.linalg.norm(vel_app)
    # vel_bridle_unit = vel_app/vel_bridle_abs
    # BRIDLE_CD = 1.1
    # BRIDLE_CD_SHEAR = 0.02
    # BRIDLE_DIAMETER = 0.04
    # length_bridle = 200
    # # angle between the bridle and the velocity vector
    # angle_bridle_vel = np.arccos(np.dot(unit_vector_bridle,vel_bridle_unit))

    # # getting the forces
    # lift_coefficient_bridle = BRIDLE_CD * (np.sin(angle_bridle_vel)**2) * (np.cos(angle_bridle_vel))
    # drag_coefficient_bridle = BRIDLE_CD * (np.sin(angle_bridle_vel)**3) + BRIDLE_CD_SHEAR

    # lift_bridle = 0.5 * RHO * (vel_bridle_abs**2) * lift_coefficient_bridle * (BRIDLE_DIAMETER * 0.5* length_bridle)
    # drag_bridle = 0.5 * RHO * (vel_bridle_abs**2) * drag_coefficient_bridle * (BRIDLE_DIAMETER * 0.5*length_bridle)
    # print(f'lift:{lift_bridle}, drag:{drag_bridle}')

    # drag_direction_unit = vel_bridle_unit
    # lift_direction = VEL_WIND - np.dot(VEL_WIND, drag_direction_unit) * drag_direction_unit
    # lift_direction_unit = lift_direction / np.linalg.norm(lift_direction)
    # drag_tether = lift_bridle*lift_direction_unit+drag_bridle*drag_direction_unit
    # print(f'Tether drag: {drag_tether}')

    # from this paper: https://www.researchgate.net/publication/329390819_A_reference_model_for_airborne_wind_energy_systems_for_optimization_and_control
    if input_tether_aero.is_with_aero_tether:
        drag_tether = (
            (1 / 8)
            * input_tether_aero.rho
            * input_tether_aero.cd_cylinder
            * input_tether_aero.tether_diameter
            * input_tether_aero.tether_length
            * (vel_app[0] ** 2)
        )
    else:
        drag_tether = 0
    return drag_tether
