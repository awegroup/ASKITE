# %% Bridle Drag
## bridle drag

import numpy as np


def calculate_force_aero_bridle_geschiere2014(points, vel_app, input_bridle_aero):
    """--- Bridle aero ---
    Aerodynamic Lift force is neglected as cylinders only have substantial lifting force
    when rotating, which they don't really do [needs source]

    Aerodynamic Drag force (acts as a damping term)
        D = 0.5 * RHO * Vb^2 * Cd_b * Sb
    RHO is the air density [ISA]
    Vb is the velocity of the bridle = Vw (QSM)
    Sb is the ref-surface area of bridle
        "The reference area is typically defined as the area of the
        orthographic projection of the object on a plane perpendicular
        to the direction of motion.”
        Sb   = diam * effective length
        leff = 1 - e_s * Vb_a / ||Vb_a||  [Geschiere 2019]
            e_s  = unit vector in the direction of the bridle
            Vb_a = velocity of the bridle in the airframe frame
            Can be derived from:
            vec_vel * vec_bridle = ||vec_vel|| * ||vec_bridle|| * cos(theta)

    Cd_b is the drag coefficient = 1.2 ##TODO: verify this
        Comes from empirical CYLINDRICAL data [Horner 1965]
        For a Reynolds of approximately ~ 2E4 (up to 8e4, according to Geschiere2019)
            1.225*30*1E-3 / (18E-6) = 2E4
        Bridles ~ 1mm diameter (Geschiere says 2mm?) ##TODO: check this
        Vb      ~ 30m/s ##TODO: needs a good source
        RHO     ~ 1.225 kg/m^3
        mu      ~ 18E-6 kg/m/s [at 15C] ##TODO: needs a good source
    """

    vel_wind = input_bridle_aero.vel_wind
    bridle_ci = input_bridle_aero.bridle_ci
    bridle_cj = input_bridle_aero.bridle_cj
    RHO = input_bridle_aero.rho
    bridle_cd = input_bridle_aero.cd_cylinder
    bridle_cd_shear = input_bridle_aero.cd_shear_cylinder
    bridle_diameter = input_bridle_aero.diameter
    kcu_cd = input_bridle_aero.kcu_cd
    kcu_index = input_bridle_aero.kcu_index
    kcu_diameter = input_bridle_aero.kcu_diameter

    force_aero_bridle = np.zeros(points.shape)

    ## determining directions
    drag_direction_unit = vel_app / np.linalg.norm(vel_app)
    if (
        vel_wind.all() == vel_app.all()
    ):  # avoid division by 0, ##TODO: remove hardcoding
        lift_direction_unit = np.array([0, 0, 1])
        side_direction_unit = np.array([0, 1, 0])
    else:
        lift_direction = (
            vel_wind - np.dot(vel_wind, drag_direction_unit) * drag_direction_unit
        )
        lift_direction_unit = lift_direction / np.linalg.norm(lift_direction)
        side_direction_unit = np.cross(lift_direction_unit, drag_direction_unit)

    for ci, cj in zip(bridle_ci, bridle_cj):
        ## Effective bridle surface area (orthographic projection)
        # getting the bridle unit vector
        node_a = points[ci]
        node_b = points[cj]
        length_bridle = np.linalg.norm(node_b - node_a)
        unit_vector_bridle = (node_b - node_a) / length_bridle
        # getting the bridle velocity vector
        vel_bridle_vec = (
            vel_app  ##TODO: this should be a func. of the V_node_a and V_node_b
        )
        vel_bridle_abs = np.linalg.norm(vel_bridle_vec)
        vel_bridle_unit = vel_bridle_vec / vel_bridle_abs

        # getting the effective surface_area
        area_bridle_eff = (
            bridle_diameter
            * length_bridle
            * (1 - np.dot(unit_vector_bridle, vel_bridle_unit))
        )

        ## Bridle drag force
        D_bridle = 0.5 * RHO * (vel_bridle_abs**2) * bridle_cd * area_bridle_eff

        ## Adding the bridle-drag force to the nodes
        force_aero_bridle[ci] += 0.5 * D_bridle * drag_direction_unit
        force_aero_bridle[cj] += 0.5 * D_bridle * drag_direction_unit

    ##TODO: include an aero-pulley force

    """ --- KCU aero ---
    Aerodynamic Lift force is neglected

    Aerodynamic Drag force (acts as a damping term)
        D = 0.5 * RHO * Vc^2 * Cd * area_kcu
    RHO is the air density [ISA]
    vel_kcu is the velocity of the KCU = Vw (QSM)

    area_kcu is the ref-surface area of KCU = 0.113 m^2
        area_kcu = pi * (kcu_diam_circle/2)^2
            kcu_diam_circle = 0.38 m [Geschiere2019] ##TODO: needs a better source
            calculate by viewing the KCU as isoscele triangle
            which has volume = 0.029m3 = (0.5*width*height*length)
            and then calculating the diameter of a circle with the same area
            kcu_diam_circle = 2*[3*volume/(4*pi)]^(1/3)

    kcu_cd is the drag coefficient = 0.47 [Geschiere2019]
        Comes from empirical SPHERICAL data [Horner 1965]
        Requires a Reynolds number = 2e5-8e5 [Geschiere2019]
            Re = RHO * vel_kcu * kcu_diam_circle / mu
            RHO = 1.225 kg/m^3
            vel_kcu = 30 m/s ##TODO: needs a good source
            mu = 18E-6 kg/m/s [at 15C] ##TODO: needs a good source
            kcu_diam_circle = 0.38 m [Geschiere2019] ##TODO: needs a better source
    """
    vel_kcu = vel_app  ##TODO: this should be calculated for non-QSM
    vel_kcu_abs = np.linalg.norm(vel_kcu)
    vel_kcu_unit = vel_kcu / vel_kcu_abs
    area_kcu = np.pi * ((kcu_diameter / 2) ** 2)
    drag_kcu = 0.5 * RHO * (vel_kcu_abs**2) * kcu_cd * area_kcu

    # adding the KCU-drag force to the nodes
    force_aero_bridle[kcu_index] += drag_kcu * vel_kcu_unit

    return force_aero_bridle


def calculate_force_aero_bridle_thedens2022(points, vel_app, input_bridle_aero):
    """--- Bridle aero ---
    Aerodynamic Lift force is neglected as cylinders only have substantial lifting force
    when rotating, which they don't really do [needs source]

    Aerodynamic Drag force (acts as a damping term)
        D = 0.5 * RHO * Vb^2 * Cd_b * Sb
    RHO is the air density [ISA]
    Vb is the velocity of the bridle = Vw (QSM)
    Sb is the surface area of bridle

    Cd_b is the drag coefficient = 1.1 ##TODO: verify this
        Comes from empirical CYLINDRICAL data [Horner 1965]
        For a Reynolds of approximately ~ 2E4 (up to 8e4, according to Geschiere2019)
            1.225*30*1E-3 / (18E-6) = 2E4
        Bridles ~ 1mm diameter (Geschiere says 2mm?) ##TODO: check this
        Vb      ~ 30m/s ##TODO: needs a good source
        RHO     ~ 1.225 kg/m^3
        mu      ~ 18E-6 kg/m/s [at 15C] ##TODO: needs a good source
    """

    vel_wind = input_bridle_aero.vel_wind
    bridle_ci = input_bridle_aero.bridle_ci
    bridle_cj = input_bridle_aero.bridle_cj
    RHO = input_bridle_aero.rho
    bridle_cd = input_bridle_aero.cd_cylinder
    bridle_cd_shear = input_bridle_aero.cd_shear_cylinder
    bridle_diameter = input_bridle_aero.diameter
    kcu_cd = input_bridle_aero.kcu_cd
    kcu_index = input_bridle_aero.kcu_index
    kcu_diameter = input_bridle_aero.kcu_diameter

    force_aero_bridle = np.zeros(points.shape)

    ## determining directions
    drag_direction_unit = vel_app / np.linalg.norm(vel_app)
    if (
        vel_wind.all() == vel_app.all()
    ):  # avoid division by 0, ##TODO: remove hardcoding
        lift_direction_unit = np.array([0, 0, 1])
        side_direction_unit = np.array([0, 1, 0])
    else:
        lift_direction = (
            vel_wind - np.dot(vel_wind, drag_direction_unit) * drag_direction_unit
        )
        lift_direction_unit = lift_direction / np.linalg.norm(lift_direction)
        side_direction_unit = np.cross(lift_direction_unit, drag_direction_unit)

    for ci, cj in zip(bridle_ci, bridle_cj):
        ## Effective bridle surface area (orthographic projection)
        # getting the bridle unit vector
        node_a = points[ci]
        node_b = points[cj]
        length_bridle = np.linalg.norm(node_b - node_a)
        unit_vector_bridle = (node_b - node_a) / length_bridle
        # getting the bridle velocity vector
        vel_bridle_vec = vel_app  ##TODO: if dynamic this should be a func. of the V_node_a and V_node_b
        vel_bridle_abs = np.linalg.norm(vel_bridle_vec)
        vel_bridle_unit = vel_bridle_vec / vel_bridle_abs

        # angle between the bridle and the velocity vector
        angle_bridle_vel = np.arccos(np.dot(unit_vector_bridle, vel_bridle_unit))

        # getting the forces
        lift_coefficient_bridle = (
            bridle_cd * (np.sin(angle_bridle_vel) ** 2) * (np.cos(angle_bridle_vel))
        )
        drag_coefficient_bridle = (
            bridle_cd * (np.sin(angle_bridle_vel) ** 3) + bridle_cd_shear
        )

        lift_bridle = (
            0.5
            * RHO
            * (vel_bridle_abs**2)
            * lift_coefficient_bridle
            * (bridle_diameter * length_bridle)
        )
        drag_bridle = (
            0.5
            * RHO
            * (vel_bridle_abs**2)
            * drag_coefficient_bridle
            * (bridle_diameter * length_bridle)
        )

        ## Adding the bridle-drag force to the nodes
        force_aero_bridle[ci] += 0.5 * (
            lift_bridle * lift_direction_unit + drag_bridle * drag_direction_unit
        )
        force_aero_bridle[cj] += 0.5 * (
            lift_bridle * lift_direction_unit + drag_bridle * drag_direction_unit
        )

    ##TODO: include an aero-pulley force

    """ --- KCU aero ---
    Aerodynamic Lift force is neglected

    Aerodynamic Drag force (acts as a damping term)
        D = 0.5 * RHO * Vc^2 * Cd * area_kcu
    RHO is the air density [ISA]
    vel_kcu is the velocity of the KCU = Vw (QSM)

    area_kcu is the ref-surface area of KCU = 0.113 m^2
        area_kcu = pi * (drag_kcu/2)^2
            drag_kcu = 0.38 m [Geschiere2019] ##TODO: needs a better source
            calculate by viewing the KCU as isoscele triangle
            which has volume = 0.029m3 = (0.5*width*height*length)
            and then calculating the diameter of a circle with the same area
            drag_kcu = 2*[3*volume/(4*pi)]^(1/3)

    kcu_cd is the drag coefficient = 0.47 [Geschiere2019]
        Comes from empirical SPHERICAL data [Horner 1965]
        Requires a Reynolds number = 2e5-8e5 [Geschiere2019]
            Re = RHO * vel_kcu * drag_kcu / mu
            RHO = 1.225 kg/m^3
            vel_kcu = 30 m/s ##TODO: needs a good source
            mu = 18E-6 kg/m/s [at 15C] ##TODO: needs a good source
            drag_kcu = 0.38 m [Geschiere2019] ##TODO: needs a better source
    """
    vel_kcu = vel_app  ##TODO: this should be calculated for non-QSM
    vel_kcu_abs = np.linalg.norm(vel_kcu)
    vel_kcu_unit = vel_kcu / vel_kcu_abs
    area_kcu = np.pi * ((kcu_diameter / 2) ** 2)
    drag_kcu = 0.5 * RHO * (vel_kcu_abs**2) * kcu_cd * area_kcu

    # adding the KCU-drag force to the nodes
    force_aero_bridle[kcu_index] += drag_kcu * vel_kcu_unit

    return force_aero_bridle
