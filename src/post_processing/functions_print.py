# importing libraries
import numpy as np

# own functions
from src.post_processing import post_processing_utils


def calculate_slack_values(
    pos, points_ini, ci, cj, ci_kite, cj_kite, springL, springL_kite, TE_extension, U_P
):
    slack_index = []

    for i in range(len(ci)):
        if i > 21:
            # LE bridles
            slack_index.append("#00000080")
        else:
            # TE bridles
            slack_index.append("#00000080")

    dL_kite_stretch_lst, dL_kite_slack_lst = np.zeros(len(ci_kite)), np.zeros(
        len(ci_kite)
    )
    for i in range(len(ci_kite)):
        sep_vec = pos[ci_kite[i]] - pos[cj_kite[i]]
        sep = np.linalg.norm(sep_vec)
        dL_perc = 100 * ((sep - springL_kite[i]) / sep)

        if dL_perc > 0:
            dL_kite_stretch_lst[i] = dL_perc
        else:
            dL_kite_slack_lst[i] = dL_perc

    dL_bridle_slack_lst, dL_bridle_stretch_lst = np.zeros(len(ci)), np.zeros(len(ci))
    for i in range(len(ci)):
        if i == 4 or i == 5:
            sep_vec_1 = pos[ci[i]] - pos[cj[i]]
            sep_1 = np.linalg.norm(sep_vec_1)
            unit_vector_1 = sep_vec_1 / sep_1

            i_n = i * 2 - 2
            sep_vec_2 = pos[ci[i_n]] - pos[cj[i_n]]
            sep_2 = np.linalg.norm(sep_vec_2)
            unit_vector_2 = sep_vec_2 / sep_2

            dL = ((sep_1 + sep_2) - (springL[i] + springL[i_n])) / (
                springL[i] + springL[i_n]
            )
            dL_perc = 100 * dL

        elif i != 6 and i != 8:
            sep_vec = pos[ci[i]] - pos[cj[i]]
            sep = np.linalg.norm(sep_vec)
            dL_perc = 100 * ((sep - springL[i]) / springL[i])
        elif i == 6:
            dL_perc = dL_bridle_stretch_lst[3]
        elif i == 8:
            dL_perc = dL_bridle_stretch_lst[5]

        if dL_perc > 0:
            dL_bridle_stretch_lst[i] = dL_perc
        else:
            dL_bridle_slack_lst[i] = dL_perc
            if -2.5 < dL_perc <= 0:
                slack_index[i] = "#1CEFCC"
            elif -5 < dL_perc <= -2.5:
                slack_index[i] = "#FFBD19"
            elif dL_perc <= -5:
                slack_index[i] = "#EF1CEC"

    return (
        pos,
        slack_index,
        dL_bridle_stretch_lst,
        dL_bridle_slack_lst,
        dL_kite_stretch_lst,
        dL_kite_slack_lst,
    )


def print_settings(vel_app, config):
    ## initialising
    U_P = config.u_p
    N_ITER = config.aero.n_iter
    ERROR_LIMIT = config.aero.error_limit
    RELAX_FACTOR = config.aero.relax_factor
    N_SPLITS = config.aero.n_splits
    AEROSTRUCTURAL_TOL = 1e-3 * config.aero_structural.tol
    STRUCTURAL_TOL = config.solver.tol
    STRUCTURAL_INTEGRATION_METHOD = config.solver.integration_method
    points_mass = config.kite.mass_points

    print(f" ")
    print(f"------------------------------------------------------------")
    print(f"    SETTINGS")
    print(f"------------------------------------------------------------")
    print(
        f"GENERAL --> U_P: {U_P}, vel_app_norm: {np.linalg.norm(vel_app):.2f} m/s, AeroStructural tol: {AEROSTRUCTURAL_TOL:.5f} mm, kite mass: {np.sum(points_mass):.2f} kg"
    )
    print(
        f"AERO    --> N_split: {N_SPLITS}, N_iter: {N_ITER}, error: {ERROR_LIMIT}, Relax_factor: {RELAX_FACTOR} "
    )
    print(
        f"STRUC   --> Tol: {str(STRUCTURAL_TOL)}, Integration Method: {STRUCTURAL_INTEGRATION_METHOD}"
    )


def print_initial_kite_dimensions(config):
    print(f" ")
    print(f"------------------------------------------------------------")
    print(f"    INITIAL KITE DIMENSIONS")
    print(f"------------------------------------------------------------")
    print(f"scaling-factor: {config.geometric_scaling_factor}")
    print(f"ref_chord: {config.kite.ref_chord:.2f}m")
    print(f"wing_span: {config.kite.span:.2f}m")
    print(f"wing height: {config.kite.height:.2f}m")
    print(f"wing area: {config.kite.area_surface:.2f}m2")
    print(f"projected_area: {config.kite.area_projected:.2f}m")
    print(
        f"At 25m/s --> Reynolds number: {1e-6*(1.225*25*config.kite.ref_chord)/(1.8e-5):.2f}e6"
    )


def calculate_trim_angle(coord_refined, N_split, points_LE_segment):
    """Returns the trim angle of the kite in degrees (hardcoded)"""

    LE_point_diff = points_LE_segment[1] - points_LE_segment[0]
    mid_point_LE = points_LE_segment[0] + 0.5 * LE_point_diff

    # # aerodynamic trim-angle, gonna change it to bridle defined trim_angle
    # mid_idx = int(len(coord_refined)/2)
    # if N_split%2 == 0: #if N_split is even
    #     mid_point_LE = coord_refined[mid_idx-1]
    # else: #if N_split is odd
    #     mid_point_LE_left = coord_refined[mid_idx-2]
    #     mid_point_LE_right = coord_refined[mid_idx]
    #     vec_left_to_right = mid_point_LE_right - mid_point_LE_left
    #     mid_point_LE = mid_point_LE_left + 0.5*vec_left_to_right

    trim_angle = np.arctan(mid_point_LE[0] / mid_point_LE[2])
    return np.rad2deg(trim_angle)


def calculate_side_angle_leading_edge_midpoint(points_LE_segment):
    LE_point_diff = points_LE_segment[1] - points_LE_segment[0]
    mid_point_LE = points_LE_segment[0] + 0.5 * LE_point_diff
    side_angle = np.arctan(mid_point_LE[1] / mid_point_LE[2])
    return np.rad2deg(side_angle)


def print_aerostructural(boolean_convergence, it, max_iter, total_time):
    ## Printing solution
    print(f" ")
    if boolean_convergence:
        print(f"------------------------------------------------------------")
        print(f"    CONVERGENCE \-,-/")
    else:
        print(f"------------------------------------------------------------")
        print(f"    NO CONVERGENCE, reached max iterations")

    print(f"------------------------------------------------------------")
    print(f"Global iterations   --> Needed {it} | MAX  {max_iter}")
    print(f"Time needed         --> {np.round(total_time,2)} s")


def print_structural(
    points,
    bridle_rest_lengths,
    wing_rest_lengths,
    config,
):
    ## initialising
    points_ini = config.kite.points_ini
    U_P = config.u_p
    bridle_stiffness = config.kite.stiffness.bridle

    stretch_values, slack_values = post_processing_utils.calculate_elongation(
        points,
        wing_rest_lengths,
        bridle_rest_lengths,
        config,
    )[:2]
    (
        dL_tube_stretch_lst,
        dL_canopy_stretch_lst,
        dL_bridle_stretch_lst,
        dL_wing_stretch_lst,
    ) = stretch_values
    (
        dL_tube_slack_lst,
        dL_canopy_slack_lst,
        dL_bridle_slack_lst,
        dL_wing_slack_lst,
    ) = slack_values

    ## Ensuring errors from empty arrays are dealt with
    if dL_tube_stretch_lst == []:
        dL_tube_stretch_lst = [0]
    if dL_canopy_stretch_lst == []:
        dL_canopy_stretch_lst = [0]
    if dL_bridle_stretch_lst == []:
        dL_bridle_stretch_lst = [0]
    if dL_wing_stretch_lst == []:
        dL_wing_stretch_lst = [0]
    if dL_tube_slack_lst == []:
        dL_tube_slack_lst = [0]
    if dL_canopy_slack_lst == []:
        dL_canopy_slack_lst = [0]
    if dL_bridle_slack_lst == []:
        dL_bridle_slack_lst = [0]
    if dL_wing_slack_lst == []:
        dL_wing_slack_lst = [0]

    print(f"")
    print(f"------------------------------------------------------------")
    print(
        f"  STRUCTURAL DEFORMATIONS in mm (BRIDLE STIFFNESS: {bridle_stiffness/1e3:.2f} kN/m)"
    )
    print(f"------------------------------------------------------------")
    print(
        f"Depower tape length  --> Diff {np.round(np.linalg.norm(points[22]-points[21]) - bridle_rest_lengths[1],2)}   | New: {np.round(np.linalg.norm(points[22]-points[21]),2)} | Initial: {np.round(bridle_rest_lengths[1],2)}, u_p: {U_P}"
    )
    print(
        f"Width TE (tip-to-tip)--> Diff {np.round(np.linalg.norm(points[1]-points[10]) - np.linalg.norm(points_ini[1]-points_ini[10]),2)} | New: {np.round(np.linalg.norm(points[1]-points[10]),2)} | Initial: {np.round(np.linalg.norm(points_ini[1]-points_ini[10]),2)} "
    )
    print(
        f"Width LE (tip-to-tip)--> Diff {np.round(np.linalg.norm(points[19]-points[20]) - np.linalg.norm(points_ini[19]-points_ini[20]),2)} | New: {np.round(np.linalg.norm(points[19]-points[20]),2)} | Initial: {np.round(np.linalg.norm(points_ini[19]-points_ini[20]),2)}"
    )
    print(
        f"Height Kite          --> Diff {np.round(np.max(points[:,2]) - np.max(points_ini[:,2]),2)}"
    )
    print(f"")
    print(
        f"SLACK Tubular frame  --> AVE {np.round(np.average(dL_tube_slack_lst),2)}% | MAX {np.round(np.min(dL_tube_slack_lst),2)}%  [should be close to 0]"
    )
    print(
        f"SLACK Canopy         --> AVE {np.round(np.average(dL_canopy_slack_lst),2)}% | MAX {np.round(np.min(dL_canopy_slack_lst),2)}% "
    )
    print(
        f"SLACK Bridle         --> AVE {np.round(np.average(dL_bridle_slack_lst),2)}% | MAX {np.round(np.min(dL_bridle_slack_lst),2)}% "
    )
    print(
        f"STRETCH Tubular frame--> AVE {np.round(np.average(dL_tube_stretch_lst),2)}% | MAX {np.round(np.max(dL_tube_stretch_lst),2)}%  [should be close to 0]"
    )
    print(
        f"STRETCH Canopy       --> AVE {np.round(np.average(dL_canopy_stretch_lst),2)}% | MAX {np.round(np.max(dL_canopy_stretch_lst),2)}%  [should be close to 0]"
    )
    print(
        f"STRETCH Bridle       --> AVE  {np.round(np.average(dL_bridle_stretch_lst),2)}% | MAX  {np.round(np.max(dL_bridle_stretch_lst),2)}%  [should be close to 0]"
    )


def calculate_lift_from_force_aero(force_aero, vel_app, VEL_WIND):
    in_plane_perpendicular_forces = []
    plane_normal = np.cross(vel_app, VEL_WIND)
    for force in force_aero:
        in_plane_perpendicular_force = np.dot(force, plane_normal) / np.linalg.norm(
            plane_normal
        )
        in_plane_perpendicular_forces.append(in_plane_perpendicular_force)
    return np.sum(in_plane_perpendicular_forces)


def calculate_drag_from_force_aero(force_aero, vel_app, VEL_WIND):
    parallel_forces = []
    for force in force_aero:
        parallel_force = np.dot(force, vel_app) / np.linalg.norm(vel_app)
        parallel_forces.append(parallel_force)
    return np.sum(parallel_forces)


def calculate_sideforce_from_force_aero(force_aero, vel_app, VEL_WIND):
    out_of_plane_perpendicular_forces = []
    for force in force_aero:
        out_of_plane_perpendicular_force = np.dot(
            force, np.cross(vel_app, VEL_WIND)
        ) / np.linalg.norm(np.cross(vel_app, VEL_WIND))
        out_of_plane_perpendicular_forces.append(out_of_plane_perpendicular_force)
    return np.sum(out_of_plane_perpendicular_forces)


# def print_aero_old(coord_refined,VEL_WIND,vel_app,vel_app_norm,force_aero,force_aero_bridle,force_aero_wing,AERO_CONFIG,REF_CHORD,MU,RHO,AREA_PROJECTED,KCU_INDEX,points_LE_segment):

#     N_SPLITS = AERO_CONFIG['N_SPLITS']
#     MODEL = AERO_CONFIG['MODEL']
#     RING_GEOMETRY = AERO_CONFIG['RING_GEOMETRY']

#     print(f'')
#     print(f'{MODEL}-AERODYNAMICS (|Va|: {vel_app_norm:.2f} m/s, Va: {vel_app}, Rey: {1e-6*((RHO*vel_app_norm*(REF_CHORD))/MU):.2f}e6, ring_geo: {RING_GEOMETRY})') ##TODO: max_chord is now hardcoded should be something like np.max(springL_wing)
#     print(f'------------------------------------------')
#     #print(f'Angle of attack(deg) --> Diff {np.round(np.rad2deg(aoa-aoa0),2)} | NEW: {np.round(np.rad2deg(aoa),2)} | Initial: {np.round(np.rad2deg(aoa0),2)}')

#     trim_angle = calculate_trim_angle(coord_refined,N_SPLITS,points_LE_segment)
#     print(f'Aero Force           --> Fx {np.sum(force_aero[:,0]):.2f} N, Fy {np.sum(force_aero[:,1]):.2f} N, Fz {np.sum(force_aero[:,2]):.2f} N')
#     print(f'Trim-angle           --> {trim_angle:.2f} deg (= from vertical to front-bridle)')

#     ## determining directions
#     drag_direction_unit = vel_app / np.linalg.norm(vel_app)
#     if VEL_WIND.all() == vel_app.all(): #avoid division by 0, ##TODO: remove hardcoding
#         lift_direction_unit = np.array([0,0,1])
#         side_direction_unit = np.array([0,1,0])
#     else:
#         lift_direction = VEL_WIND - np.dot(VEL_WIND, drag_direction_unit) * drag_direction_unit
#         lift_direction_unit = lift_direction / np.linalg.norm(lift_direction)
#         side_direction_unit = np.cross(lift_direction_unit,drag_direction_unit)

#     ## determining forces
#     force_lift = np.sum(force_aero*lift_direction_unit)
#     force_drag = np.sum(force_aero*drag_direction_unit)
#     force_side = np.sum(force_aero*side_direction_unit)
#     ## determining coefficients
#     C_L = np.sum(force_lift)/ (0.5*RHO*(vel_app_norm**2)*AREA_PROJECTED)
#     C_D = np.sum(force_drag)/ (0.5*RHO*(vel_app_norm**2)*AREA_PROJECTED)
#     C_S = np.sum(force_side)/ (0.5*RHO*(vel_app_norm**2)*AREA_PROJECTED)

#     print(f'Lift                 --> {force_lift:.2f} N, CL: {C_L:.3f}, CL/CD: {C_L/C_D:.2f} ("glide-ratio E")')

#     force_drag_bridle = calculate_drag_from_force_aero(force_aero_bridle, vel_app, VEL_WIND)
#     force_drag_kcu = calculate_drag_from_force_aero(force_aero[KCU_INDEX], vel_app, VEL_WIND)
#     force_drag_bridle_minus_kcu = force_drag_kcu - force_drag_bridle
#     force_drag_wing = calculate_drag_from_force_aero(force_aero_wing, vel_app, VEL_WIND)

#     print(f'Drag                 --> {force_drag:.2f} N, CD: {C_D:.3f} (Bridle: {force_drag_bridle_minus_kcu:.2f} N, KCU: {force_drag_kcu:.2f} N, Wing: {force_drag_wing:.2f} N)')

#     # ## computing one half's side force
#     # F_side_half = []
#     # wing_right_side_indices = [6,7,8,9,20,1,0,11,12,13,14]
#     # for i,force in enumerate(force_aero):
#     #     if i in wing_right_side_indices: ##TODO: remove hardcode
#     #         F_side_half.append(force[1])
#     print(f'Side-Force           --> {force_side:.2f} N, CS: {C_S:.3f}')  #(Half-wing: {np.max(F_side_half):.2f} N)')

#     ##TODO: change force_aero --> F_aero

#     return


def print_forces(
    f_internal,
    f_external,
    force_aero,
    residual_f_including_fixed_nodes,
    residual_f,
    f_tether_drag,
    f_gravity,
    points,
):
    print(f" ")
    print(f"------------------------------------------------------------")
    print(f"    FORCES")
    print(f"------------------------------------------------------------")
    print(
        f"F_aero                          : Fx: {np.sum(force_aero[:,0]):.2f}N, Fy: {np.sum(force_aero[:,1]):.2f}N, Fz: {np.sum(force_aero[:,2]):.2f}N"
    )
    print(
        f"F_gravity                       : Fx: {np.sum(f_gravity[:,0]):.2f}N, Fy: {np.sum(f_gravity[:,1]):.2f}N, Fz: {np.sum(f_gravity[:,2]):.2f}N"
    )
    print(f"_________+")
    f_external_reshaped = np.reshape(f_external, (len(points), 3))
    print(
        f"F_external                      : Fx: {np.sum(f_external_reshaped[:,0]):.2f}N, Fy: {np.sum(f_external_reshaped[:,1]):.2f}N, Fz: {np.sum(f_external_reshaped[:,2]):.2f}N"
    )
    f_internal_reshaped = np.reshape(f_internal, (len(points), 3))
    print(
        f"F_internal                      : Fx: {np.sum(f_internal_reshaped[:,0]):.2f}N, Fy: {np.sum(f_internal_reshaped[:,1]):.2f}N, Fz: {np.sum(f_internal_reshaped[:,2]):.2f}N"
    )
    print(f"__________+")
    residual_f_incl = np.reshape(residual_f_including_fixed_nodes, (len(points), 3))
    print(
        f"F_residual (incl. fixed nodes)  : Fx: {np.sum(residual_f_incl[:,0]):.2f}N, Fy: {np.sum(residual_f_incl[:,1]):.2f}N, Fz: {np.sum(residual_f_incl[:,2]):.2f}N"
    )
    # TODO: this now assumes that there is only 1 fixed node
    residual_f_excl = np.reshape(residual_f, (len(points) - 1, 3))
    print(
        f"F_residual (excl. fixed nodes)  : Fx: {np.sum(residual_f_excl[:,0]):.2f}N, Fy: {np.sum(residual_f_excl[:,1]):.2f}N, Fz: {np.sum(residual_f_excl[:,2]):.2f}N"
    )
    print(f"__________-")
    f_bridle_point = [
        np.sum(residual_f_incl[:, 0]) - np.sum(residual_f_excl[:, 0]),
        np.sum(residual_f_incl[:, 1]) - np.sum(residual_f_excl[:, 1]),
        np.sum(residual_f_incl[:, 2]) - np.sum(residual_f_excl[:, 2]),
    ]
    print(
        f"F_bridle_point (fixed node)     : Fx: {f_bridle_point[0]:.2f}N, Fy: {f_bridle_point[1]:.2f}N, Fz: {f_bridle_point[2]:.2f}N"
    )
    print(
        f"F_tether_drag                   : Fx: {f_tether_drag[0]:.2f}N, Fy: {f_tether_drag[1]:.2f}N, Fz: {f_tether_drag[2]:.2f}N"
    )
    print(f"__________+")
    f_tether_excl = f_bridle_point + f_tether_drag
    print(
        f"F_tether (due to kite forces)   : Fx: {f_tether_excl[0]:.2f}N, Fy: {f_tether_excl[1]:.2f}N, Fz: {f_tether_excl[2]:.2f}N"
    )


def print_aero(points, vel_app, force_aero_wing, force_aero_bridle, config):
    ## initialisation
    REF_CHORD = config.kite.ref_chord

    KCU_INDEX = config.kite.kcu.index
    AREA_PROJECTED = config.kite.area_projected
    plate_point_indices = config.kite.connectivity.plate_point_indices
    n_segments = config.kite.n_segments
    VEL_WIND = config.vel_wind
    N_SPLITS = config.aero.n_splits
    MODEL = config.aero.model
    RING_GEOMETRY = config.aero.ring_geometry

    RHO = config.rho
    MU = config.mu

    points_LE_segment = points[plate_point_indices[int(n_segments / 2)]][:2]
    points_TE_segment = points[plate_point_indices[int(n_segments / 2)]][2:]
    LE_point_diff = points_LE_segment[1] - points_LE_segment[0]
    TE_point_diff = points_TE_segment[1] - points_TE_segment[0]
    mid_point_LE = points_LE_segment[0] + 0.5 * LE_point_diff
    mid_point_TE = points_TE_segment[0] + 0.5 * TE_point_diff
    trim_angle = np.arctan(mid_point_LE[0] / mid_point_LE[2])
    delta_z_chord = mid_point_LE[2] - mid_point_TE[2]
    delta_x_chord = np.abs(mid_point_LE[0] - mid_point_TE[0])
    angle_of_mid_chord = np.rad2deg(np.arctan(delta_z_chord / delta_x_chord))
    angle_va_with_xy = np.rad2deg(np.arctan(vel_app[2] / vel_app[0]))
    AoA = angle_of_mid_chord + angle_va_with_xy
    ## determining directions
    drag_direction_unit = vel_app / np.linalg.norm(vel_app)
    side_direction_unit = np.array([0, 1, 0])
    lift_direction_unit = np.cross(drag_direction_unit, side_direction_unit)

    ## determining forces
    force_aero = force_aero_wing + force_aero_bridle
    force_lift = np.sum(force_aero * lift_direction_unit)
    force_drag = np.sum(force_aero * drag_direction_unit)
    force_side = np.sum(force_aero * side_direction_unit)
    ## determining coefficients
    C_L = np.sum(force_lift) / (
        0.5 * 1.225 * (np.linalg.norm(vel_app) ** 2) * AREA_PROJECTED
    )
    C_D = np.sum(force_drag) / (
        0.5 * 1.225 * (np.linalg.norm(vel_app) ** 2) * AREA_PROJECTED
    )
    C_S = np.sum(force_side) / (
        0.5 * 1.225 * (np.linalg.norm(vel_app) ** 2) * AREA_PROJECTED
    )

    # force_drag_kcu = calculate_drag_from_force_aero(
    #     force_aero_wing[KCU_INDEX], vel_app, VEL_WIND
    # )
    # force_drag_bridle_minus_kcu = force_drag_kcu - force_drag_bridle
    force_lift_wing = np.sum(force_aero_wing * lift_direction_unit)
    force_drag_wing = np.sum(force_aero_wing * drag_direction_unit)
    force_side_wing = np.sum(force_aero_wing * side_direction_unit)

    force_lift_bridle = np.sum(force_aero_bridle * lift_direction_unit)
    force_drag_bridle = np.sum(force_aero_bridle * drag_direction_unit)
    force_side_bridle = np.sum(force_aero_bridle * side_direction_unit)

    # force_drag_wing = calculate_drag_from_force_aero(force_aero_wing, vel_app, VEL_WIND)
    # force_drag_bridle = calculate_drag_from_force_aero(
    #     force_aero_bridle, vel_app, VEL_WIND
    # )
    force_lift_kcu = force_lift - (force_lift_bridle + force_lift_wing)
    force_drag_kcu = force_drag - (force_drag_bridle + force_drag_wing)
    force_side_kcu = force_side - (force_side_bridle + force_side_wing)

    ## CALCULATING MOMENTS
    # Full Wing Moments
    points_wing = np.copy(points)
    evaluation_point = np.array([max(points_wing[:, 0]) / 2, 0, max(points_wing[:, 2])])
    evaluation_point = np.zeros(3)
    for i, point in enumerate(points_wing):
        points_wing[i] = point - evaluation_point
        if i not in config.kite.connectivity.plate_point_indices:
            points_wing[i] = np.array([0, 0, 0])
            force_aero[i] = np.array([0, 0, 0])

    moments_wing = np.cross(points_wing, force_aero)

    # Half Wing Moments
    for i, point in enumerate(points_wing):
        if points_wing[i][1] > 0:  # taking only the right wing
            points_wing[i] = np.array([0, 0, 0])
            force_aero[i] = np.array([0, 0, 0])

    moments_half_wing = np.cross(points_wing, force_aero)

    ##TODO: max_chord is now hardcoded should be something like np.max(springL_wing)
    print(f" ")
    print(f"------------------------------------------------------------")
    print(
        f"  AERODYNAMICS -{MODEL}- ( Va: {vel_app} |{np.linalg.norm(vel_app):.2f}|, N_splits:{N_SPLITS}, Rey: {1e-6*(RHO*np.linalg.norm(vel_app*(REF_CHORD))/MU):.2f}e6, ring_geo: {RING_GEOMETRY})"
    )
    print(f"------------------------------------------------------------")
    # print(
    #     f"Aero Force   --> Fx {np.sum(force_aero_wing[:,0]):.2f} N, Fy {np.sum(force_aero_wing[:,1]):.2f} N, Fz {np.sum(force_aero_wing[:,2]):.2f} N"
    # )
    print(
        f"Angles       --> AoA: {AoA:.2f}deg , trim-angle: {np.rad2deg(trim_angle):.2f}deg (defined from z-axis to LE attachment points)"
    )
    print(
        f"Ratios       --> CL: {C_L:.3f}, CD: {C_D:.3f}, CS: {C_S:.3f}, CL/CD: {C_L/C_D:.2f} ('glide-ratio E') "
    )
    print(
        f"Lift         --> {force_lift:.2f} N (Wing: {force_lift_wing:.2f}N, Bridle: {force_lift_bridle:.2f}N, KCU: {force_lift_kcu:.2f}N)"
    )
    print(
        f"Drag         --> {force_drag:.2f} N (Wing: {force_drag_wing:.2f}N, Bridle: {force_drag_bridle:.2f}N, KCU: {force_drag_kcu:.2f}N)"
    )
    print(
        f"Side-Force   --> {force_side:.2f} N (Wing: {force_side_wing:.2f}N, Bridle: {force_side_bridle:.2f}N, KCU: {force_side_kcu:.2f}N)"
    )  # (Half-wing: {np.max(F_side_half):.2f} N)')
    print(
        f"Moments Full Wing on mid-span, mid-chord point --> Mx: {moments_wing[:,0].sum():.3f} My: {moments_wing[:,1].sum():.3f} Mz: {moments_wing[:,2].sum():.2f}"
    )
    print(
        f"Moments Half Right Wing on mid-span, mid-chord point --> Mx: {moments_half_wing[:,0].sum():.3f} My: {moments_half_wing[:,1].sum():.3f} Mz: {moments_half_wing[:,2].sum():.2f}"
    )
    return


def print_kite_matplotlib_old_uri():
    # Plotting solution, geometry

    plt.rcParams.update({"font.size": 10})
    elev = 0
    azim = 180

    fig = plt.figure(figsize=(10, 10))
    # set up the axes for the first plot
    ax = fig.add_subplot(1, 1, 1, projection="3d")
    # functions_plot.plot_kite_matplotlib(points/1000,bridle_ci,bridle_cj,wing_ci,wing_cj,plate_point_indices,ax,True,'black',elev)
    functions_plot.plot_kite_matplotlib_pretty(
        points / 1000,
        bridle_ci,
        bridle_cj,
        ax,
        True,
        "black",
        elev,
        tube_line_indices,
        wing_ci,
        wing_cj,
        plate_point_indices,
    )
    # functions_plot.plot_kite_matplotlib_pretty(points_ini/1000,bridle_ci,bridle_cj,ax,True,'red',elev,tube_line_indices,wing_ci,wing_cj,plate_point_indices)
    ax.grid(False)
    ax.axis("off")
    ax.view_init(elev=elev, azim=azim)
    bbox = fig.bbox_inches.from_bounds(2, 2, 8, 6)

    # plt.savefig('./plots/functions_plot.pdf',bbox_inches = bbox)

    # BOTH KITES
    plt.rcParams.update({"font.size": 10})
    elev = 0
    azim = 90

    fig = plt.figure(figsize=(10, 10))
    # set up the axes for the first plot
    ax = fig.add_subplot(1, 1, 1, projection="3d")
    # functions_plot.plot_kite_matplotlib(points/1000,bridle_ci,bridle_cj,wing_ci,wing_cj,plate_point_indices,ax,True,'black',elev)
    functions_plot.plot_kite_matplotlib_pretty(
        points / 1000,
        bridle_ci,
        bridle_cj,
        ax,
        True,
        "black",
        elev,
        tube_line_indices,
        wing_ci,
        wing_cj,
        plate_point_indices,
    )
    # functions_plot.plot_kite_matplotlib_pretty(points_ini/1000,bridle_ci,bridle_cj,ax,True,'red',elev,tube_line_indices,wing_ci,wing_cj,plate_point_indices)
    functions_plot.plot_kite_matplotlib_pretty(
        points_struc / 1000,
        bridle_ci,
        bridle_cj,
        ax,
        True,
        "green",
        elev,
        tube_line_indices,
        wing_ci,
        wing_cj,
        plate_point_indices,
    )
    ax.grid(False)
    ax.axis("off")
    ax.view_init(elev=elev, azim=azim)
    bbox = fig.bbox_inches.from_bounds(2, 2, 8, 6)

    # plt.savefig('./plots/functions_plot.pdf',bbox_inches = bbox)

    ##%% OLD plots
    # fig = plt.figure(figsize= (10,8))
    ## set up the axes for the first plot
    # elev =  0
    # azim = 180
    # ax = fig.add_subplot(1, 2, 1, projection='3d')
    # functions_plot.plot_kite_matplotlib(points/1000,bridle_ci,bridle_cj,wing_ci,wing_cj,plate_point_indices,slack_index,'black',ax,True,elev)
    ## functions_plot.plot_kite_matplotlib(points_ini/1000,bridle_ci,bridle_cj,wing_ci,plate_point_indices,wing_cj,slack_index,'#FF000080',ax,False,elev)
    # ax.grid(False)
    # ax.axis('off')
    # ax.legend().set_visible(False)
    # ax.view_init(elev = elev, azim = azim)
    # plt.show
    # bbox = fig.bbox_inches.from_bounds(3.3, 3, 5.6, 5)
    #
    ##%% AERO SPECIFIC PLOTS
    # elev = 10
    # azim = 165
    # fig = plt.figure(figsize= (10,10))
    ## set up the axes for the first plot
    # ax = fig.add_subplot(1, 1, 1, projection='3d')
    # VSM.plot_geometry(wingpanels,controlpoints,rings,F_rel,coord_L,ax,'True')
    # ax.grid(False)
    # ax.axis('off')
    # ax.set_xlim((-4,4))
    # ax.set_ylim((-4,4))
    ## ax.set_zlim((0,9))
    # ax.view_init(elev = elev, azim = azim)
    # bbox = fig.bbox_inches.from_bounds(1.5, 3, 7, 4.5)
    ## plt.savefig('./plots/VSM.pdf',bbox_inches = bbox)
    #
    ##%% AERO SPECIFIC PLOTS
    #
    # elev = 10
    # azim = 165
    # fig = plt.figure(figsize= (10,10))
    ## set up the axes for the first plot
    # ax = fig.add_subplot(1, 1, 1, projection='3d')
    # VSM.plot_geometry(wingpanels,controlpoints,rings,F_rel,coord_L,ax,'True')
    # ax.grid(False)
    # ax.axis('off')
    # ax.set_xlim((-4,4))
    # ax.set_ylim((-4,4))
    ## ax.set_zlim((0,9))
    # ax.view_init(elev = elev, azim = azim)
    # bbox = fig.bbox_inches.from_bounds(1.5, 3, 7, 4.5)
    ## plt.savefig('./plots/VSM.pdf',bbox_inches = bbox)
    #
    ##%% AERO SPECIFIC PLOTS
    #
    # plt.rcParams.update({'font.size': 10})
    # elev = 0
    # azim = 180
    #
    # fig = plt.figure(figsize= (10,10))
    ## set up the axes for the first plot
    # ax = fig.add_subplot(1, 1, 1, projection='3d')
    ## functions_plot.plot_kite_matplotlib(points/1000,bridle_ci,bridle_cj,wing_ci,wing_cj,plate_point_indices,ax,True,'black',elev)
    # functions_plot.plot_kite_matplotlib_pretty(points/1000,bridle_ci,bridle_cj,ax,True,'black',elev,tube_line_indices,wing_ci,wing_cj,plate_point_indices)
    ##functions_plot.plot_kite_matplotlib_pretty(points_ini/1000,bridle_ci,bridle_cj,ax,True,'red',elev,tube_line_indices,wing_ci,wing_cj,plate_point_indices)
    # ax.grid(False)
    # ax.axis('off')
    # ax.view_init(elev = elev, azim = azim)
    # bbox = fig.bbox_inches.from_bounds(2, 2, 8, 6)
    #
    ## plt.savefig('./plots/functions_plot.pdf',bbox_inches = bbox)
    #
    #
    ##%% AERO SPECIFIC PLOTS
    # ax.set_title(r'$U_P = $ '+str(U_P)+r' $u_s = $ '+ str(u_s))
    # fig = plt.figure(figsize=plt.figaspect(0.5))
    # force_aero = np.zeros(points.shape) #Defining a matrix, filled with zeros of similar shape as the points matrix
    # force_aero = functions_plot.get_force_aeros2nodes(points,bridle_ci,bridle_cj,plate_point_indices,F_rel, F_mag[:,2],ringvec,controlpoints,force_aero)  # Put in the new pointsitions of the points
    # elev = 0
    # azim = 180
    # ax = fig.add_subplot(1, 2, 1, projection='3d')
    # functions_plot.plot_kite_matplotlibforces(points,bridle_ci,bridle_cj,wing_ci,wing_cj,plate_point_indices,force_aero,ax,elev)
    ## functions_plot.plot_kite_matplotlib(points,bridle_ci,bridle_cj,wing_ci,wing_cj,plate_point_indices,ax,True,'black',elev)
    # ax.grid(False)
    # ax.axis('off')
    # ax.view_init(elev = elev, azim = azim)
    #

    return
