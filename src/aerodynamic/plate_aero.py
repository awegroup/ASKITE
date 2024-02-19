import numpy as np

# TODO: below is the code used to call this function

# Plate-aero static aero
# AERO_CONFIG = {'RING_GEOMETRY': config.aero.ring_geometry,
#  'MODEL': config.aero.model,
#  'N_ITER': config.aero.n_iter,
#  'ERROR_LIMIT': config.aero.error_limit,
#  'RELAX_FACTOR': config.aero.relax_factor,
#  'N_SPLITS': config.aero.n_splits,
#  'BRIDLE_CD': config.aero.cd_cylinder,
#  'BRIDLE_CD_SHEAR': config.aero.cd_shear_cylinder,
#  'CD_MULTIPLIER': config.aero.cd_multiplier,}

# KCU_DATA = {'CD': config.kite.kcu.drag_coefficient,
#             'DIAMETER': config.kite.kcu.diameter,
#             'INDEX': config.kite.kcu.index,
#             'MASS': config.kite.kcu.mass}

# force_aero_wing = plate_aero.calculate_force_aero_plate(config.kite.connectivity.plate_point_indices,config.kite.points_ini,vel_app_norm,config.kite.area_projected,config.rho,equal_boolean=False)
# force_aero_bridle = bridle_line_system_aero.calculate_force_aero_bridle_thedens2022(points, vel_app, config)
# force_aero = force_aero_wing + force_aero_bridle
# print(f' ')
# print(f'    plate-aero')
# coord_refined, controlpoints, rings, wingpanels, ringvec, coord_L, n_chordwise_elements_refined, config.kite.n_segments_refined = coupling_struc2aero.struc2aero(
# config.kite.points_ini,vel_app,n_chordwise_elements,config.kite.connectivity.plate_point_indices,BILLOWING_ANGLES,AERO_CONFIG)
# functions_print.print_aero(coord_refined,config.vel_wind,vel_app,vel_app_norm,force_aero,force_aero_bridle,
#    force_aero_wing,AERO_CONFIG,config.kite.ref_chord,MU,config.rho,config.kite.area_projected,KCU_DATA['INDEX'],
#    config.kite.points_ini[config.kite.connectivity.plate_point_indices[int(config.kite.n_segments/2)]][:2])


def calculate_force_aero_plate(
    plate_point_indices, pos, Va_norm, A_projected, rho, equal_boolean=False
):
    # Specifying tubular frame + plate diagonals.
    # Reasons for incorporating plate diagonals is that otherwise the quadrilaterals can deform into other shapes
    # This phenomena is best explained by looking at a rectangle that turns into a diamond/kite shape

    aero_force = np.zeros(pos.shape)

    for i in np.arange(0, len(plate_point_indices)):  # looping through each panel
        ### Calculating the area of the plate ##TODO: this could be part of a class
        # calculating the sides
        side_1 = np.linalg.norm(
            pos[plate_point_indices[i][0]] - pos[plate_point_indices[i][1]]
        )
        side_2 = np.linalg.norm(
            pos[plate_point_indices[i][1]] - pos[plate_point_indices[i][2]]
        )
        side_3 = np.linalg.norm(
            pos[plate_point_indices[i][2]] - pos[plate_point_indices[i][3]]
        )
        side_4 = np.linalg.norm(
            pos[plate_point_indices[i][3]] - pos[plate_point_indices[i][0]]
        )

        # Calculating the semi-perimeter (s) of the given quadilateral
        s = (side_1 + side_2 + side_3 + side_4) / 2

        # Applying Brahmagupta's formula to #https://en.wikipedia.org/wiki/Brahmagupta%27s_formula
        # get maximum area of quadrilateral
        area_p = np.sqrt((s - side_1) * (s - side_2) * (s - side_3) * (s - side_4))
        # a,b,c,d = side_1,side_2,side_3,side_4
        # area_p_other = (1/4) * np.sqrt((-a+b+c+d)*(a-b+c+d)*(a+b-c+d)*(a+b+c-d))

        ### Calculating the angle of attack
        middle_LE_point = (
            np.subtract(pos[plate_point_indices[i][1]], pos[plate_point_indices[i][0]])
            / 2
        ) + pos[plate_point_indices[i][0]]
        middle_TE_point = (
            np.subtract(pos[plate_point_indices[i][3]], pos[plate_point_indices[i][2]])
            / 2
        ) + pos[plate_point_indices[i][2]]
        middle_vec = np.subtract(middle_TE_point, middle_LE_point)
        middle_vec_unit = middle_vec / np.linalg.norm(middle_vec)
        vw_vec_unit = np.array([1, 0, 0]) / np.linalg.norm(np.array([1, 0, 0]))
        aoa_p = np.arccos(
            np.dot(middle_vec_unit, vw_vec_unit)
            / (np.linalg.norm(middle_vec_unit) * np.linalg.norm(vw_vec_unit))
        )

        ###  Define the perpendicular unit vector to each segment
        # Defining 2 vectors tangential to the plane, by using the diagonals
        diagonal_1 = np.subtract(
            pos[plate_point_indices[i][0]], pos[plate_point_indices[i][2]]
        )
        diagonal_2 = np.subtract(
            pos[plate_point_indices[i][1]], pos[plate_point_indices[i][3]]
        )

        # Finding the cross product of these two vectors, should give a perpendicular vector / two solutions actually.
        perpendicular_vector_1 = np.cross(diagonal_1, diagonal_2)
        unit_perp_vec = perpendicular_vector_1 / np.linalg.norm(perpendicular_vector_1)

        ### Find the direction of lift and drag
        L_vec_unit = unit_perp_vec

        # correcting the angle of attack for the orientation of the plate
        if unit_perp_vec[0] > 0:  # the plate is tilted backwards
            aoa_p = aoa_p
        elif unit_perp_vec[0] < 0:  # the plate is tilted forwards
            aoa_p = -aoa_p

        # AoA is equal to aoa_p
        Cl_p = 2 * np.pi * np.sin(aoa_p)
        L_p = 0.5 * rho * (Va_norm**2) * Cl_p * (area_p)

        ## Splitting lift
        Fx_p = L_p * L_vec_unit[0]
        Fy_p = L_p * L_vec_unit[1]  # side-force
        Fz_p = L_p * L_vec_unit[2]

        ## Apply 1/4 of perpendicular unit vector to each node of the respective panel
        for k, j in enumerate(
            plate_point_indices[i]
        ):  # loop Through all the indexes for the defined plate
            if equal_boolean == True:
                aero_force[j, 0] += 0.25 * Fx_p  # Drag only works parallel to the wind
                aero_force[j, 1] += 0.25 * Fy_p
                aero_force[j, 2] += 0.25 * Fz_p
            else:
                if k == 0 or k == 1:  # LEADING EDGE
                    aero_force[j, 0] += (
                        0.25 * Fx_p
                    )  # Drag only works parallel to the wind
                    aero_force[j, 1] += 0.5 * 0.75 * Fy_p
                    aero_force[j, 2] += 0.5 * 0.75 * Fz_p
                elif k == 2 or k == 3:  # TRAILING EDGE
                    aero_force[j, 0] += (
                        0.25 * Fx_p
                    )  # Drag only works parallel to the wind
                    aero_force[j, 1] += 0.5 * 0.25 * Fy_p
                    aero_force[j, 2] += 0.5 * 0.25 * Fz_p

    return np.array(aero_force)
