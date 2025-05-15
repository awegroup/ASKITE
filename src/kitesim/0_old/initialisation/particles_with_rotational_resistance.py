import numpy as np
from kitesim.structural import structural_model


def extract_points_between_dict(config):
    # acquiring the attachment points that lie in between
    points_attachment_inbetween = []
    points_index_attachment_inbetween = []
    for idx in config.kite.connectivity.wing_ci[
        config.kite.connectivity.tube_line_indices
    ]:
        # if it is not a TE point
        if (
            idx
            not in config.kite.connectivity.wing_ci[
                config.kite.connectivity.te_line_indices
            ]
            and idx
            not in config.kite.connectivity.wing_cj[
                config.kite.connectivity.te_line_indices
            ]
            and idx not in config.kite.connectivity.plate_point_indices
        ):  # and not a plate point
            points_attachment_inbetween.append(config.kite.points_ini[idx])
            points_index_attachment_inbetween.append(idx)

    # need to sort these points per strut
    points_between_dict = {}
    tol_linepoint = 0.2
    for plate_index in config.kite.connectivity.plate_point_indices:
        left_le = config.kite.points_ini[plate_index[0]]
        left_te = config.kite.points_ini[plate_index[3]]
        for idx, point in enumerate(points_attachment_inbetween):
            if (
                abs(
                    structural_model.distance_point_to_line(point, [left_le, left_te])[
                        0
                    ]
                )
                < tol_linepoint
            ):
                points_between_dict[str(points_index_attachment_inbetween[idx])] = [
                    plate_index[0],
                    plate_index[3],
                ]

    # treating last strut separately
    right_le = config.kite.points_ini[
        config.kite.connectivity.plate_point_indices[-1][1]
    ]
    right_te = config.kite.points_ini[
        config.kite.connectivity.plate_point_indices[-1][2]
    ]
    for idx, point in enumerate(points_attachment_inbetween):
        if (
            abs(structural_model.distance_point_to_line(point, [right_le, right_te])[0])
            < tol_linepoint
        ):
            points_between_dict[str(points_index_attachment_inbetween[idx])] = [
                config.kite.connectivity.plate_point_indices[-1][1],
                config.kite.connectivity.plate_point_indices[-1][2],
            ]

    # the indices are sorted from left-to-right
    # the indices_mid are even sorted from TE to LE
    indices_le = [int(value[0]) for value in points_between_dict.values()]
    indices_mid = [int(key) for key in points_between_dict.keys()]
    indices_te = [int(value[1]) for value in points_between_dict.values()]

    return points_between_dict


def extract_rotational_resistances_dicts(points_between_dict, config):
    ##TODO: should actually be able to split up the LE and Struts
    # the indices are sorted from left-to-right
    # the indices_mid are even sorted from TE to LE
    # and all have the same length
    indices_le = [int(value[0]) for value in points_between_dict.values()]
    indices_mid = [int(key) for key in points_between_dict.keys()]
    indices_te = [int(value[1]) for value in points_between_dict.values()]

    le_rotational_resistance_dict = {}
    # making unique copies
    indices_le_unique = np.unique(np.copy(indices_le))
    indices_le_unique = sorted(
        indices_le_unique, key=lambda idx: -config.kite.points_ini[idx][1]
    )
    indices_mid_unique = np.unique(np.copy(indices_mid))
    indices_te_unique = np.unique(np.copy(indices_te))
    indices_te_unique = sorted(
        indices_te_unique, key=lambda idx: -config.kite.points_ini[idx][1]
    )

    # first point
    le_rotational_resistance_dict[str(indices_le_unique[0])] = [
        indices_le_unique[1],
        indices_te_unique[0],
    ]
    # middle points
    for idx in range(1, len(indices_le_unique) - 1):
        le_rotational_resistance_dict[str(indices_le_unique[idx])] = [
            indices_le_unique[idx - 1],
            indices_le_unique[idx + 1],
        ]
    # last point
    le_rotational_resistance_dict[str(indices_le_unique[-1])] = [
        indices_le_unique[-2],
        indices_te_unique[-1],
    ]

    # TODO: continue with the middle-points
    indices_le_sorted = sorted(
        indices_le, key=lambda idx: -config.kite.points_ini[idx][1]
    )
    indices_te_sorted = sorted(
        indices_te, key=lambda idx: -config.kite.points_ini[idx][1]
    )
    strut_rotational_resistance_dict = {}

    ## loop through all the le points, knowing that when there is no more te point, the strut is finished
    # the indices_mid are aranged from left-to-right and from TE-to-LE
    for i, index in enumerate(indices_le_sorted[:-1]):
        # TODO: implement some logic, something like below
        # if has another point back and front
        if (
            indices_le_sorted[i] == indices_le_sorted[i + 1]
            and indices_le_sorted[i] == indices_le_sorted[i - 1]
        ):
            strut_rotational_resistance_dict[str(indices_mid[i])] = [
                indices_mid[i - 1],
                indices_mid[i + 1],
            ]
        elif indices_le_sorted[i] == indices_le_sorted[i + 1]:
            strut_rotational_resistance_dict[str(indices_mid[i])] = [
                indices_te_sorted[i],
                indices_mid[i + 1],
            ]
        elif indices_le_sorted[i] == indices_le_sorted[i - 1]:
            strut_rotational_resistance_dict[str(indices_mid[i])] = [
                indices_mid[i - 1],
                indices_le_sorted[i],
            ]
        elif (
            indices_le_sorted[i] != indices_le_sorted[i + 1]
            and indices_le_sorted[i] != indices_le_sorted[i - 1]
        ):
            strut_rotational_resistance_dict[str(indices_mid[i])] = [
                indices_te_sorted[i],
                indices_le_sorted[i],
            ]

    return le_rotational_resistance_dict, strut_rotational_resistance_dict


def initialize_bending_spring(k_bend, init_cond, params, b_conn, rotational_dict):
    # Identify indices of the last spring
    # to be able to add these as keys to the new bending_params
    idx_last_normal_spring = len(params["k"]) - 1
    idx_of_this_spring = idx_last_normal_spring
    bending_params = {}

    # Define a rotational particle list
    # each node with resistance is an index/key
    # each value represents the links
    rotational_particles = rotational_dict

    # Looping over each key, and defining the surroundings
    # print(f'rotational_particles: {len(rotational_particles)} {rotational_particles}')
    for key in rotational_particles.keys():
        # finding the particles nearby
        p1_idx, p2_idx = rotational_particles[key]
        pe_idx = int(key)
        # print(f"pe_idx:{pe_idx}: [p1_idx:{p1_idx}, p2_idx:{p2_idx}]")
        # print(f"init_cond[p1_idx]:{init_cond[p1_idx]}")
        # extracting only the initial position
        p1 = np.array(init_cond[p1_idx][0])
        # extracting only the initial position
        p2 = np.array(init_cond[p2_idx][0])

        # defining the particles
        Pa = p1
        Pb = p2
        Pe = np.array(init_cond[int(key)][0])
        # print(f"Pa: {p1_idx} - {Pa}")
        # print(f"Pb: {p2_idx} - {Pb}")
        # print(f"Pe: {pe_idx} - {Pe}")

        # calculating the projection vector
        P_proj = Pa + np.dot(Pe - Pa, Pb - Pa) / np.dot(Pb - Pa, Pb - Pa) * (Pb - Pa)
        # distance from A and B to the projected point
        h_a = np.linalg.norm(P_proj - Pa)
        h_b = np.linalg.norm(P_proj - Pb)

        # defining the alpha's
        alpha_a = h_b / (h_a + h_b)
        alpha_b = h_a / (h_a + h_b)
        alpha_e = -1

        # making a dict of the bending_parameters
        # print(f'idx_of_this_spring: {idx_of_this_spring}')
        bending_params[f"{idx_of_this_spring}"] = [
            alpha_a,
            alpha_b,
            alpha_e,
            pe_idx,
        ]
        idx_of_this_spring += 1

        # TODO: this essentially reduces k_bend to mu by a factor of 10
        thickness_of_rotational_resistance_beam = (
            0.1  # represents the thickness of the strut tubes
        )
        # calculating lambda_bend, fed in as k[idx] for the rotational True springs
        lambda_bend = (
            (2 / 3)
            * ((h_a + h_b) / ((h_a * h_b) ** 2))
            * k_bend
            * thickness_of_rotational_resistance_beam
        )

        # increasing the length of the variables, by 1 - done for each additional rotational spring element
        b_conn = np.append(b_conn, [[p1_idx, p2_idx]], axis=0)
        params["k"] = np.append(params["k"], lambda_bend)
        params["l0"] = np.append(params["l0"], np.linalg.norm(p1 - p2))
        params["is_compression"] = np.append(params["is_compression"], False)
        params["is_tension"] = np.append(params["is_tension"], False)
        params["is_pulley"] = np.append(params["is_pulley"], False)
        params["is_rotational"] = np.append(params["is_rotational"], True)

    # adding the alpha_dict to params
    # print(f'params[l0]: {params["l0"]}')
    params["bending_params"] = bending_params
    # print(f'bending_params: {bending_params}')
    # print(f'b_conn: {b_conn}')

    return params
