import numpy as np


def extract_pulley_connectivity(points, bridle_ci, bridle_cj, pulley_data):
    PULLEY_point_indices = pulley_data["point_indices"]
    number_of_pulleys_in_back_lines = pulley_data["number_of_pulleys_in_back_lines"]

    # print(f'PULLEY_point_indices {len(PULLEY_point_indices)}, number_of_pulleys_in_back_lines: {number_of_pulleys_in_back_lines}')
    # Sorting pulleys based on z-coordinate
    sorted_pulley_point_indices = sorted(
        PULLEY_point_indices, key=lambda index: points[index][2]
    )

    # TODO: this is where you input the number of pulleys in the back-lines, that have a different feature
    # for pulleys where the line is below the pulley point
    pulley_point_indices_line_below = sorted_pulley_point_indices[
        :number_of_pulleys_in_back_lines
    ]
    if (
        len(PULLEY_point_indices) > number_of_pulleys_in_back_lines
    ):  # if there are pulleys in the front lines
        pulley_point_indices_line_above = sorted_pulley_point_indices[
            number_of_pulleys_in_back_lines:
        ]
    else:
        pulley_point_indices_line_above = []  # this should be the case for the V3.25

    all_possible_line_indices = []

    for idx, (idx_bridle_node_i, idx_bridle_node_j) in enumerate(
        zip(bridle_ci, bridle_cj)
    ):  # loop through each bridle line
        # if the current line its index i OR j  is a pulley point
        # AND the other line is LOWER than the pulley_point, i.e. pulley line is BELOW the pulley
        if (
            idx_bridle_node_i in pulley_point_indices_line_below
            and points[idx_bridle_node_i][2] > points[idx_bridle_node_j][2]
        ) or (
            idx_bridle_node_j in pulley_point_indices_line_below
            and points[idx_bridle_node_j][2] > points[idx_bridle_node_i][2]
        ):
            all_possible_line_indices.append(idx)

        # if the current line its index i OR j  is a pulley point
        # AND the other line is HIGHER than the pulley_point, i.e. pulley line is ABOVE the pulley
        if (
            idx_bridle_node_i in pulley_point_indices_line_above
            and points[idx_bridle_node_i][2] < points[idx_bridle_node_j][2]
        ) or (
            idx_bridle_node_j in pulley_point_indices_line_above
            and points[idx_bridle_node_j][2] < points[idx_bridle_node_i][2]
        ):
            all_possible_line_indices.append(idx)

    # loop through all_possible_line_indices, twice to try and find the matching line pair
    pulley_line_pair_indices = {}
    pulley_line_indices = []

    for line_index_1 in all_possible_line_indices:
        for line_index_2 in all_possible_line_indices:
            # break if the same line is compared
            if line_index_1 == line_index_2:
                break

            # list the point_indices of the two lines
            point_indices = [
                bridle_ci[line_index_1],
                bridle_cj[line_index_1],
                bridle_ci[line_index_2],
                bridle_cj[line_index_2],
            ]
            # sort the point_indices based on the z-coordinate of the points[index]
            sorted_indices = sorted(point_indices, key=lambda index: points[index][2])

            # IF the pulley line is BELOW the pulley point, i.e. IF LAST two indices correspond to the same point
            # AND the index is in the pulley_index list (i.e. it is actually a pulley)
            # AND the line is not already used, i.e. not in the used_pulley_line_indices set
            if (
                sorted_indices[2] == sorted_indices[3]
                and sorted_indices[3] in pulley_point_indices_line_below
                and line_index_1 not in pulley_line_indices
                and line_index_2 not in pulley_line_indices
            ):
                # Append the line indices, to the pulley_line_indices
                pulley_line_indices.append(line_index_1)
                pulley_line_indices.append(line_index_2)
                # Make new key indices for the pulley line pair indices
                pulley_line_pair_indices[str(line_index_1)] = line_index_2

            # IF the pulley line is ABOVE the pulley point, i.e. IF FIRST indices correspond to the same point
            # AND the index is in the pulley_index list (i.e. it is actually a pulley)
            # AND the line is not already used, i.e. not in the used_pulley_line_indices set
            elif (
                sorted_indices[0] == sorted_indices[1]
                and sorted_indices[0] in pulley_point_indices_line_above
                and line_index_1 not in pulley_line_indices
                and line_index_2 not in pulley_line_indices
            ):
                # Append the line indices, to the pulley_line_indices
                pulley_line_indices.append(line_index_1)
                pulley_line_indices.append(line_index_2)
                # Make new key indices for the pulley line pair indices
                pulley_line_pair_indices[str(line_index_1)] = line_index_2

    # i know that line key and line value make up a pulley
    # so i want to plot both line key and line value
    # extract them separately and append
    pulley_ci_key = [
        bridle_ci[int(key_index)] for key_index in pulley_line_pair_indices.keys()
    ]
    pulley_cj_key = [
        bridle_cj[int(key_index)] for key_index in pulley_line_pair_indices.keys()
    ]
    pulley_ci_value = [
        bridle_ci[int(value_index)] for value_index in pulley_line_pair_indices.values()
    ]
    pulley_cj_value = [
        bridle_cj[int(value_index)] for value_index in pulley_line_pair_indices.values()
    ]
    pulley_ci = np.concatenate((pulley_ci_key, pulley_ci_value))
    pulley_cj = np.concatenate((pulley_cj_key, pulley_cj_value))

    pulley_data["line_indices"] = np.ndarray.flatten(np.array(pulley_line_indices))
    pulley_data["line_pair_indices"] = pulley_line_pair_indices
    pulley_data["ci"] = pulley_ci
    pulley_data["cj"] = pulley_cj

    additional_dict = {}
    for i in range(0, len(pulley_data["line_indices"]), 2):
        # line 1, key: line 1 and value: line 2 data
        line_1_key = str(pulley_data["line_indices"][i])
        line_1_idx_p3 = bridle_ci[pulley_data["line_indices"][i + 1]]
        line_1_idx_p4 = bridle_cj[pulley_data["line_indices"][i + 1]]
        line_1_rest_length_p3p4 = np.linalg.norm(
            points[line_1_idx_p3] - points[line_1_idx_p4]
        )
        additional_dict[line_1_key] = np.array(
            [line_1_idx_p3, line_1_idx_p4, line_1_rest_length_p3p4]
        )

        # line 2, key: line 2 and value: line 1 data
        line_2_key = str(pulley_data["line_indices"][i + 1])
        line_2_idx_p3 = bridle_ci[pulley_data["line_indices"][i]]
        line_2_idx_p4 = bridle_cj[pulley_data["line_indices"][i]]
        line_2_rest_length_p3p4 = np.linalg.norm(
            points[line_2_idx_p3] - points[line_2_idx_p4]
        )
        additional_dict[line_2_key] = np.array(
            [line_2_idx_p3, line_2_idx_p4, line_2_rest_length_p3p4]
        )

    pulley_data["other_line_pair"] = additional_dict

    return pulley_data
