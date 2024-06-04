import numpy as np


def calculate_edge_lengths(ci, cj, pos):
    """returns the edge lengths between the nodes with index ci and cj
    for the given positions pos
    input : ci,cj,pos
    output: springL"""
    springL = np.zeros(ci.shape)
    for idx, (ci, cj) in enumerate(zip(ci, cj)):
        springL[idx] = np.linalg.norm(pos[cj, :] - pos[ci, :])
    return springL


def update_for_billowing(points_ini, springL_wing, u_p):
    ##TODO: hardcoded, ideally this would be resolved by the structural solver
    """adjusts the spring lengths of the wing for billowing (hardcoded)
    for the given u_p, for reference visit J.Poland MSc Thesis:
    http://resolver.tudelft.nl/uuid:39d67249-53c9-47b4-84c0-ddac948413a5"""

    def billowing_up(u_p):
        """defines ratio's to multiply the TE lengths by
        from J.Poland MSc Thesis - http://resolver.tudelft.nl/uuid:39d67249-53c9-47b4-84c0-ddac948413a5
        """
        bill_plate_1 = -0.056
        bill_plate_2 = -0.0216
        bill_plate_3 = +0.01
        bill_plate_4 = +0  # outlier, thus 0

        return np.array(
            [
                1 + bill_plate_1 * (1 - u_p),
                1 + bill_plate_2 * (1 - u_p),
                1 + bill_plate_3 * (1 - u_p),
                1 + bill_plate_4 * (1 - u_p),
            ]
        )

    ##TODO: fix hardcoding of indices
    ## Adjusts the spring lengths (multiply by ratios) of the wing for billowing, hardcoded indices
    springL_wing[8] = springL_wing[8] * billowing_up(u_p)[3]  # 17 -> 18
    springL_wing[14] = springL_wing[14] * billowing_up(u_p)[2]  # 16->17
    springL_wing[20] = springL_wing[20] * billowing_up(u_p)[1]  # 15- >16
    springL_wing[26] = springL_wing[26] * billowing_up(u_p)[0]  # 14- >15
    springL_wing[32] = springL_wing[32] * billowing_up(u_p)[1]  # 13 to 14
    springL_wing[38] = springL_wing[38] * billowing_up(u_p)[2]  # 12 to 13
    springL_wing[44] = springL_wing[44] * billowing_up(u_p)[3]  # 11 to 12

    # also adjusts the inner points to make things symmetrical again
    points_ini[14, 1] = points_ini[14, 1] * (
        billowing_up(u_p)[0]
    )  # Move points too to avoid asymmetry
    points_ini[15, 1] = points_ini[15, 1] * (
        billowing_up(u_p)[0]
    )  # Move points too to avoid asymmetry

    return points_ini, springL_wing


# def make_symmetrical(points, tolerance=1e-5):

#     # loop through each point
#     for i,point_left in enumerate(points):

#         if point_left[1]>0: # if the point is on the left-side
#             # loop through each other point
#             for j,point_right in enumerate(points):

#                 # if the point is on the right-side and the mirror point of i
#                 # mirror it, by making it equal to the left-point, but negative in the y
#                 if point_right[1]<0 and i!=j and \
#                     np.isclose(point_left[0], point_right[0], atol=tolerance) and \
#                     np.isclose(point_left[2], point_right[2], atol=tolerance) and \
#                     np.isclose(point_left[1], -point_right[1], atol=tolerance):
#                         points[j] = [point_left[0],-point_left[1],point_left[2]]

#     return points

# def make_symmetrical(points, tolerance=1e-5):
#     mirrored_points = {}

#     for i, point in enumerate(points):
#         x, y, z = point
#         if y > 0:
#             mirrored_points[i] = [x, -y, z]
#         else:
#             mirrored_points[i] = point

#     symmetrical_points = [mirrored_points[i] for i in range(len(points))]
#     return np.array(symmetrical_points)


def find_symmetrical_pairs(points, tolerance=1e-5):
    left_points_indices = [i for i, point in enumerate(points) if point[1] > 0]
    right_points_indices = [i for i, point in enumerate(points) if point[1] < 0]

    symmetrical_pairs = []

    for left_index in left_points_indices:
        mirrored_left_point = [
            points[left_index][0],
            -points[left_index][1],
            points[left_index][2],
        ]
        matched_indices = [
            right_index
            for right_index in right_points_indices
            if np.allclose(points[right_index], mirrored_left_point, atol=tolerance)
        ]
        symmetrical_pairs.extend(
            [(left_index, right_index) for right_index in matched_indices]
        )

    return symmetrical_pairs


def make_symmetrical(points_ini, points, tolerance=1e-5):

    pair_indices = find_symmetrical_pairs(points_ini, tolerance=1e-5)

    for (
        pair_index
    ) in (
        pair_indices
    ):  # loop through each pair, and make [1] equal to the mirror of [0], using -y coordinate
        point_to_be_mirrored = points[pair_index[0]]
        points[pair_index[1]] = [
            point_to_be_mirrored[0],
            -point_to_be_mirrored[1],
            point_to_be_mirrored[2],
        ]

    return points
