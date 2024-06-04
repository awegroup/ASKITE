import numpy as np


def reorder_wing_nodes(points_ini):
    """Takes points_ini (initial nodes), and re-orders them
    from LE-to-TE and from right-to-left, such that:
    [0] and [1] form the most-right chordwise line
    (hardcoded)
    """

    coord = np.empty((20, 3))

    coord[0, :] = points_ini[20, :]
    coord[1, :] = points_ini[10, :]

    coord[2, :] = points_ini[9, :]
    coord[3, :] = points_ini[11, :]

    coord[4, :] = points_ini[8, :]
    coord[5, :] = points_ini[12, :]

    coord[6, :] = points_ini[7, :]
    coord[7, :] = points_ini[13, :]

    coord[8, :] = points_ini[6, :]
    coord[9, :] = points_ini[14, :]

    coord[10, :] = points_ini[5, :]
    coord[11, :] = points_ini[15, :]

    coord[12, :] = points_ini[4, :]
    coord[13, :] = points_ini[16, :]

    coord[14, :] = points_ini[3, :]
    coord[15, :] = points_ini[17, :]

    coord[16, :] = points_ini[2, :]
    coord[17, :] = points_ini[18, :]

    coord[18, :] = points_ini[19, :]
    coord[19, :] = points_ini[1, :]

    return coord


# #defining the coordinates of the corner points of the structural mesh panels
# def order_struc_nodes_right_to_left(points_ini,plate_point_indices):
#     # defined from LE-to-TE and from right-to-left ##TODO: implement this from left-to-right
#     ##TODO: plates should be named segments...
#     # we can get these coordinates from plate_points_indices
#     # orderded from left-to-right like: [idx_left_LE, idx_right_LE, idx_right_TE,idx_left_TE]

#     plate_point_indices_R_to_L = np.flipud(plate_point_indices) #flip up-down

#     coordinates_plates_LE_to_TE_R_to_L = []
#     for segment_indices in plate_point_indices_R_to_L: # loop through each segment
#         # append right-LE point (segment_indices[1]) and right-TE point [segment_indices[2]]
#         coordinates_plates_LE_to_TE_R_to_L.append(points_ini[segment_indices[1]])
#         coordinates_plates_LE_to_TE_R_to_L.append(points_ini[segment_indices[2]])

#     # appending the last left edge (as now all right segment edges have been appended)
#     coordinates_plates_LE_to_TE_R_to_L.append(points_ini[plate_point_indices[0][0]]) #LE
#     coordinates_plates_LE_to_TE_R_to_L.append(points_ini[plate_point_indices[0][3]]) #TE

#     # rename and make into array
#     coordinates_struc = np.array(coordinates_plates_LE_to_TE_R_to_L)

#     return coordinates_struc


# defining the coordinates of the corner points of the structural mesh panels
def order_struc_nodes_right_to_left(points_ini, plate_point_indices):
    # defined from LE-to-TE and from right-to-left ##TODO: implement this from left-to-right
    ##TODO: plates should be named segments...
    # we can get these coordinates from plate_points_indices
    # orderded from left-to-right like: [idx_left_LE, idx_right_LE, idx_right_TE,idx_left_TE]

    plate_point_indices_R_to_L = np.flipud(plate_point_indices)

    coordinates_plates_LE_to_TE_R_to_L = np.empty((0, 3))  # Initialize an empty array

    for segment_indices in plate_point_indices_R_to_L:
        coordinates_plates_LE_to_TE_R_to_L = np.vstack(
            (
                coordinates_plates_LE_to_TE_R_to_L,
                points_ini[segment_indices[1]],
                points_ini[segment_indices[2]],
            )
        )

    coordinates_plates_LE_to_TE_R_to_L = np.vstack(
        (
            coordinates_plates_LE_to_TE_R_to_L,
            points_ini[plate_point_indices[0][0]],
            points_ini[plate_point_indices[0][3]],
        )
    )

    return coordinates_plates_LE_to_TE_R_to_L


def create_geometry_LEI(coordinates, vel_app, N, RING_GEOMETRY, MODEL):

    filaments = []
    controlpoints = []
    rings = []
    wingpanels = []
    ringvec = []
    coord_L = []

    ##TODO: rename these to something more descriptive (assesement point, orientation determination point etc.)
    chord_1_4 = 1 / 4
    chord_3_4 = 3 / 4

    for i in range(N - 1):

        section = {
            "p1": coordinates[2 * i, :],
            "p2": coordinates[2 * i + 2, :],
            "p3": coordinates[2 * i + 3, :],
            "p4": coordinates[2 * i + 1, :],
        }
        wingpanels.append(section)
        di = np.linalg.norm(
            coordinates[2 * i, :] * chord_3_4
            + coordinates[2 * i + 1, :] * chord_1_4
            - (
                coordinates[2 * i + 2, :] * chord_3_4
                + coordinates[2 * i + 3, :] * chord_1_4
            )
        )
        if i == 0:
            diplus = np.linalg.norm(
                coordinates[2 * (i + 1), :] * chord_3_4
                + coordinates[2 * (i + 1) + 1, :] * chord_1_4
                - (
                    coordinates[2 * (i + 1) + 2, :] * chord_3_4
                    + coordinates[2 * (i + 1) + 3, :] * chord_1_4
                )
            )
            ncp = di / (di + diplus)
        elif i == N - 2:
            dimin = np.linalg.norm(
                coordinates[2 * (i - 1), :] * chord_3_4
                + coordinates[2 * (i - 1) + 1, :] * chord_1_4
                - (
                    coordinates[2 * (i - 1) + 2, :] * chord_3_4
                    + coordinates[2 * (i - 1) + 3, :] * chord_1_4
                )
            )
            ncp = dimin / (dimin + di)
        else:
            dimin = np.linalg.norm(
                coordinates[2 * (i - 1), :] * chord_3_4
                + coordinates[2 * (i - 1) + 1, :] * chord_1_4
                - (
                    coordinates[2 * (i - 1) + 2, :] * chord_3_4
                    + coordinates[2 * (i - 1) + 3, :] * chord_1_4
                )
            )
            diplus = np.linalg.norm(
                coordinates[2 * (i + 1), :] * chord_3_4
                + coordinates[2 * (i + 1) + 1, :] * chord_1_4
                - (
                    coordinates[2 * (i + 1) + 2, :] * chord_3_4
                    + coordinates[2 * (i + 1) + 3, :] * chord_1_4
                )
            )
            ncp = chord_1_4 * (dimin / (dimin + di) + di / (di + diplus) + 1)

        ncp = 1 - ncp
        chord = np.linalg.norm(
            (section["p2"] + section["p1"]) / 2 - (section["p3"] + section["p4"]) / 2
        )
        LLpoint = (section["p2"] * (1 - ncp) + section["p1"] * ncp) * chord_3_4 + (
            section["p3"] * (1 - ncp) + section["p4"] * ncp
        ) * chord_1_4
        VSMpoint = (section["p2"] * (1 - ncp) + section["p1"] * ncp) * chord_1_4 + (
            section["p3"] * (1 - ncp) + section["p4"] * ncp
        ) * chord_3_4
        coord_L.append(LLpoint)

        # Define bound vortex filament
        bound = {
            "id": "bound",
            "x1": section["p1"] * chord_3_4 + section["p4"] * chord_1_4,
            "x2": section["p2"] * chord_3_4 + section["p3"] * chord_1_4,
            "gamma": 0,
        }
        filaments.append(bound)

        x_airf = np.cross(VSMpoint - LLpoint, section["p2"] - section["p1"])
        x_airf = x_airf / np.linalg.norm(x_airf)
        y_airf = VSMpoint - LLpoint
        y_airf = y_airf / np.linalg.norm(y_airf)
        z_airf = bound["x2"] - bound["x1"]
        # z_airf[0] = 0
        z_airf = z_airf / np.linalg.norm(z_airf)
        airf_coord = np.column_stack([x_airf, y_airf, z_airf])

        normal = x_airf
        tangential = y_airf
        if MODEL == "VSM":
            cp = {
                "coordinates": VSMpoint,
                "chord": chord,
                "normal": normal,
                "tangential": tangential,
                "airf_coord": airf_coord,
                "coordinates_aoa": LLpoint,
            }
            controlpoints.append(cp)
        elif MODEL == "LLT":

            cp = {
                "coordinates": LLpoint,
                "chord": chord,
                "normal": normal,
                "tangential": tangential,
                "airf_coord": airf_coord,
            }
            controlpoints.append(cp)

        temp = {
            "r0": bound["x2"] - bound["x1"],
            "r1": cp["coordinates"] - bound["x1"],
            "r2": cp["coordinates"] - bound["x2"],
            "r3": cp["coordinates"] - (bound["x2"] + bound["x1"]) / 2,
        }
        ringvec.append(temp)

        temp = vel_app / np.linalg.norm(vel_app)
        if RING_GEOMETRY == "3fil":
            # create trailing filaments, at x1 of bound filament
            temp1 = {"dir": temp, "id": "trailing_inf1", "x1": bound["x1"], "gamma": 0}
            filaments.append(temp1)

            # create trailing filaments, at x2 of bound filament
            temp1 = {"x1": bound["x2"], "dir": temp, "id": "trailing_inf2", "gamma": 0}
            filaments.append(temp1)
        elif RING_GEOMETRY == "5fil":
            temp1 = {
                "x1": section["p4"],
                "x2": bound["x1"],
                "gamma": 0,
                "id": "trailing1",
            }
            filaments.append(temp1)

            temp1 = {
                "dir": temp,
                "id": "trailing_inf1",
                "x1": section["p4"],
                "gamma": 0,
            }
            filaments.append(temp1)

            # create trailing filaments, at x2 of bound filament
            temp1 = {
                "x2": section["p3"],
                "x1": bound["x2"],
                "gamma": 0,
                "id": "trailing1",
            }
            filaments.append(temp1)

            temp1 = {
                "x1": section["p3"],
                "dir": temp,
                "id": "trailing_inf2",
                "gamma": 0,
            }
            filaments.append(temp1)

        #

        rings.append(filaments)
        filaments = []

    coord_L = np.array(coord_L)
    return controlpoints, rings, wingpanels, ringvec, coord_L


def refine_LEI_mesh(coord, N_sect, N_SPLITS):
    refined_coord = []

    for i_sec in range(N_sect):
        temp_coord = np.empty((int(N_SPLITS * 2), 3))
        for i_spl in range(N_SPLITS):
            temp_coord[2 * i_spl] = (
                coord[2 * i_sec, :] * (N_SPLITS - i_spl) / N_SPLITS
                + coord[2 * (i_sec + 1), :] * (i_spl) / N_SPLITS
            )
            temp_coord[2 * i_spl + 1] = (
                coord[2 * i_sec + 1, :] * (N_SPLITS - i_spl) / N_SPLITS
                + coord[2 * (i_sec + 1) + 1, :] * (i_spl) / N_SPLITS
            )
        if i_sec == 0:
            refined_coord = temp_coord
        else:
            refined_coord = np.append(refined_coord, temp_coord, axis=0)

    refined_coord = np.append(
        refined_coord, [coord[2 * N_sect, :], coord[2 * N_sect + 1, :]], axis=0
    )

    return refined_coord


import numpy as np


def refine_LEI_mesh(coord, N_sect, N_SPLITS):
    refined_coord = np.empty((0, 3))

    for i_sec in range(N_sect):
        temp_coord = np.empty((int(N_SPLITS * 2), 3))
        for i_spl in range(N_SPLITS):
            temp_coord[2 * i_spl] = (
                coord[2 * i_sec, :] * (N_SPLITS - i_spl) / N_SPLITS
                + coord[2 * (i_sec + 1), :] * (i_spl) / N_SPLITS
            )
            temp_coord[2 * i_spl + 1] = (
                coord[2 * i_sec + 1, :] * (N_SPLITS - i_spl) / N_SPLITS
                + coord[2 * (i_sec + 1) + 1, :] * (i_spl) / N_SPLITS
            )
        refined_coord = np.append(refined_coord, temp_coord, axis=0)

    refined_coord = np.append(
        refined_coord, [coord[2 * N_sect, :], coord[2 * N_sect + 1, :]], axis=0
    )

    return refined_coord


def refine_LEI_mesh_billowing(wingpanels, bill_angles, N_SPLITS):
    refined_coord = []
    for i_sec in range(len(wingpanels)):
        angle = np.deg2rad(bill_angles[i_sec])
        L_sec1 = np.linalg.norm(wingpanels[i_sec]["p2"] - wingpanels[i_sec]["p1"])
        R1 = L_sec1 / 2 / np.sin(angle)
        L_sec2 = np.linalg.norm(wingpanels[i_sec]["p3"] - wingpanels[i_sec]["p4"])
        R2 = L_sec2 / 2 / np.sin(angle)
        zvec = (wingpanels[i_sec]["p2"] + wingpanels[i_sec]["p1"]) / 2 - (
            wingpanels[i_sec]["p4"] + wingpanels[i_sec]["p3"]
        ) / 2
        zvec = zvec / np.linalg.norm(zvec)

        xvec1 = wingpanels[i_sec]["p2"] - wingpanels[i_sec]["p1"]
        xvec1 = xvec1 / np.linalg.norm(xvec1)
        yvec1 = np.cross(zvec, xvec1)
        yvec1 = yvec1 / np.linalg.norm(yvec1)

        xvec2 = wingpanels[i_sec]["p3"] - wingpanels[i_sec]["p4"]
        xvec2 = xvec2 / np.linalg.norm(xvec2)
        yvec2 = np.cross(zvec, xvec2)
        yvec2 = yvec2 / np.linalg.norm(yvec2)

        if i_sec > 4:
            xvec1 = wingpanels[i_sec]["p1"] - wingpanels[i_sec]["p2"]
            xvec1 = xvec1 / np.linalg.norm(xvec1)
            xvec2 = wingpanels[i_sec]["p4"] - wingpanels[i_sec]["p3"]
            xvec2 = xvec2 / np.linalg.norm(xvec2)

        xloc1 = np.linspace(-L_sec1 / 2, L_sec1 / 2, N_SPLITS)
        yloc01 = np.sqrt(R1**2 - (L_sec1 / 2) ** 2)
        yloc1 = -np.sqrt(R1**2 - xloc1**2) + yloc01
        zloc1 = np.zeros(N_SPLITS)

        xloc2 = np.linspace(-L_sec2 / 2, L_sec2 / 2, N_SPLITS)
        yloc02 = np.sqrt(R2**2 - (L_sec2 / 2) ** 2)
        yloc2 = -np.sqrt(R2**2 - xloc2**2) + yloc02
        zloc2 = np.zeros(N_SPLITS)

        vec1 = np.array([xvec1, yvec1, zvec]).T
        vec2 = np.array([xvec2, yvec2, zvec]).T
        ax_pos1 = (wingpanels[i_sec]["p2"] + wingpanels[i_sec]["p1"]) / 2
        ax_pos2 = (wingpanels[i_sec]["p3"] + wingpanels[i_sec]["p4"]) / 2
        temp_coord = np.empty((int(N_SPLITS * 2), 3))
        for i_spl in range(N_SPLITS):
            coord_loc1 = np.array([xloc1[i_spl], yloc1[i_spl], zloc1[i_spl]])
            coord_loc2 = np.array([xloc2[i_spl], yloc2[i_spl], zloc2[i_spl]])
            coord1 = np.matmul(vec1, coord_loc1) + ax_pos1
            coord2 = np.matmul(vec2, coord_loc2) + ax_pos2

            if i_sec > 4:
                ind = 2 * N_SPLITS - 1 - (2 * i_spl + 1)
                temp_coord[ind] = coord1
                ind = 2 * N_SPLITS - 1 - 2 * i_spl
                temp_coord[ind] = coord2
            else:
                temp_coord[2 * i_spl] = coord1
                temp_coord[2 * i_spl + 1] = coord2

        if i_sec == 0:
            refined_coord = temp_coord
        else:
            refined_coord = np.append(refined_coord, temp_coord[2::, :], axis=0)

    return refined_coord


def struc2aero_old(
    points_ini,
    vel_app,
    n_chordwise_elements,
    plate_point_indices,
    BILLOWING_ANGLES,
    AERO_CONFIG,
):
    """
    This function transform the structural mesh into an aerodynamic mesh
    and prepares the input parameters for the VSM (/LLT)

    It uses the following functions:
        reorder_wing_nodes
        create_geometry_LEI
        refine_LEI_mesh_billowing

    input:
        points_ini  : initial coordinates used in simulation
        vel_app        : freestream velocity
        RING_GEOMETRY    : 3fil or 5fil ##TODO: what are these?
        MODEL       : VSM or LLT
        N_SPLITS     : number of times each segment is splitted for the aero-refinement

    output:
        coord_refined                   : refined coordinates of the wing ordered from LE-to-TE and right-to-left
        controlpoints                   : ##TODO: what are these?
        rings                           : ##TODO: what are these?
        wingpanels                      : ##TODO: what are these?
        ringvec                         : ##TODO: what are these?
        coord_L                         : ##TODO: what are these?
        tubediams_refined               : the tube diameters of each refined element
        canopyheights_refined           : the canopy height of each refined element
        n_chordwise_elements_refined    : the number of chordwise elements (each chordwise_element consists of 2 nodes)
        N_segments                      : the number of segments (canopy pieces between struts)
    """
    RING_GEOMETRY = AERO_CONFIG["RING_GEOMETRY"]
    MODEL = AERO_CONFIG["MODEL"]
    N_SPLITS = AERO_CONFIG["N_SPLITS"]

    # Find the coordinates of the wing in the right-order (LE-to-TE, right-to-left)
    # hardcoded for V3 coord = reorder_wing_nodes(points_ini)
    coord = order_struc_nodes_right_to_left(points_ini, plate_point_indices)

    ##TODO: understand this bit of code
    # ? what is coord_L
    # ? what is ringvec
    # ? what is wingpanels
    # "define a system of vorticity"
    controlpoints, rings, wingpanels, ringvec, coord_L = create_geometry_LEI(
        coord, vel_app, n_chordwise_elements, RING_GEOMETRY, MODEL
    )

    # Find the coordinates of the refined mesh, include billowing
    coord_refined = refine_LEI_mesh_billowing(
        wingpanels, BILLOWING_ANGLES, N_SPLITS + 1
    )

    # define the number of refined chordwise elements (each chordwise_element consists of 2 nodes)
    n_chordwise_elements_refined = int(len(coord_refined) / 2)

    # define the number of refined segments {type: int}
    N_segments_refined = int(n_chordwise_elements_refined - 1)

    # define the refined system of vorticity
    controlpoints, rings, wingpanels, ringvec, coord_L = create_geometry_LEI(
        coord_refined, vel_app, n_chordwise_elements_refined, RING_GEOMETRY, MODEL
    )

    return (
        coord_refined,
        controlpoints,
        rings,
        wingpanels,
        ringvec,
        coord_L,
        n_chordwise_elements_refined,
        N_segments_refined,
    )


def struc2aero(
    points_ini,
    vel_app,
    n_chordwise_elements,
    plate_point_indices,
    BILLOWING_ANGLES,
    AERO_CONFIG,
):
    """
    This function transform the structural mesh into an aerodynamic mesh
    and prepares the input parameters for the VSM (/LLT)

    It uses the following functions:
        reorder_wing_nodes
        create_geometry_LEI
        refine_LEI_mesh_billowing

    input:
        points_ini  : initial coordinates used in simulation
        vel_app        : freestream velocity
        RING_GEOMETRY    : 3fil or 5fil ##TODO: what are these?
        MODEL       : VSM or LLT
        N_SPLITS     : number of times each segment is splitted for the aero-refinement

    output:
        coord_refined                   : refined coordinates of the wing ordered from LE-to-TE and right-to-left
        controlpoints                   : ##TODO: what are these?
        rings                           : ##TODO: what are these?
        wingpanels                      : ##TODO: what are these?
        ringvec                         : ##TODO: what are these?
        coord_L                         : ##TODO: what are these?
        tubediams_refined               : the tube diameters of each refined element
        canopyheights_refined           : the canopy height of each refined element
        n_chordwise_elements_refined    : the number of chordwise elements (each chordwise_element consists of 2 nodes)
        N_segments                      : the number of segments (canopy pieces between struts)
    """
    RING_GEOMETRY = AERO_CONFIG["RING_GEOMETRY"]
    MODEL = AERO_CONFIG["MODEL"]
    N_SPLITS = AERO_CONFIG["N_SPLITS"]

    # Find the coordinates of the wing in the right-order (LE-to-TE, right-to-left)
    # hardcoded for V3 coord = reorder_wing_nodes(points_ini)
    coord = order_struc_nodes_right_to_left(points_ini, plate_point_indices)

    ##TODO: understand this bit of code
    # ? what is coord_L
    # ? what is ringvec
    # ? what is wingpanels
    # "define a system of vorticity"
    controlpoints, rings, wingpanels, ringvec, coord_L = create_geometry_LEI(
        coord, vel_app, n_chordwise_elements, RING_GEOMETRY, MODEL
    )

    # Find the coordinates of the refined mesh, include billowing
    coord_refined = refine_LEI_mesh_billowing(
        wingpanels, BILLOWING_ANGLES, N_SPLITS + 1
    )

    # define the number of refined chordwise elements (each chordwise_element consists of 2 nodes)
    n_chordwise_elements_refined = int(len(coord_refined) / 2)

    # define the number of refined segments {type: int}
    N_segments_refined = int(n_chordwise_elements_refined - 1)

    # define the refined system of vorticity
    controlpoints, rings, wingpanels, ringvec, coord_L = create_geometry_LEI(
        coord_refined, vel_app, n_chordwise_elements_refined, RING_GEOMETRY, MODEL
    )

    return (
        coord_refined,
        controlpoints,
        rings,
        wingpanels,
        ringvec,
        coord_L,
        n_chordwise_elements_refined,
        N_segments_refined,
    )


### NEW METHOD post 2024/02
def extract_wingpanel_corners_aero_orderded(points_ini, plate_point_indices):
    # defined from LE-to-TE and from right-to-left
    # TODO: implement this from left-to-right
    # TODO: plates should be named segments...
    # We can get these coordinates from plate_point_indices
    # Ordered from left-to-right like: [idx_left_LE, idx_right_LE, idx_right_TE,idx_left_TE]

    plate_point_indices_R_to_L = np.flipud(plate_point_indices)
    coordinates_plates_LE_to_TE_R_to_L = np.empty((0, 3))  # Initialize an empty array
    index_transformation = []

    for segment_indices in plate_point_indices_R_to_L:
        coordinates_plates_LE_to_TE_R_to_L = np.vstack(
            (
                coordinates_plates_LE_to_TE_R_to_L,
                points_ini[segment_indices[1]],
                points_ini[segment_indices[2]],
            )
        )
        # Keep track of the indices transformation
        index_transformation.extend([segment_indices[1], segment_indices[2]])

    coordinates_plates_LE_to_TE_R_to_L = np.vstack(
        (
            coordinates_plates_LE_to_TE_R_to_L,
            points_ini[plate_point_indices[0][0]],
            points_ini[plate_point_indices[0][3]],
        )
    )
    # Add the indices of the last two points
    index_transformation.extend([plate_point_indices[0][0], plate_point_indices[0][3]])

    return coordinates_plates_LE_to_TE_R_to_L, np.array(index_transformation)
