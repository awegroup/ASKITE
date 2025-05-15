import numpy as np
from kitesim.initialisation.path_functions import load_module_from_path

# def extract_points(simulation_directory):
#     # """ reads out the surfplan file and if not found, raises an error """
#     filepath = f"{simulation_directory}/src/initialisation/V9_60C/results_bridle_only/"

#     if os.path.exists(filepath):
#         discretized_kite = np.load(filepath)
#     else:
#         raise Exception(f"Error: no file found, with filepath: {filepath}.")
#     return discretized_kite


def extract_te_line_indices(plate_point_indices, wing_ci, wing_cj):
    ##TODO: same definition is used in other analyze_surfplan file, should thus be in separate file
    """Extracts a list of TE_indices from the plate_point_indices
    using its order: [left_LE, right_LE, right_TE,left_TE]"""
    TE_point_indices = []
    for plate in plate_point_indices:
        TE_point_indices.append(plate[2])
        TE_point_indices.append(plate[3])

    TE_line_indices, tube_line_indices = [], []
    for idx, (ci, cj) in enumerate(zip(wing_ci, wing_cj)):
        if ci in TE_point_indices and cj in TE_point_indices:
            TE_line_indices.append(idx)
    # return [2,8,14,20,26,32,38,44,50] #old hardcoded V3
    return list(set(TE_line_indices))


def extract_points_and_connectivity(folder_path_kite_data, surfplan_file):
    points_struc = np.load(f"{folder_path_kite_data}/points.npy")
    bridle_ci = np.load(f"{folder_path_kite_data}/bridle_ci.npy")
    bridle_cj = np.load(f"{folder_path_kite_data}/bridle_cj.npy")
    plate_point_indices = np.load(f"{folder_path_kite_data}/plate_point_indices.npy")
    wing_ci = np.load(f"{folder_path_kite_data}/conn_wing_i.npy")
    wing_cj = np.load(f"{folder_path_kite_data}/conn_wing_j.npy")
    te_line_indices = extract_te_line_indices(plate_point_indices, wing_ci, wing_cj)
    tube_line_indices = np.load(f"{folder_path_kite_data}/tube_line_indices.npy")
    tube_ci = np.load(f"{folder_path_kite_data}/conn_tubular_frame_i.npy")
    tube_cj = np.load(f"{folder_path_kite_data}/conn_tubular_frame_j.npy")

    return (
        points_struc,
        bridle_ci,
        bridle_cj,
        plate_point_indices,
        wing_ci,
        wing_cj,
        te_line_indices,
        tube_line_indices,
    )


def extract_airfoil_geometry_data_from_ribs(
    rib_db_whole_model_struts, folder_path_initialisation, KITE_NAME
):
    """
    Returns the tubediam and max-height (b) of the LEI airfoil as defined by the
    Surfplan model. The thickness and camber are defined at the struts.

    Returns:
        tubediam_struc: tube thickness of the LEI airfoil
        canopyheight_struc: maximum height of the canopy wrt chord of the LEI airfoil

    {{ SideNote: Oriol Cayon states in his MSc thesis:
        "The camber of the LEI airfoil is defined as the distance
        from the LE tube to the maximum height
        divided by the chord (c)"
        Which would be: camb_list_struc= (b-D_tube)/c
        This is not implemented, instead the camb_list_struc= b/c is used }}

    """
    # print(f"folder_path_initialisation: {folder_path_initialisation}")
    # sys.path.append(folder_path_initialisation)  # adding the initialisation path
    # get_centroid = importlib.import_module(f"{KITE_NAME}.functions_wing").get_centroid
    # get_discretized_canopy = importlib.import_module(
    #     f"{KITE_NAME}.functions_wing"
    # ).get_discretized_canopy

    # # Define the path to the module
    # functions_wing_path = f"{folder_path_initialisation}/{KITE_NAME}/functions_wing.py"

    # # Load the module from the specified path
    # functions_wing_module = load_module_from_path(
    #     f"{KITE_NAME}.functions_wing", functions_wing_path
    # )

    # Extract the required functions
    get_centroid = load_module_from_path(
        KITE_NAME, f"{folder_path_initialisation}/functions_wing.py"
    ).get_centroid
    get_discretized_canopy = load_module_from_path(
        KITE_NAME, f"{folder_path_initialisation}/functions_wing.py"
    ).get_discretized_canopy

    tubediam_struc, canopyheight_struc = [], []
    D_tube, c, b = [], [], []
    for i in range(len(rib_db_whole_model_struts)):
        #   ##############################
        #   ##  c --> distance between TE and LE tube
        #
        #   # getting the LE most front points
        #   LE_front_i  = np.array(rib_db_whole_model_struts[i][4][0])
        #   LE_front_i1 = np.array(rib_db_whole_model_struts[i+1][4][0])
        #
        #   # getting the TE most back points
        #   TE_i = rib_db_whole_model_struts[i][4][5]
        #   TE_i1 = rib_db_whole_model_struts[i+1][4][5]
        #
        #   # getting the chord
        #   c_i = np.linalg.norm(np.array(LE_front_i) - np.array(TE_i))
        #   c_i1 = np.linalg.norm(np.array(LE_front_i1) - np.array(TE_i1))
        #   c_section = (c_i + c_i1)/2

        ##############################
        ## tubediam_struc

        # getting the LE tube points
        LE_tube_i = np.array(rib_db_whole_model_struts[i][1])
        #   LE_tube_i1 = np.array(rib_db_whole_model_struts[i+1][1])

        # getting the centers
        center_i = get_centroid(LE_tube_i)
        #   center_i1 =get_centroid(LE_tube_i1)

        # getting the diameter
        tubediam_i = 2 * np.mean(np.linalg.norm(LE_tube_i - center_i, axis=1))
        #   tubediam_i1   = 2* np.mean(np.linalg.norm(LE_tube_i1 - center_i1, axis=1))

        # getting the thickness and appending it
        #   D_tube_section = (tubediam_i + tubediam_i1)/2
        #   thicc_list_struc.append(D_tube_section/c_section)
        tubediam_struc.append(tubediam_i)

        ##############################
        ## b = "max height canopy wrt chord-line"

        # getting the canopy's
        canopy_i = np.array(rib_db_whole_model_struts[i][2])
        #   canopy_i1   = np.array(rib_db_whole_model_struts[i+1][2])

        # Find most front LE point and most aft TE point
        LE_point_i = np.array(rib_db_whole_model_struts[i][4][0])
        TE_point_i = rib_db_whole_model_struts[i][4][5]

        # add LE and TE
        canopy_i = np.vstack(([LE_point_i, TE_point_i], canopy_i))
        #   canopy_i1 = np.vstack(([LE_front_i1,TE_i1],canopy_i1))

        # defining the sorting condition, which is x over yz
        def get_angle_TE_to_LE(point):
            return np.arctan(point[0] / np.sqrt((point[1] ** 2) + ((point[2]) ** 2)))

        # sort from LE to TE
        canopy_i_sorted = np.array(
            sorted(canopy_i, key=get_angle_TE_to_LE, reverse=True)
        )
        #   canopy_i1_sorted = np.array(sorted(canopy_i1, key=get_angle_TE_to_LE, reverse=True))

        # interpolate to find more points
        interp_value_list = np.linspace(0, 1, 300)
        interp_values_i = np.array(
            get_discretized_canopy(canopy_i_sorted, interp_value_list)
        )
        #   interp_values_i1 = np.array(get_discretized_canopy(canopy_i1_sorted,interp_value_list))

        # defining a distance to line function
        def distance_point_to_line(point, line):
            # Calculate the equation of the line
            line_direction = line[1] - line[0]
            line_direction = line_direction / np.linalg.norm(line_direction)
            line_normal = np.cross(line_direction, [0, 0, 1])
            line_normal = line_normal / np.linalg.norm(line_normal)

            # Calculate the distance from the point to the line
            distance = np.dot(point - line[0], line_normal)

            return np.abs(distance)

        # get the max height
        canopyheight_i = np.amax(
            distance_point_to_line(
                interp_values_i,
                [LE_point_i, TE_point_i],
            )
        )
        #   b_i1 = np.amax(distance_point_to_line(interp_values_i1, [LE_front_i1,TE_i1]))
        canopyheight_struc.append(canopyheight_i)

        # get the camber and append it
        ##TODO: remove the correction factor
        #   b_corr_for_not_taking_mid = 1.2 #multiplying with a mid-canopy-billow correction factor
        #   b_section = ((b_i + b_i1)/2 ) * b_corr_for_not_taking_mid
        #   camb_list_struc.append(b_section/c_section)

        ## for checking
        #   D_tube.append((tubediam_i + tubediam_i1)/2)
        #   c.append((c_i + c_i1)/2)
        #   b.append((b_i + b_i1)/2)

    ## BILLOWING ANGLE
    ## From Oriol Cayon's MSc thesis
    # theta_e (billowing angle) = angle at which canopy leaves the strut (for V3 assumed = LE tube angle)
    #
    # I -  "If the pressure difference ∆P is assumed constant along the spanwise direction
    #       of the canopy, then the shape of the canopy billowing is defined by a semi-circular arc."
    # II - "If the membrane is assumed inelastic, the radius of curvature is not dependent on the magnitude of ∆P"
    #
    # I + II --> the spanwise shape of the canopy is assumed to be semi-circular
    # Under this assumption; the spanwise shape is found as f(ini_billowing_angle & the TE length)
    #
    # Which also means that when TE length increases, but canopy doesn't stretch, the canopy will be flatter
    #
    # TODO: WHAT ABOUT CHORDWISE VARIATIONS OF BILLOWING ANGLE?

    ##TODO: calc. instead of assume & define somewhere the N_segments_struc  instead of len(k)
    # Assume all ini billowing angles to be 5 deg, make a list as long as camb_list_struc
    billowing_angles = 5 * np.ones(
        len(tubediam_struc) - 1
    )  ##TODO: this should be computed instead

    return tubediam_struc, canopyheight_struc, billowing_angles


# def extract_rib_db_whole(folder_path_data_input, surfplan_file):
#     surfplan_discretized_path = f"{folder_path_data_input}/results_bridle_only/"
#     # conn_canopy_i = np.load(f'{surfplan_discretized_path}conn_canopy_i.npy')
#     # conn_canopy_j = np.load(f'{surfplan_discretized_path}conn_canopy_j.npy')
#     rib_db_whole_model_struts = np.load(
#         f"{surfplan_discretized_path}rib_db_whole_model_struts.npy", allow_pickle=True
#     )
#     return rib_db_whole_model_struts


# def extract_bridle_line_system_data(simulation_directory, surfplan_file):
#     """Extracts the bridle line system data from the surfplan_file
#     This should be done from the surfplan_file, but is instead harcoded
#     """
#     surfplan_discretized_path = (
#         f"{simulation_directory}/src/initialisation/V9_60C/results_bridle_only/"
#     )
#     PULLEY_POINT_INDICES = np.load(
#         f"{surfplan_discretized_path}pulley_point_indices.npy", allow_pickle=True
#     )
#     KCU_POINT_INDICES = np.load(
#         f"{surfplan_discretized_path}kcu_point_indices.npy", allow_pickle=True
#     )
#     KCU_LINE_INDICES = np.load(
#         f"{surfplan_discretized_path}kcu_line_indices.npy", allow_pickle=True
#     )
#     KCU_PLATE_INDICES = np.load(
#         f"{surfplan_discretized_path}kcu_plate_indices.npy", allow_pickle=True
#     )
#     BRIDLE_POINT_INDEX = np.load(
#         f"{surfplan_discretized_path}bridlepoint_index.npy", allow_pickle=True
#     )

#     ##TODO: this should be adjusted from V3 to V9 data
#     BRIDLE_DATA = {
#         "DIAMETER": float(0.002),  # m
#         "RHO": float(230),
#         "BRIDLE_POINT_INDEX": BRIDLE_POINT_INDEX,
#     }  # kg/m3
#     PULLEY_DATA = {
#         "POINT_INDICES": PULLEY_POINT_INDICES,  # [-]
#         "MASS": float(0.1),
#         "NUMBER_OF_PULLEYS_IN_BACK_LINES": int(2),
#     }  # kg
#     KCU_DATA = {
#         "CD": float(0.47),  # [-]
#         "DIAMETER": float(0.38),  # m
#         "MASS": float(8.4),  # kg
#         "INDEX": int(BRIDLE_POINT_INDEX),
#         "POINT_INDICES": KCU_POINT_INDICES,  # [-]
#         "LINE_INDEX": KCU_LINE_INDICES,
#         "KCU_PLATE_INDICES": KCU_PLATE_INDICES,
#     }

#     return BRIDLE_DATA, PULLEY_DATA, KCU_DATA
