import numpy as np


## defining a distance to line function
def distance_point_to_line(point, line):
    """calculates distance between point and line
    return distance, vector and angle"""
    # defining p1, p2, p3
    p1, p2, p3 = line[0], line[1], point
    # defining the vector from p1 to p3
    p1_p3 = p3 - p1
    p1_p3_norm = np.linalg.norm(p1_p3)
    # defining the line vector from p1 to p2
    p1_p2 = p2 - p1
    p1_p2_norm = np.linalg.norm(p1_p2)
    p1_p2_unit = p1_p2 / p1_p2_norm

    # calculating the dot product and realize its equal to p1 to intersection point ps
    # dot(A,B) = |A|*|B|*cos(theta)
    # |B|*cos(theta) = dot(A,B) / |A|
    # take A as the line vector (p1_p2) and B as the diagonal (p1_p3)
    # |p1_p3|*cos(theta) = dot(p1_p2,p1_p3) / |p1_p2|

    # if we think trigonometry
    # cos(theta) = adjacent / diagonal
    # diagonal*cos(theta) = adjacent
    # |p1_p3| *cos(theta) = p1_ps_norm

    # ---> thus p1_ps_norm = dot(p1_p2,p1_p3) / |p1_p2|

    p1_ps_norm = np.dot(p1_p2, p1_p3) / p1_p2_norm
    # find the intersection point, by knowing it lies on the line
    ps = p1 + p1_p2_unit * p1_ps_norm
    # calculate the distance vector
    p3_ps = ps - p3
    # calculate the distance
    p3_ps_norm = np.linalg.norm(p3_ps)
    # calculate the angle between p1_p3 and p1_ps
    cos_angle_p1_p3_and_p1_ps = np.dot(p1_p3, p1_p2) / (p1_p2_norm * p1_p3_norm)
    # avoiding arccos input errors, it only takes -1< x <1
    if cos_angle_p1_p3_and_p1_ps > (1 - 1e-10):
        cos_angle_p1_p3_and_p1_ps = 1 - 1e-10
    elif cos_angle_p1_p3_and_p1_ps < (-1 + 1e-10):
        cos_angle_p1_p3_and_p1_ps = -1 + 1e-10
    angle_p1_p3_and_p1_ps = np.arccos(cos_angle_p1_p3_and_p1_ps)

    distance, vector, angle = p3_ps_norm, p3_ps, angle_p1_p3_and_p1_ps
    return distance, vector, angle


def calculate_force_spring(
    points, bridle_rest_lengths, wing_rest_lengths, input_calculate_force_spring
):
    ## initialise
    bridle_ci = input_calculate_force_spring.bridle_ci
    bridle_cj = input_calculate_force_spring.bridle_cj
    wing_ci = input_calculate_force_spring.wing_ci
    wing_cj = input_calculate_force_spring.wing_cj
    te_line_indices = input_calculate_force_spring.te_line_indices
    tube_line_indices = input_calculate_force_spring.tube_line_indices

    bridle_stiffness = input_calculate_force_spring.bridle_stiffness
    tube_stiffness = input_calculate_force_spring.tube_stiffness
    te_stiffness = input_calculate_force_spring.te_stiffness
    canopy_stiffness = input_calculate_force_spring.canopy_stiffness
    pulley_line_indices = input_calculate_force_spring.pulley_line_indices
    pulley_line_pair_indices = input_calculate_force_spring.pulley_line_pair_indices

    is_with_elongation_limit = input_calculate_force_spring.is_with_elongation_limit
    elongation_limit = input_calculate_force_spring.elongation_limit
    is_with_compression_limit = input_calculate_force_spring.is_with_compression_limit
    compression_limit = input_calculate_force_spring.compression_limit
    limit_stiffness_factor = input_calculate_force_spring.limit_stiffness_factor

    def calculate_fs(rest_length, K, node_i, node_j, is_only_elongation_resistance):
        """Calculates the spring force between two nodes
        input:  rest_length, K, node_i, node_j, is_only_elongation_resistance
        output: F_spring, unit_vector"""

        ## Calculate the spring force
        vector = (
            node_i - node_j
        )  # Vector(ci --> cj) separating the points, indicating direction
        length = np.linalg.norm(
            vector
        )  # Absolute magnitude of the vector (works both ways), indicating strength
        elongation = (
            length - rest_length
        ) / rest_length  # SpringL is defined on a range(0,len(ci)) loop (works both ways)
        unit_vector = vector / length  # Define the unit_vector (ci --> cj)
        F_spring = K * elongation

        if is_only_elongation_resistance:  # if the spring only works in elongation
            if elongation < 0:
                F_spring = 0  # then: if compressed, set F_spring to zero

            # if elongated, with elongation limit and above the elongation limit
            elif (
                is_with_elongation_limit and elongation > elongation_limit
            ):  # if elongation is above the set limit
                F_spring *= limit_stiffness_factor

        # if the spring works in both directions
        else:
            # if with compression limit, compressed and below the compression limit
            if is_with_compression_limit and elongation < compression_limit:
                F_spring *= limit_stiffness_factor

        return F_spring, unit_vector

    # Define a spring force matrix of the right size
    spring_force = np.zeros(
        points.shape
    )  # Initialising with zero matrix in same shape as points

    ## Bridle lines
    for idx, (idx_bridle_node_i, idx_bridle_node_j) in enumerate(
        zip(bridle_ci, bridle_cj)
    ):  # loop through each bridle line
        ## set parameters to bridle
        K = bridle_stiffness
        is_only_elongation_resistance = True
        springL = bridle_rest_lengths[idx]
        node_i = points[idx_bridle_node_i]
        node_j = points[idx_bridle_node_j]

        ## If the bridle line is a pulley-line
        # pulley_line_pair_indices is a dict, with each key being index_1 of pulley line and the value = index_2
        if str(idx) in pulley_line_pair_indices:
            ## Dealing with the first bit of line --> idx_1
            idx_1 = idx  # Vector(ci --> cj) separating the points, indicating direction
            vector_1 = points[bridle_ci[idx_1]] - points[bridle_cj[idx_1]]
            length_1 = np.linalg.norm(vector_1)
            rest_length_1 = bridle_rest_lengths[idx]
            unit_vector_1 = vector_1 / length_1  # Define the unit_vector (ci --> cj)

            ## Dealing with the second bit of line --> idx_2
            idx_2 = pulley_line_pair_indices[str(idx)]
            vector_2 = (
                points[bridle_ci[idx_2]] - points[bridle_cj[idx_2]]
            )  # Vector(ci --> cj) separating the points, indicating direction
            length_2 = np.linalg.norm(vector_2)
            rest_length_2 = bridle_rest_lengths[idx_2]
            unit_vector_2 = vector_2 / length_2  # Define the unit_vector (ci --> cj)

            ## Finding the combined elongation & spring-force
            elongation = ((length_1 + length_2) - (rest_length_1 + rest_length_2)) / (
                rest_length_1 + rest_length_2
            )
            F_spring = K * elongation

            ## Check for compression
            if (
                not is_only_elongation_resistance
            ):  # if the spring only works in elongation
                if elongation < 0:
                    F_spring = 0  # then: if compressed, set F_spring to zero

            ## Apply spring force
            # to first node pair
            spring_force[bridle_ci[idx_1]] += -F_spring * np.array(
                [unit_vector_1[0], unit_vector_1[1], unit_vector_1[2]]
            )  # to node_i
            spring_force[bridle_cj[idx_1]] += F_spring * np.array(
                [unit_vector_1[0], unit_vector_1[1], unit_vector_1[2]]
            )  # OPpointsITE DIRECTION to node-j
            # to second node pair
            spring_force[bridle_ci[idx_2]] += -F_spring * np.array(
                [unit_vector_2[0], unit_vector_2[1], unit_vector_2[2]]
            )  # to node_i
            spring_force[bridle_cj[idx_2]] += F_spring * np.array(
                [unit_vector_2[0], unit_vector_2[1], unit_vector_2[2]]
            )  # OPpointsITE DIRECTION to node-j

        ## If instead we are dealing with ordinary bridle lines between knots
        elif idx not in pulley_line_indices:  # IF not a pulley line
            ## Finding the spring-force: with only elongation resistance (i.e. compression resistance is zero)
            F_spring, unit_vector = calculate_fs(
                springL, K, node_i, node_j, is_only_elongation_resistance
            )

            ## Apply spring force
            spring_force[idx_bridle_node_i] += -F_spring * np.array(
                [unit_vector[0], unit_vector[1], unit_vector[2]]
            )  # to node_i
            spring_force[idx_bridle_node_j] += F_spring * np.array(
                [unit_vector[0], unit_vector[1], unit_vector[2]]
            )  # OPpointsITE DIRECTION to node-j

    ## Wing lines
    for idx, (idx_wing_node_i, idx_wing_node_j) in enumerate(
        zip(wing_ci, wing_cj)
    ):  # loop through each bridle line
        ## Set parameters to wing
        springL = wing_rest_lengths[idx]
        node_i = points[idx_wing_node_i]
        node_j = points[idx_wing_node_j]

        ## If the line is a tube, spring force works both ways (compression & elongation)
        if idx in tube_line_indices:
            K = tube_stiffness
            is_only_elongation_resistance = (
                False  # spring force works in both directions
            )

        ## If the line is a TE line, spring force only works in elongation
        elif idx in te_line_indices:
            K = te_stiffness
            is_only_elongation_resistance = (
                True  # spring force only works in elongation
            )

        ## If the line is a diagonal canopy line, spring force only works in elongation
        elif idx not in tube_line_indices and idx not in te_line_indices:
            K = canopy_stiffness
            is_only_elongation_resistance = (
                True  # spring force only works in elongation
            )

        ## Finding the spring-force: with only elongation resistance (i.e. compression resistance is zero)
        F_spring, unit_vector = calculate_fs(
            springL, K, node_i, node_j, is_only_elongation_resistance
        )

        ## Apply spring force
        spring_force[idx_wing_node_i] += -F_spring * np.array(
            [unit_vector[0], unit_vector[1], unit_vector[2]]
        )  # to node_i
        spring_force[idx_wing_node_j] += F_spring * np.array(
            [unit_vector[0], unit_vector[1], unit_vector[2]]
        )  # OPpointsITE DIRECTION to node-j

    return spring_force


# def calculate_elongation(points,bridle_rest_lengths,wing_rest_lengths,bridle_ci,bridle_cj,wing_ci,wing_cj,te_line_indices,tube_line_indices,STIFFNESS_DATA):

#     BRIDLE_STIFFNESS = STIFFNESS_DATA['BRIDLE']
#     tube_stiffness = STIFFNESS_DATA['TUBE']
#     te_stiffness = STIFFNESS_DATA['TE']
#     canopy_stiffness = STIFFNESS_DATA['CANOPY']

#     def calculate_fs(rest_length,K,node_i,node_j,is_only_elongation_resistance = True):
#         ''' Calculates the spring force between two nodes
#             input:  rest_length, K, node_i, node_j, is_only_elongation_resistance
#             output: F_spring, unit_vector'''

#         ## Calculate the spring force
#         vector = node_i - node_j                                # Vector(ci --> cj) separating the points, indicating direction
#         length = np.linalg.norm(vector)                 # Absolute magnitude of the vector (works both ways), indicating strength
#         elongation = (length - rest_length)/rest_length # SpringL is defined on a range(0,len(ci)) loop (works both ways)
#         unit_vector = vector / length           # Define the unit_vector (ci --> cj)
#         F_spring  = K*elongation

#         if is_only_elongation_resistance:      # if the spring only works in elongation
#             if elongation < 0: F_spring  = 0    # then: if compressed, set F_spring to zero

#         return F_spring,unit_vector,elongation

#     # Define a spring force matrix of the right size
#     spring_force = np.zeros(points.shape) #Initialising with zero matrix in same shape as points

#     bridle_elongation_values,wing_elongation_values = np.zeros(len(bridle_ci)),np.zeros(len(wing_ci))

#     ## Bridle lines
#     for idx, (idx_bridle_node_i,idx_bridle_node_j) in enumerate(zip(bridle_ci,bridle_cj)): # loop through each bridle line

#         ## set parameters to bridle
#         K       = BRIDLE_STIFFNESS
#         springL = bridle_rest_lengths[idx]
#         node_i  = points[idx_bridle_node_i]
#         node_j  = points[idx_bridle_node_j]

#         ## If the bridle line is a pulley-line
#         #if (idx_bridle_node_i and idx_bridle_node_j) in pulley_indices: ##TODO: this would be better than hardcoding
#         if idx == 4 or idx == 5:

#             ## Dealing with the first bit of line --> idx_1
#             idx_1 = idx
#             vector_1        = points[bridle_ci[idx_1]] - points[bridle_cj[idx_1]]                   # Vector(ci --> cj) separating the points, indicating direction
#             length_1        = np.linalg.norm(vector_1)
#             rest_length_1   = bridle_rest_lengths[idx_1]
#             unit_vector_1   = vector_1 / length_1              # Define the unit_vector (ci --> cj)

#             ## Dealing with the second bit of line --> idx_2
#             idx_2 = idx_1*2 -2 #4 --> 6, 5 --> 8
#             vector_2        = points[bridle_ci[idx_2]] - points[bridle_cj[idx_2]]    # Vector(ci --> cj) separating the points, indicating direction
#             length_2        = np.linalg.norm(vector_2)
#             rest_length_2    = bridle_rest_lengths[idx_2]
#             unit_vector_2   = vector_2 / length_2              # Define the unit_vector (ci --> cj)

#             ## Finding the combined elongation & spring-force
#             elongation = ((length_1+ length_2) - (rest_length_1+rest_length_2)) /(rest_length_1 + rest_length_2)
#             F_spring  = K*elongation

#             ## Check for compression
#             if elongation < 0: F_spring  = 0    #if compressed, set F_spring to zero

#             ## Apply spring force
#             # to first node pair
#             spring_force[bridle_ci[idx_1]] += -F_spring * np.array([unit_vector_1[0],unit_vector_1[1],unit_vector_1[2]]) # to node_i
#             spring_force[bridle_cj[idx_1]] +=  F_spring * np.array([unit_vector_1[0],unit_vector_1[1],unit_vector_1[2]]) # OPpointsITE DIRECTION to node-j
#             # to second node pair
#             spring_force[bridle_ci[idx_2]] += -F_spring * np.array([unit_vector_2[0],unit_vector_2[1],unit_vector_2[2]]) # to node_i
#             spring_force[bridle_cj[idx_2]] +=  F_spring * np.array([unit_vector_2[0],unit_vector_2[1],unit_vector_2[2]]) # OPpointsITE DIRECTION to node-j


#         ## If instead we are dealing with ordinary bridle lines between knots
#         elif idx != 6 and idx != 8: # making sure they are not pulley points ##TODO: remove this hardcoding

#             ## Finding the spring-force: with only elongation resistance (i.e. compression resistance is zero)
#             F_spring,unit_vector,elongation  = calculate_fs(springL,K,node_i,node_j,is_only_elongation_resistance = False)

#             ## Apply spring force
#             spring_force[idx_bridle_node_i] += -F_spring * np.array([unit_vector[0],unit_vector[1],unit_vector[2]]) # to node_i
#             spring_force[idx_bridle_node_j] +=  F_spring * np.array([unit_vector[0],unit_vector[1],unit_vector[2]]) # OPpointsITE DIRECTION to node-j

#         bridle_elongation_values[idx] = elongation

#     ## Wing lines
#     for idx, (idx_wing_node_i,idx_wing_node_j) in enumerate(zip(wing_ci,wing_cj)): # loop through each bridle line

#         ## Set parameters to wing
#         springL = wing_rest_lengths[idx]
#         node_i  = points[idx_wing_node_i]
#         node_j  = points[idx_wing_node_j]

#         ## If the line is a tube, spring force works both ways (compression & elongation)
#         if idx in tube_line_indices:
#             K = tube_stiffness
#             is_only_elongation_resistance = True # spring force works in both directions

#         ## If the line is a TE line, spring force only works in elongation
#         elif idx in te_line_indices:
#             K = te_stiffness
#             is_only_elongation_resistance = False # spring force only works in elongation

#         ## If the line is a diagonal canopy line, spring force only works in elongation
#         elif idx not in tube_line_indices and idx not in te_line_indices:
#             K = canopy_stiffness
#             is_only_elongation_resistance = False # spring force only works in elongation

#         ## Finding the spring-force: with only elongation resistance (i.e. compression resistance is zero)
#         F_spring,unit_vector,elongation  = calculate_fs(springL,K,node_i,node_j,is_only_elongation_resistance)

#         ## Apply spring force
#         spring_force[idx_wing_node_i] += -F_spring * np.array([unit_vector[0],unit_vector[1],unit_vector[2]]) # to node_i
#         spring_force[idx_wing_node_j] +=  F_spring * np.array([unit_vector[0],unit_vector[1],unit_vector[2]]) # OPpointsITE DIRECTION to node-j

#         wing_elongation_values[idx] = elongation

#     return bridle_elongation_values,wing_elongation_values


def get_symmetrical(points, tolerance=1e-5):
    # loop through each point
    for i, point_left in enumerate(points):
        if point_left[1] > 0:  # if the point is on the left-side
            # loop through each other point
            for j, point_right in enumerate(points):
                # if the point is on the right-side and the mirror point of i
                # mirror it, by making it equal to the left-point, but negative in the y
                if (
                    point_right[1] < 0
                    and i != j
                    and np.isclose(point_left[0], point_right[0], atol=tolerance)
                    and np.isclose(point_left[2], point_right[2], atol=tolerance)
                    and np.isclose(point_left[1], -point_right[1], atol=tolerance)
                ):
                    points[j] = [point_left[0], -point_left[1], point_left[2]]

    return points


# %% OLD FUNCTIONS

# def get_spring_force_kite_no_damp(ci,cj,points,springL,K_kite,canopy_stiffness,spring_force,tube_idx):

# 	for i in range(0,len(ci)): #i loops over each (bridle or kite) line

# 		sep_vec= points[ci[i]] - points[cj[i]] 	# Vector(ci --> cj) separating the points, indicating direction
# 		sep = np.linalg.norm(sep_vec)		# Absolute magnitude of the vector (works both ways), indicating strength
# 		unit_vector = sep_vec/sep 			# Define the unit_vector (ci --> cj)
# 		dL = (sep - springL[i])/1e3 			# SpringL is defined on a range(0,len(ci)) loop (works both ways)

# 		K_factor = K_kite
# 		TE = [2,8,14,20,26,32,38,44,50]
# 		sec = 0
# 		if i not in tube_idx and i not in TE:
# 			if i%3 == 0:
# 				sec +=1
# 			K_factor = canopy_stiffness
# 			if dL<0:
# 				K_factor = 0
# # 		if i in TE:
# #  			K_factor = canopy_stiffness
# #  			if dL<0: K_factor = 0
# 		#Apply spring force to ci[i]
# 		spring_force[ci[i], 0] += -K_factor * dL * unit_vector[0]  # fill x_coord of point in spring_force with: K*dL*unit_vectir
# 		spring_force[ci[i], 1] += -K_factor * dL * unit_vector[1]
# 		spring_force[ci[i], 2] += -K_factor * dL * unit_vector[2]

# 		#Apply spring force in oppointsite direction to cj[i]
# 		spring_force[cj[i], 0] += -K_factor * dL * -unit_vector[0]   # fill x_coord of point in spring_force with: K*dL*unit_vectir
# 		spring_force[cj[i], 1] += -K_factor * dL * -unit_vector[1]
# 		spring_force[cj[i], 2] += -K_factor * dL * -unit_vector[2]


# 	return spring_force


# def get_spring_force_no_damp(ci,cj,points,springL,K,spring_force,tube_idx):

# 	for i in range(0,len(ci)): #i loops over each (bridle or kite) line

# 		if i == 4 or i == 5: #Hardcode the pulley occurences
# 			sep_vec_1 = points[ci[i]] - points[cj[i]]  # Vector(ci --> cj) separating the points, indicating direction
# 			sep_1 = np.linalg.norm(sep_vec_1)  # Absolute magnitude of the vector (works both ways), indicating strength
# 			unit_vector_1 = sep_vec_1 / sep_1  # Define the unit_vector (ci --> cj)

# 			i_n = i*2 -2 #4 --> 6, 5 --> 8
# 			sep_vec_2 = points[ci[i_n]] - points[cj[i_n]]  # Vector(ci --> cj) separating the points, indicating direction
# 			sep_2 = np.linalg.norm(sep_vec_2)  # Absolute magnitude of the vector (works both ways), indicating strength
# 			unit_vector_2 = sep_vec_2 / sep_2  # Define the unit_vector (ci --> cj)


# 			dL =((sep_1+sep_2) - (springL[i]+springL[i_n]))/1e3 # SpringL is defined on a range(0,len(ci)) loop (works both ways)


# 			### Making a smooth decrease in spring_force, for stability
# 			if dL >= 0 : K_factor = K*(springL[i]+springL[i_n])*dL/1e3
# 			elif dL <0	: K_factor = 0  # Think mm's

# 			# Apply spring force to ci[i]
# 			spring_force[ci[i], 0] += -K_factor * unit_vector_1[0]  # fill x_coord of point in spring_force with: K*dL*unit_vectir
# 			spring_force[ci[i], 1] += -K_factor * unit_vector_1[1]
# 			spring_force[ci[i], 2] += -K_factor * unit_vector_1[2]

# 			# Apply spring force in oppointsite direction to cj[i]
# 			spring_force[cj[i], 0] += -K_factor * -unit_vector_1[0]  # fill x_coord of point in spring_force with: K*dL*unit_vectir
# 			spring_force[cj[i], 1] += -K_factor * -unit_vector_1[1]
# 			spring_force[cj[i], 2] += -K_factor * -unit_vector_1[2]

# 			# Apply spring force to ci[i_n]
# 			spring_force[ci[i_n], 0] += -K_factor * unit_vector_2[0]   # fill x_coord of point in spring_force with: K*dL*unit_vectir
# 			spring_force[ci[i_n], 1] += -K_factor * unit_vector_2[1]
# 			spring_force[ci[i_n], 2] += -K_factor * unit_vector_2[2]

# 			# Apply spring force in oppointsite direction to cj[i_n]
# 			spring_force[cj[i_n], 0] += -K_factor * -unit_vector_2[0]  # fill x_coord of point in spring_force with: K*dL*unit_vectir
# 			spring_force[cj[i_n], 1] += -K_factor * -unit_vector_2[1]
# 			spring_force[cj[i_n], 2] += -K_factor * -unit_vector_2[2]

# 		elif i != 6 and i != 8: #Knot bridle points
# 			sep_vec= points[ci[i]] - points[cj[i]] 	# Vector(ci --> cj) separating the points, indicating direction
# 			sep = np.linalg.norm(sep_vec)		# Absolute magnitude of the vector (works both ways), indicating strength
# 			unit_vector = sep_vec/sep 			# Define the unit_vector (ci --> cj)
# 			dL = (sep - springL[i])/1e3  			# SpringL is defined on a range(0,len(ci)) loop (works both ways)

# 			### Making a smooth decrease in spring_force, for stability
# 			if dL >= 0 : K_factor = K*springL[i]*dL/1e3
# 			elif dL <0	: K_factor = 0  # Think mm's

# 			#Apply spring force to ci[i]
#             # fill x_coord of point in spring_force with: K*dL*unit_vectir
# 			spring_force[ci[i], 0] += -K_factor * unit_vector[0]
# 			spring_force[ci[i], 1] += -K_factor * unit_vector[1]
# 			spring_force[ci[i], 2] += -K_factor * unit_vector[2]

# 			#Apply spring force in oppointsite direction to cj[i]
# 			spring_force[cj[i], 0] += -K_factor * -unit_vector[0]  # fill x_coord of point in spring_force with: K*dL*unit_vectir
# 			spring_force[cj[i], 1] += -K_factor * -unit_vector[1]
# 			spring_force[cj[i], 2] += -K_factor * -unit_vector[2]

# 	return spring_force

# def get_symmetrical_V3(points):

# 	points[22] = [points[22][0],0,points[22][2]]
# 	left_kite  = [5,4,3,2,19,15,16,17,18,1 ]
# 	right_kite = [6,7,8,9,20,14,13,12,11,10]

# 	left_bridle_LE  = [24,31,32,33]
# 	right_bridle_LE = [28,34,35,36]

# 	left_bridle_TE = [23,24,25,26]
# 	right_bridle_TE= [27,28,29,30]

# 	def mirroring(points,left,right):
# 		for i in range(0,len(left)):
# 			points[right[i]] = [points[left[i]][0],-points[left[i]][1],points[left[i]][2]]
# 		return points

# 	points = mirroring(points,left_kite,right_kite)
# 	points = mirroring(points,left_bridle_LE,right_bridle_LE)
# 	points = mirroring(points,left_bridle_TE,right_bridle_TE)

# 	return points
