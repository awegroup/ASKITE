#%% Functions

#importing libraries
import re
import numpy as np
import math
from scipy.interpolate import interp1d # Functions necessary for finding each point

def read_obj_file_wing(file_path):
    ''' Reads the .obj file and returns the vertices, normals, and faces'''
    
    # libraries
    import re

    vertices = []
    normals = []
    faces = []

                                
    vertex_pattern = re.compile(r'v\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)')  #  v 2669.676921 661.911471 25.857700
    normal_pattern = re.compile(r'vn\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)')
    face_pattern = re.compile(r'f\s + (-?\d+)/(-?\d+)/(-?\d+)\s + (-?\d+)/(-?\d+)/(-?\d+)\s + (-?\d+)/(-?\d+)/(-?\d+)\s ') # f 18003/11004/21081 18009/11011/21088 18250/11012/21089 18244/11005/21082
    f_pattern = re.compile(r'f\s+.')

    flag = True
    with open(file_path, 'r') as f:
        #for line in range(len(f)):
        for index, line in enumerate(f):

            vertex_match = vertex_pattern.match(line)
            if vertex_match and flag == True:
                # This line describes a vertex
                x, y, z = map(float, vertex_match.groups())
                vertices.append((x, y, z))

            elif flag == True:
                normal_match = normal_pattern.match(line)
                if normal_match:
                    # This line describes a vertex normal
                    x, y, z = map(float, normal_match.groups())
                    normals.append((x, y, z))

                elif f_pattern.match(line): #when not vertex or normal
                    faces.append([num.split('/')[0] for num in line.split()[1:]])

            # Setting the flag = False once the bridle line pattern is found        
            # This is done by checking if the word bridle is present in the splitted line
            
            if 'bridle' in re.split(r"-|,|_|\n", line):
                flag = False
                break

    return vertices, normals, faces

def get_planar_in_order(point_list, tol_distance_off_planar, tol_min_points_for_planar, tol_max_points_for_planar):
    '''
    This function takes a list of points and returns a list of planar shapes.
    A planar shape is a list of points that lie within a certain tolerance of a plane.

    input:  point_list, 
            tol_distance_off_planar       ~ 1 [mm] 
            min_points_TubeCircle_Pulley  ~ 30 (close to number of points in a circle)
    output: min_points_TubeCircle_Pulley '''
    
    
    point_list = np.array(point_list)

    # Initialize a list to store the points that lie within the tolerance of the plane
    planar_shape_list,planar_shape = [],[point_list[0],point_list[1],point_list[2]]
    distance_list = []
    planar_shape = [] #initialize in case the first 3 points are not planar
    not_planar_shapes = []
    # Iterate through the elements of the list
    for i in range(3,len(point_list)): # skip the first 2 points as you need 3 to define a surface

        # Calculate the vectors formed by the last 2 points
        vec_1 = point_list[i-2] - point_list[i-3]
        vec_2 = point_list[i-1] - point_list[i-3]

        # Calculate the cross product of the vectors to get the normal vector of the plane
        cross = np.cross(vec_1, vec_2)

        # Calculate the distance from the plane
        distance = np.abs(np.dot(point_list[i] - point_list[i-2], cross)) / np.linalg.norm(cross)
        distance_list.append(distance)

        # Check if the distance is within the tolerance
        if distance <= tol_distance_off_planar:
            planar_shape.append(point_list[i]) # append the point to the circle
        else: # if not planar, within the defined tolerance
            #not_planar_shapes.append(point_list[i])
            if len(planar_shape) > tol_min_points_for_planar and len(planar_shape) < tol_max_points_for_planar:
                planar_shape_list.append(planar_shape) # append the circle to the circle list
                planar_shape = [] # reset the circle
            elif len(planar_shape)==0:
                not_planar_shapes.append(point_list[i])
                planar_shape = []
            else: #if planar_shape to short or to long
                for point in planar_shape:
                    not_planar_shapes.append(point)
                planar_shape = []
    return planar_shape_list,not_planar_shapes



def get_circle_shapes(planar_shape_list, tol_diameter, min_counter_for_circle):
    '''
    This function takes a list of planar shapes and returns a list of circles and a list of radii

    input:  planar_shape_list 
            tol_diameter            ~ 1 [mm]
            min_counter_for_circle  ~ 15 (close to number of points in a circle)
                (number of points that must lie within the diameter)

    output: circle_list
            radius_list
            center_list    
    '''

    def find_max_distance(points, otherPoint):
        ''' This function takes a list of points and calculates the maximum distance
            between the point and any other point in the list of points.

            input:      points = list of points
                        otherPoint = point that you want to find the maximum distance to

            output:     maxDistance
                        index_point_furthest_away = index of the point that is furthest away from the otherPoint
        '''
        ## TODO: can probably also do this more neatly with scipy distance matix

        maxDistance = 0 #initialize to 0
        index_point_furthest_away = 0 
        # Iterate over the list of points
        for i in range(len(points)):
            point = points[i]
            # Calculate the distance between the point and the other point
            distance = math.sqrt((point[0] - otherPoint[0])**2 + (point[1] - otherPoint[1])**2 + (point[2] - otherPoint[2])**2)
            # If the distance is greater than the maximum distance, update the farthest point and the maximum distance
            if distance > maxDistance:
                maxDistance = distance
                index_point_furthest_away = i #each time the max distance is updated, the index of the point is also updated

        return maxDistance, index_point_furthest_away

    def findMedian(numbers): # defining the median extractor function
        ''' This function takes a list of numbers and returns the median of the list of numbers

            input: list of numbers
            output: median of the list of numbers'''
        # Sort the list of numbers
        numbers.sort()

        # If the length of the list is odd, return the middle number
        if len(numbers) % 2 == 1:
            return numbers[len(numbers) // 2]
        # If the length of the list is even, return the average of the two middle numbers
        else:
            return (numbers[len(numbers) // 2 - 1] + numbers[len(numbers) // 2]) / 2



    circle_list,radius_list,center_list,non_circle_list = [],[],[],[]
    for planar_shape in planar_shape_list:

        ### TODO: check the effect this removing of duplicates has         
        ## making sure the lists is unique
        # Use a list comprehension and the allclose method to remove duplicates
        #points_list = planar_shape
        #planar_shape = [point for i, point in enumerate(points_list) if not any(np.allclose(point, other_point) for other_point in points_list[:i])]


        # find a list with farthest distances for each point to any other point in the planar shape
        farthest_distance_list,mid_point_list = [],[]
        for point in planar_shape:
            farthest_distance,index_point_furthest_away = find_max_distance(planar_shape, point) # find the farthest distance
            farthest_distance_list.append(farthest_distance) # append the farthest distance to the list

            # find the mid-point between the two points of max_distance
            point1 = point
            point2 = planar_shape[index_point_furthest_away]
            mid_point = ((point1[0] + point2[0]) / 2, (point1[1] + point2[1]) / 2, (point1[2] + point2[2]) / 2)
            mid_point_list.append(mid_point)

        # count the number of points that match the farthest_distance to the median farthest distance
        centers,radii =[],[]
        counter = 0
        median_distance = findMedian(farthest_distance_list)
        for i in range(len(farthest_distance_list)): # iterate through the list of farthest distances
            distance = farthest_distance_list[i]
            if abs(distance-median_distance) < tol_diameter: # if the distance close to the median distance
                counter = counter + 1 # increase the counter
                centers.append(mid_point_list[i]) # append the mid_point to the list of centers
                radii.append(distance/2)


        # if enough max_distances (read "diameter") match, it must be a circle
        if counter > min_counter_for_circle: # if the counter is larger than the minimum number of points that must lie within the diameter
            circle_list.append(planar_shape) 
            radius_list.append(np.mean(radii))
            # finding the best approximation of the center
            center_x_mean = np.mean([center[0] for center in centers])
            center_y_mean = np.mean([center[1] for center in centers])
            center_z_mean = np.mean([center[2] for center in centers])
            center_list.append([center_x_mean,center_y_mean,center_z_mean])
        
        else: #if not deemed a circle
            non_circle_list.append(planar_shape)

    return circle_list,radius_list,center_list,non_circle_list

def get_strut_radius(circles,circles_radii):
    '''
    Finding the strut diameter, by checking which circle has the smallest z-value

    input:  circle_list
            radius_list

    output: strut_radius
    '''
    min_z_list =[]
    for i in range(len(circles)):
        min_z_list.append(min(z[2] for z in circles[i]))

    minIndex = min_z_list.index(min(min_z_list))
    strut_radius = circles_radii[minIndex]
    return strut_radius

def split_LE_strut_ends_wing(almost_tubular_frame, almost_tubular_frame_radii,almost_tubular_frame_centroids,LE_radius):
  
  new_tubular_frame = []
  new_tubular_frame_radii = []
  new_tubular_frame_centroids = []

  for i,circular_shape in enumerate(almost_tubular_frame):

    if circular_shape[0][2] > 0 and circular_shape[0][1] > 1600:
      if almost_tubular_frame_radii[i] > LE_radius:
        new_tubular_frame.append(circular_shape)
        new_tubular_frame_radii.append(almost_tubular_frame_radii[i])
        new_tubular_frame_centroids.append(almost_tubular_frame_centroids[i])
      
    else:
      new_tubular_frame.append(circular_shape)
      new_tubular_frame_radii.append(almost_tubular_frame_radii[i])
      new_tubular_frame_centroids.append(almost_tubular_frame_centroids[i])

  return new_tubular_frame, new_tubular_frame_radii, new_tubular_frame_centroids


def split_circles_on_radius(circles,circles_raddi,circles_centroids,split_radius):
    '''
    Splitting lists of circles that are larger than the strut diameter

    input:  circle_list
            radius_list
            center_list
            split_radius ~ 0.8*strut_radius (just make it slightly smaller)

    output: big_circles
            big_circles_centroids
            small_circles
            small_circles_centroids
    '''
    big_circles_radii = []
    # Based on the strut diameter splitting the tubular frame from the pulleys and tube-valves
    big_circles,big_circles_centroids,small_circles,small_circles_centroids = [],[],[],[]
    for i in range(len(circles)):
        if circles_raddi[i] > (split_radius):
            big_circles.append(circles[i])
            big_circles_centroids.append(circles_centroids[i])
            big_circles_radii.append(circles_raddi[i])
        else:
            small_circles.append(circles[i])
            small_circles_centroids.append(circles_centroids[i])

    return big_circles,big_circles_centroids,big_circles_radii,small_circles,small_circles_centroids


def split_LE_strut_ends(almost_tubular_frame, almost_tubular_frame_radii,almost_tubular_frame_centroids,LE_radius):
  
    new_tubular_frame = []
    new_tubular_frame_radii = []
    new_tubular_frame_centroids = []

    for i,circular_shape in enumerate(almost_tubular_frame):

        if circular_shape[0][2] > 0 and circular_shape[0][1] > 1600:
            if almost_tubular_frame_radii[i] > LE_radius:
                new_tubular_frame.append(circular_shape)
                new_tubular_frame_radii.append(almost_tubular_frame_radii[i])
                new_tubular_frame_centroids.append(almost_tubular_frame_centroids[i])
            
        else:
            new_tubular_frame.append(circular_shape)
            new_tubular_frame_radii.append(almost_tubular_frame_radii[i])
            new_tubular_frame_centroids.append(almost_tubular_frame_centroids[i])

        return new_tubular_frame, new_tubular_frame_radii, new_tubular_frame_centroids

def split_plate_corners_from_tubular_frame(tubular_frame_centroids,tol_plate_corners_x,tol_plate_corners_y,min_z_ratio):
    '''
    This function splits the corners of the plate from the rest of the tubular frame
    
    input:  tubular_frame_centroids: list of centroids of the tubular frame
            tol_plate_corners_x: ~150 mm (this tolerance is used to indicate that two points could form a strut)
            tol_plate_corners_y: ~100 mm (this tolerance is used to indicate that two points could form a strut)
            min_z_ratio:         ~0.5 (this ratio is used to separate the strut ends (and tips) from the LE)

    output: plate_corners: list of centroids of the plate corners
            tubular_frame_centroids: list of centroids of the tubular frame without the plate corners
            N_corners: number of corners of the plate
            N_struts:  number of struts of the plate
            N_panels:  number of panels of the plate
    '''

    min_z = min(z[2] for z in tubular_frame_centroids)
    plate_corners,index_list_plate_corners= [],[]
    for i in range(len(tubular_frame_centroids)):
        # if the z-point is close to the bottom of the frame, i.e. either tip or strut end
        if tubular_frame_centroids[i][2]< min_z*min_z_ratio: 
            for j in range(len(tubular_frame_centroids)):
                if i!=j:

                    if tubular_frame_centroids[i][0] <=0:
                        x_diff = abs(tubular_frame_centroids[i][0] - tubular_frame_centroids[j][0])
                        y_diff = abs(tubular_frame_centroids[i][1] - tubular_frame_centroids[j][1])

                        if x_diff < tol_plate_corners_x and y_diff < tol_plate_corners_y:
                            plate_corners.append(tubular_frame_centroids[i])
                            plate_corners.append(tubular_frame_centroids[j])
                            index_list_plate_corners.append(i)
                            index_list_plate_corners.append(j)
                            break #for each TE (or tip) point, there should be only one LE point
    
    length_plate_corners = len(plate_corners)
    for i in range(length_plate_corners):
        plate_corners.append([-plate_corners[i][0],plate_corners[i][1],plate_corners[i][2]])

    # Defining some numbers
    N_corners    = int(len(plate_corners))
    N_struts     = int(N_corners/2)
    N_panels     = int(N_struts -1)

    # Defining those points, that are not plate corners
    tubular_frame_centroids_non_plate_corners = []
    for i in range(len(tubular_frame_centroids)):
        if i not in index_list_plate_corners:
            tubular_frame_centroids_non_plate_corners.append(tubular_frame_centroids[i])

    return plate_corners,tubular_frame_centroids_non_plate_corners,N_corners,N_struts,N_panels

def transform_coordinate_system_Surfplan_to_Model(point_list):
    #####################################################################
    # Surfplan defines 
        # x positive to the right 
        # y positive upwards
        # z positive against the wind, from TE to LE
    #####################################################################

    #####################################################################
    # Model defines 
        # x positive with the wind, from LE to TE 
        # y positive to the left
        # z positive upwards
    #####################################################################

    point_list_new = []
    for point in point_list:
        x_new = -point[2]
        y_new = -point[0]
        z_new = point[1]
        point_list_new.append([x_new,y_new,z_new])
    return point_list_new

def get_kite_plate_connectivity(points_ini,N_corners):
	

    # sorting the plate corners on the y-coordinate from left-to-right (high to low)
    sorted_plate_corners = sorted(points_ini[0:N_corners], key=lambda point: point[1], reverse=True)

    # filling the plates list in the appropriate order left-right looking from below 
    plates_points = []
    N_panels = int(N_corners/2-1)
    for i in range(N_panels): #loop through each panel
        strut1 = [sorted_plate_corners[i*2]  ,sorted_plate_corners[i*2+1]] #strut1 is first strut from left-to-right of the panel
        strut2 = [sorted_plate_corners[i*2+2],sorted_plate_corners[i*2+3]] #strut2 is second strut from left-to-right of the panel
        sorted_strut1 = sorted(strut1, key=lambda point: point[0]) #sorting to have the LE point first
        sorted_strut2 = sorted(strut2, key=lambda point: point[0]) #sorting to have the LE point first
        #plates are defined clock-wise as: [LE1,LE2,TE2,TE1]
        plates_points.append([sorted_strut1[0],sorted_strut2[0],sorted_strut2[1],sorted_strut1[1]])

    # Getting from points to indices
    plates = []
    for plate in plates_points:
        plate_index = []
        for point in plate:
            index = points_ini.index(point)
            plate_index.append(index)
        plates.append(plate_index)

    # Using the structure of the plates to get the right connectivity matrix structure
    ci_kite, cj_kite = [], []
    for i in np.arange(0, len(plates)):
        # The 4 lines describing the tubular frame
        ci_kite.append(plates[i][0]) #LE
        cj_kite.append(plates[i][1]) #LE

        ci_kite.append(plates[i][1]) #Strut right
        cj_kite.append(plates[i][2]) #Strut right

        ci_kite.append(plates[i][2]) #TE
        cj_kite.append(plates[i][3]) #TE

        ci_kite.append(plates[i][3]) #Strut left
        cj_kite.append(plates[i][0]) #Strut left

        # Them diagonals
        ci_kite.append(plates[i][0])
        cj_kite.append(plates[i][2])

        ci_kite.append(plates[i][1])
        cj_kite.append(plates[i][3])

    ci_kite = np.reshape(ci_kite,len(ci_kite))
    cj_kite = np.reshape(cj_kite, len(cj_kite))
    
    return ci_kite, cj_kite, plates


def struct2aero_geometry(points_ini,N_corners):


    points_ini = np.array(points_ini)

    # sorting the plate corners on the y-coordinate (from right to left)
    sorted_plate_corners = sorted(points_ini[1:N_corners+1], key=lambda point: point[1])
    N_struts = int(len(sorted_plate_corners)/2)
    
    # filling the coord list in the appropriate order right-left looking from below 
    coord = np.empty((N_corners,3))
    
    for i in range(N_struts): #loop through each strut
        strut = [sorted_plate_corners[i*2]  ,sorted_plate_corners[i*2+1]] #strut1 = right tip, strut2 = 1st right strut, etc.
        sorted_strut = sorted(strut, key=lambda point: point[0]) #sorting to have the LE point first
        #struts are defined as: [LE,TE]
        coord[i*2,:] = sorted_strut[0]
        coord[i*2+1,:] = sorted_strut[1]

    return coord


def get_centroid(vertices):
    ''' INPUT: vertices = list of points (np.array)
        OUTPUT: centroid = np.array of centroid'''

    def euclidian_distance(point1,point2):
        return abs(np.linalg.norm((np.array(point1)-np.array(point2))))

    def centroid_algorithm(points):
        ''' INPUT: points = list of points (np.array)
            OUTPUT: centroid = np.array of centroid'''
        
        points = np.array(points)
        distance = 0
        # loop through all possible combinations and whenever one finds a new larger distance, make that the basis for centroid
    
        for i in range(0,len(points)): 
            for j in range(i,len(points)):
                distance_new = euclidian_distance(points[i],points[j])
                if distance_new > distance: #find those points that are the furthest apart
                    distance = distance_new
                    centroid = (points[i]+points[j])/2 #the centroid must be midpoint between these two points
        return np.array(centroid)

    centroids = []
    if len(vertices[0]) == 3: #2 dimensional, array in list
        centroids = centroid_algorithm(vertices)

    else: #3 dimensional, array in list in list
        for i in range(0,len(vertices)):
            centroids.append(centroid_algorithm(vertices[i]))

    return centroids

# Define the function to calculate the cumulative distance along the curve
def cum_dist(points):
    # finds the cumulative distance along the curve
    cum_dist = np.cumsum(np.sqrt(np.sum(np.diff(points, axis=0)**2, axis=1)))
    return np.concatenate([[0], cum_dist])

def get_discretized_canopy(points,percentages):
    # Calculate the cumulative distance along the curve
    distances = cum_dist(points)
    points = np.array(points)

    # Interpolate the x, y, and z values separately
    interp_func_x = interp1d(distances, points[:,0])
    interp_func_y = interp1d(distances, points[:,1])
    interp_func_z = interp1d(distances, points[:,2])

    # Define the distances at which to interpolate
    total_dist = distances[-1]
    interp_distances = []
    for perc in percentages:
        interp_distances.append(perc*total_dist)

    # Interpolate the values at the desired distances
    interp_values = []
    for d in interp_distances:
        interp_values.append([interp_func_x(d), interp_func_y(d), interp_func_z(d)])
    
    return interp_values

########
########
########

#%% Maybe useful later on..

def main_wing(filepath,planar_tol=2,planar_min_points=10,circle_diam_tol=1,circle_min_points=15,circle_split_radius_factor=0.9,LE_radius=50,plate_corners_tol_x=150,plate_corners_tol_y=100,plate_corners_min_z_ratio=0.5):


    # (1) Reading the obj file
    vertices = np.array(read_obj_file_wing(filepath)[0])

    # (2) Splitting the vertices into planar shapes and non_planar shapes
    planars_in_order_list,non_planars = get_planar_in_order(vertices, planar_tol, planar_min_points)
    #get_3D_plot(planars_in_order_list, 'blue',False,1)

    # (3) Extracting the circles from the planar shapes
    circles,circles_radii,circles_centroids,non_circles = get_circle_shapes(planars_in_order_list, circle_diam_tol, circle_min_points)
    #get_3D_plot(circle_list, 'blue',False,1)

    # (3) Getting the strut radius
    strut_radius = get_strut_radius(circles,circles_radii)

    # (4) Splitting the circles based on radius
    big_circles,big_circles_centroids,big_circles_radii,small_circles,small_circles_centroids = split_circles_on_radius(circles,circles_radii,circles_centroids,circle_split_radius_factor*strut_radius)
    #get_3D_plot(almost_tubular_frame, 'blue',False,1)

    # (5) Splitting the LE strut ends away from the rest of the tubular frame
    tubular_frame, tubular_frame_radii, tubular_frame_centroids_surfplan = split_LE_strut_ends(big_circles, big_circles_radii, big_circles_centroids, LE_radius)
    #get_3D_plot(tubular_frame, 'blue',False,1)

    # (6) Seperate the plate corners from the rest of the tubular frame
    plate_corners_surfplan,tubular_frame_centroids_non_plate_corners,N_corners,N_struts,N_panels = split_plate_corners_from_tubular_frame(tubular_frame_centroids_surfplan,plate_corners_tol_x,plate_corners_tol_y,plate_corners_min_z_ratio)
    #get_3D_plot(plate_corners, 'blue',False,1)

    # (7) Transforming the coordinate system
    tubular_frame_centroids_model = transform_coordinate_system_Surfplan_to_Model(tubular_frame_centroids_surfplan)
    plate_corners_model = transform_coordinate_system_Surfplan_to_Model(plate_corners_surfplan)
    #get_3D_plot(plate_corners_model, 'blue',False,1)

    # (8) Getting the connectivity
    wing_conn_i, wing_conn_j, plates = get_kite_plate_connectivity(plate_corners_model,N_corners)
    #get_3D_plot_with_lines(plate_corners_model,wing_conn_i,wing_conn_j,show_legend=False)

    return wing_conn_i, wing_conn_j, plates, tubular_frame_centroids_model, plate_corners_model

