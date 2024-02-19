#%% Notes to self

'''
Important to realize the data structure and its properties here

rib_db_whole is a 2D array that describes all the ribs (used for defining supervector & connectivity)
    The first 14 ribs are the struts
    The last 4 are mid-tips
    Ranked from left-to-right
    Data structure (lst): [rib1, rib2, rib3 ..., ribN]
    Rib structure (lst):  ['strut',LE,Canopy,TE,Particles, Particle Indices]
        LE, Canopy, TE, Particles --> all lists with points as np.arrays(x,y,z)
        Particles Indices --> gives, the indices of supervector, of the particles in that rib
 
supervector is a 1D array of all the particles, filled with np.arrays(x,y,z)
    the first len(supervector_wing) are wing particles
    the last len(supervector_bridle) are bridle particles
    these make up all the particles in the system, and will become the --> points_ini
    Data structure: [(x1,y1,z1),(x2,y2,z2)]

Connectivity is given by
    row_ind: filled with indices (supervector) that are connected to the col_ind ones
    col_ind: filled with indices (supervector) that are connected to the row_ind ones
    conn_indicator: filled with 1s, indicating that the row_ind and col_ind are connected
    
    ci = row_ind, cj = col_ind
'''

# for more tips on plotly, check out: 
# https://plotly.com/python-api-reference/generated/plotly.graph_objects.layout.html#plotly.graph_objects.layout.Title



#%% WING
# (1) variables and reading

filename = 'V9_60C_4_53_M_lines'
#filename = 'V9_10_3d'
filepath = '/home/jellepoland/surfdrive/phd/code/phd_Jelle_Poland/Simulations/data/surfplan_files/kitepower_confidential/V9_60C_4_53_M_lines.obj'

#LE part
tol_planar_distance = 3
tol_min_points = 30

planar_tol=2
planar_min_points=10
circle_diam_tol=1
circle_min_points=15
circle_split_radius_factor=0.9
LE_radius=50
plate_corners_tol_x=150
plate_corners_tol_y=100
plate_corners_min_z_ratio=0.5

# libraries
import re
import numpy as np
import math
import plotly.graph_objects as go
import plotly.offline as offline
import copy
import re # for regular expressions
from sklearn.cluster import DBSCAN # for clustering
import pandas as pd # for dataframes
import plotly.graph_objects as go
import plotly.offline as offline

# functions
from initialisation.V9_60C.functions_wing import *
from initialisation.V9_60C.functions_plotting import *
from initialisation.V9_60C.functions_reading import *
from initialisation.V9_60C.functions_bridles import *

### nice if you want to read out the obj file properties, to decide on rest of analysis
#analyze_obj_file(filepath)

## Reading out the obj-files (selected using the "analyze_obj_file" function)
vertices_LE = read_obj_file_part(filepath,'LE','strut')[0]
vertices_strut = read_obj_file_part(filepath,'strut','16')[0] #also negelcted, as it read out as separete ribs
#vertices_tubes = read_obj_file_part(filepath,'16','upper')[0] #this is neglected entirely
vertices_canopy = read_obj_file_part(filepath,'upper','bridle')[0]
vertices_wing = np.concatenate((vertices_LE,vertices_strut,vertices_canopy),axis=0)  
vertices_bridle = read_obj_file_part(filepath,'bridle','end')[0]
vertices_all = np.concatenate((vertices_wing,vertices_bridle),axis=0)

#get_3D_plot_N_D(2,vertices_wing, 'blue',1)
#get_3D_plot_N_D(3,[vertices_LE,vertices_canopy], 'green',1)


## Getting only the left-part, as the right-part should be exactly the same (and contains some weird problems)
def split_whole_to_left(vertices):
    left_vertices = []
    for i in range(len(vertices)):
        if vertices[i][0] < 0:
            left_vertices.append(np.array(vertices[i]))
    return left_vertices

left_vertices_wing = split_whole_to_left(vertices_wing)
left_vertices_LE = split_whole_to_left(vertices_LE)
left_vertices_canopy =  split_whole_to_left(vertices_canopy)
left_vertices_bridle = split_whole_to_left(vertices_bridle)
left_vertices_strut = split_whole_to_left(vertices_strut)

plot_LE     = [2,vertices_LE, "red", 3]
plot_canopy = [2,vertices_canopy, "black", 1.5]
plot_strut  = [2,vertices_strut, "blue", 2]
plot_title  = "red is LE, black is canopy, blue is strut"
get_3D_plot_N_D_triple(plot_LE,plot_canopy,plot_strut,plot_title)

#%% (2) LE split tubes of LE
# (2) LE

'''
## plotting the result, blue should be tip, red should be mid
plot_LE     = [2,left_vertices_LE, "black", 1]
plot_LE_tip = [2,left_vertices_LE[0:5], "blue", 4]
plot_LE_mid = [2,left_vertices_LE[-5:-1], "red", 4]
plot_title = "Should show the left part of the LE, with blue at the tip and red at mid-span"
get_3D_plot_N_D_triple(plot_LE,plot_LE_tip,plot_LE_mid,plot_title)
'''

def get_planar_in_order(point_list, tol_distance_off_planar):
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

        else:
            if len(planar_shape)>0:
                planar_shape_list.append(planar_shape)
                planar_shape = []

            not_planar_shapes.append(point_list[i])
        
        '''
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
        '''

    return planar_shape_list,not_planar_shapes

# Defining the planars present in left_vertices_LE basedon on the defined tol_distance_off_planar
tol_distance_off_planar = 0.008
left_vertices_LE_planars,LE_other_planar = get_planar_in_order(left_vertices_LE, tol_distance_off_planar)

# plotting
plot_left_vertices_LE_planars = [3,left_vertices_LE_planars, "black", 1]
plot_LE_other = [2,LE_other_planar, "red", 3]
plot_title = "Planar filter: Black is kept, red is filtered out"
#get_3D_plot_N_D_double(plot_left_vertices_LE_planars,plot_LE_other,plot_title)

## Removing the mini tubes from the LE circles
def get_LE_tube_from_planars(left_vertices_LE_planars,LE_other_planar, tol_min_points_for_planar,tol_max_points_for_planar):
    left_vertices_LE_splitted = []
    left_other_splitted = [LE_other_planar]
    for i in range(0,len(left_vertices_LE_planars)):
        if len(left_vertices_LE_planars[i]) > tol_min_points_for_planar and len(left_vertices_LE_planars[i]) < tol_max_points_for_planar:
            left_vertices_LE_splitted.append(left_vertices_LE_planars[i])
        else:
            left_other_splitted.append(left_vertices_LE_planars[i])
            #np.concatenate((LE_other_planar,left_vertices_LE_planars[i]))

    return left_vertices_LE_splitted,left_other_splitted

tol_min_points_for_planar = 30
tol_max_points_for_planar = 80
left_vertices_LE_splitted,left_other_splitted = get_LE_tube_from_planars(left_vertices_LE_planars,LE_other_planar,tol_min_points_for_planar,tol_max_points_for_planar)

# plotting
plot_left_vertices_LE_planars = [3,left_vertices_LE_splitted, "black", 2]
plot_LE_other = [3,left_other_splitted, "red", 3]
plot_title = "Number of points filter: Black is kept, red is filtered out"
get_3D_plot_N_D_double(plot_left_vertices_LE_planars,plot_LE_other,plot_title)

#%% Arrange LE circles from tip to mid

def euclidian_distance(point1,point2):
    return abs(np.linalg.norm((np.array(point1)-np.array(point2))))

## defining a functions that finds the angle (TIP to MID) of a point on the wing, with respect to the mid_point
def get_angle_for_point_on_wing(point_on_wing):
    min_y = min(point_i[1] for point_i in left_vertices_wing)
    angle = np.arctan( (point_on_wing[0]) / (point_on_wing[1]-min_y) )
    return angle

# define sorting on the angle condition
def condition_angle_tip_to_mid(rib_element): #defining the sorting condition ##TODO: this can also be written as a lambda function
    rib_element_centroid = get_centroid(rib_element)
    return get_angle_for_point_on_wing(rib_element_centroid)

## sorting on the set angle condition from tip to mid-span
left_vertices_LE_splitted.sort(key=condition_angle_tip_to_mid)

## plotting to check
# tip should be red, all the tubes should be filtered out
plot_LE = [3,left_vertices_LE_splitted, "black", 1]
plot_LE_tip = [2,left_vertices_LE_splitted[0], "red", 1]
plot_LE_mid = [2,left_vertices_LE_splitted[-1], "blue", 1]
plot_canopy = [2,left_vertices_canopy, "black", 0.5]
plot_title = "LE tubes: Red shows tip, blue shows mid-span"
get_3D_plot_N_D_4times(plot_LE,plot_LE_tip,plot_LE_mid,plot_canopy,plot_title)

#%% defining the database
## defining the database
# ordered from tip-to-mid-span
# db = [type="strut", LE , canopy , TE ]
rib_db = []
## looping through LEs to append them to the rib_db

# define the number of tip LE circles
N_tip_LE_circles = 3 ##TODO: should figure this out 
# concatenate them together ##TODO: remove the hardcoding
#LE_tip_LE_circles = np.concatenate((left_vertices_LE_splitted[1],left_vertices_LE_splitted[2],left_vertices_LE_splitted[3]),axis=0)
##TODO: the above takes into account all LE_tip_circles (better)
# the below only considers one of the circles, leaving out important details
LE_tip_LE_circles = left_vertices_LE_splitted[1]

# append the tip to the rib database
rib_db.append(["strut",LE_tip_LE_circles,[],left_vertices_LE_splitted[0]])

# loop through the rest of the LE circles
for i in range((N_tip_LE_circles+1),len(left_vertices_LE_splitted)):
    if i%2 == 0: # if even
        rib_db.append(["mid",left_vertices_LE_splitted[i],[],[]])
    else: # if odd
        rib_db.append(["strut",left_vertices_LE_splitted[i],[],[]])

### separate function to read out the rib_db
def get_rib_db_index_from_rib_db(type,index,rib_db):
    rib_db_index = []
    for i in range(0,len(rib_db)):
        if rib_db[i][0]==type:
            rib_db_index.append(rib_db[i][index])
    return rib_db_index

## plotting struts stuff
rib_strut_LE = get_rib_db_index_from_rib_db("strut",1,rib_db)
rib_strut_TE = get_rib_db_index_from_rib_db("strut",3,rib_db)

plot_LE = [3,left_vertices_LE_splitted, "black", 1]
plot_strut_LE = [3,rib_strut_LE, 'red',2]
plot_strut_TE = [3,rib_strut_TE, 'blue',2]
plot_canopy = [2,left_vertices_canopy, "black", 0.5]
plot_title = "LE tubes connected to strut: red is LE, blue is TE, black is canopy"
get_3D_plot_N_D_4times(plot_LE,plot_strut_TE,plot_strut_LE,plot_canopy,plot_title)

# Order still checks out here

#%% (3) dealing with struts
# (3) struts

## this list is list with lists, so it needs an additional loop
vertices_strut = read_obj_file_rib(filepath,'strut','16')[3]
left_vertices_strut_splitted = []
for i in range(0,len(vertices_strut)):
    left_vertices_strut_local = split_whole_to_left(vertices_strut[i])
    if len(left_vertices_strut_local) > 0: #only take the ribs that have a left part
        left_vertices_strut_splitted.append(left_vertices_strut_local)

# plot shows that indeed the order is from mid to tip here
# get_3D_plot_N_D(3,left_vertices_strut_splitted[0:5], 'black',1)

## order must be switched around 
## sorting on the set angle condition from tip to mid-span
left_vertices_strut_splitted.sort(key=condition_angle_tip_to_mid)

# new plot shows order is now from tip to mid
plot_strut_vertices_idx_0   = [2,left_vertices_strut_splitted[0], 'red',1]
plot_strut_vertices_idx_1   = [2,left_vertices_strut_splitted[1], 'blue',1]
plot_strut_vertices_idx_rest= [3,left_vertices_strut_splitted[2:], 'black',1]
get_3D_plot_N_D_triple(plot_strut_vertices_idx_0,plot_strut_vertices_idx_1,plot_strut_vertices_idx_rest,"Strut order: red is idx:0, blue is idx:1, black is rest")

#%% Splitting off the Canopy Surface and TEs

def get_canopy_surface_from_strut_rib(rib,tol_TE_split,tol_perc_tubes,tol_point_distance):
        # this code takes one rib, splits it up into multiple "planar shapes"
        # sorts each planar shape on the z-axis
        # using the chord-wise distance it finds which shape corresponds to the canopy and lower surface
        # the list of these points still contains some outliers, which must be filtered out using a "polynomial function fit filter"
        # the function returns the canopy surface points, the lower surface points and the points that are not part of the canopy or lower surface

        ### sort on z-axis (from TE to LE)
        rib_sorted_on_z = sorted(rib, key=lambda point: point[2])

        ### loop over each point from TE to LE (low z to high z)
        rib_TE_circle,rib_canopy,rib_other = [],[],[]
        flag = True
        for point in rib_sorted_on_z:
            ## Make flag False when z/|z| > tol_TE_split
            if (point[2]-rib_sorted_on_z[0][2]) > tol_TE_split: ##TODO: could change this to planar constraint instead of a defined length, both require input tho...
                flag = False

            point = np.array(point)
            ## if still on the TE_circle
            if flag: 
                rib_TE_circle.append(point)
            ## if above the TE_circle and below tol_perc_tubes: "0.6" z/|z| (rib length)
            elif point[2] < tol_perc_tubes*rib_sorted_on_z[-1][2]: ##TODO: should ideally be without a hardcoded edge, but not sure how yet
                rib_canopy.append(point)
            ## if above the TE_circle and above tol_perc_tubes: "0.6" z/|z| (rib length), append to other
            else: ##TODO: better would be to split these points here based on a line constraint, but lets' focus on other things first
                rib_other.append(point)

        j = len(rib_canopy)-1 #start at the end of the list
        rib_other_new = []

        for i in range(0,len(rib_other)):
            ## check if point is close to last, most upper point on the canopy surface ##TODO: should be a line constraint instead of a point constraint
            if euclidian_distance(rib_canopy[j],rib_other[i]) < tol_point_distance:
                rib_canopy.append(rib_other[i])
                j = j+1
            else:
                rib_other_new.append(rib_other[i])

        # for those ribs that don't have the canopy surface somehow
        # we check if the length is atleast 0.6 of the total rib length
        ## if canopy_z_min - canopy_z_max < 0.6*rib_z_max, than we can assume that the canopy surface is not there ##TODO: also a bit vague
        if abs(rib_canopy[0][2] - rib_canopy[-1][2]) < tol_perc_tubes*rib_sorted_on_z[-1][2]:
            rib_other = np.concatenate((rib_other,rib_canopy),axis=0)
            rib_canopy = []

        return rib_TE_circle,rib_canopy,rib_other_new

# setting the tolerances to split the strut
tol_TE_split = 0.02
tol_perc_tubes = 0.6
tol_point_distance = 0.02

strut_TE_list,strut_canopy_list,strut_other_list = [],[],[]
for i in range(0,len(left_vertices_strut_splitted)):
    strut_TE,strut_canopy,strut_other = get_canopy_surface_from_strut_rib(left_vertices_strut_splitted[i],tol_TE_split,tol_perc_tubes,tol_point_distance)
    
    if len(strut_TE) >0: strut_TE_list.append(strut_TE)
    if len(strut_canopy) >0: strut_canopy_list.append(strut_canopy)
    if len(strut_other) >0:strut_other_list.append(strut_other)

strut_TE_list.sort(key=condition_angle_tip_to_mid)
strut_canopy_list.sort(key=condition_angle_tip_to_mid)

# new plot shows order is now from tip to mid
plot_strut_TE_idx_0   = [2,strut_TE_list[0], 'red',1]
plot_strut_TE_idx_1   = [2,strut_TE_list[1], 'blue',1]
plot_strut_TE_idx_rest= [3,strut_TE_list[2:], 'black',1]
#get_3D_plot_N_D_triple(plot_strut_TE_idx_0,plot_strut_TE_idx_1,plot_strut_TE_idx_rest,"TE order: red is idx:0, blue is idx:1, black is rest")

plot_strut_canopy_idx_0   = [2,strut_canopy_list[0], 'red',1]
plot_strut_canopy_idx_1   = [2,strut_canopy_list[1], 'blue',1]
plot_strut_canopy_idx_rest= [3,strut_canopy_list[2:], 'black',1]
strut_LE = get_rib_db_index_from_rib_db("strut",1,rib_db)
plot_strut_LE = [3,strut_LE, "green", 1]
get_3D_plot_N_D_4times(plot_strut_LE,plot_strut_canopy_idx_0,plot_strut_canopy_idx_1,plot_strut_canopy_idx_rest,"Canopy order: red is idx:0, blue is idx:1, black is rest")

#%% Add the canopy's and TEs to the rib_db

##TODO: remove the hard-coding
# loop through each rib_db
for i in range(0,len(rib_db)):
    # only when rib_db is a strut 
    if rib_db[i][0] == "strut":
        # only if not dealing with the tip-strut /most outer strut
        if i > 1:
            # add the TE to the rib_db 
            idx_TE = int(i/2)-1
            rib_db[i][3] = strut_TE_list[idx_TE]

        # only when not the outer, nor the second to outer strut
        if i > 3:
            # add the canopy to the rib_db
            idx_canopy = int(i/2)-2
            rib_db[i][2] = strut_canopy_list[idx_canopy]

plot_canopy         = [2,left_vertices_canopy, 'black',1]
plot_strut_LE       = [3,get_rib_db_index_from_rib_db("strut",1,rib_db), "red", 3]
plot_strut_canopy   = [3,get_rib_db_index_from_rib_db("strut",2,rib_db), 'green',3]
plot_strut_TE       = [3,get_rib_db_index_from_rib_db("strut",3,rib_db), 'blue',3]
get_3D_plot_N_D_4times(plot_canopy,plot_strut_TE,plot_strut_LE,plot_strut_canopy,"LE: red, covered canopy: green (don't need perfect coverage), TE: blue, rest canopy is black")

def get_check_rib_db(rib_db,name="strut"):
    for i,rib in enumerate(rib_db):
        if rib_db[i][0] == name:
            print(rib_db[i][0])
            print("LE     :",len(rib_db[i][1]))
            print("canopy :",len(rib_db[i][2]))
            print("TE     :",len(rib_db[i][3]))
    return

get_check_rib_db(rib_db)

#%% (4) FILLING THE CANOPY OF THE CANOPY-EMPTY STRUTS (TWO MOST OUTER STRUTS)) 
## (4) FILLING THE CANOPY OF THE CANOPY-EMPTY STRUTS (TWO MOST OUTER STRUTS)) 

## Defining a distance to plane function
def distance_point_to_plane(point, plane):
    px1, py1, pz1 = plane[0]
    px2, py2, pz2 = plane[1]
    px3, py3, pz3 = plane[2]
    nx, ny, nz = np.cross((px1 - px2, py1 - py2, pz1 - pz2), (px1 - px3, py1 - py3, pz1 - pz3))
    ax, ay, az = point
    distance = abs((ax - px1) * nx + (ay - py1) * ny + (az - pz1) * nz) / np.sqrt(nx**2 + ny**2 + nz**2)
    return distance

## Defining a distance to line function 
#def distance_point_to_line(point, line):
#    # Calculate the equation of the line
#    line_direction = line[1] - line[0]
#    line_direction = line_direction / np.linalg.norm(line_direction)
#    line_normal = np.cross(line_direction, [0, 0, 1])
#    line_normal = line_normal / np.linalg.norm(line_normal)
#
#    # Calculate the distance from the point to the line
#    distance = np.dot(point - line[0], line_normal)
#
#    return np.abs(distance)

## Defining a distance to line function 
def lineseg_dists(p, a, b):
    # Handle case where p is a single point, i.e. 1d array.
    p = np.atleast_2d(p)

    ##TODO for you: consider implementing @Eskapp's suggestions
    if np.all(a == b):
        return np.linalg.norm(p - a, axis=1)

    # normalized tangent vector
    d = np.divide(b - a, np.linalg.norm(b - a))

    # signed parallel distance components
    s = np.dot(a - p, d)
    t = np.dot(p - b, d)

    # clamped parallel distance
    h = np.maximum.reduce([s, t, np.zeros(len(p))])

    # perpendicular distance component, as before
    # note that for the 3D case these will be vectors
    c = np.cross(p - a, d)

    # use hypot for Pythagoras to improve accuracy
    return np.hypot(h, c)

## defining the mid point
min_y = min(point_i[1] for point_i in left_vertices_wing)
ave_z = ((max(point[2] for point in left_vertices_canopy) + min(point[2] for point in left_vertices_canopy))/2) #chord
mid_point = np.array([0,min_y,ave_z])

## define sorting the angle calculations (TE to LE)
def get_angle_TE_to_LE(point): #defining the sorting condition, which is z over xy
    return np.arctan( point[2] / np.sqrt( (point[0]**2)+((point[1]-min_y)**2) ) )

tol_angle_canopy_to_strut   = 0.5    #0.5 #tolerance angle in degrees
tol_off_line = 0.015 #0.015 #% tolerance off line (for non filled struts)

new_canopy_lst = [] #storing the runs for plotting purposes
##################
## loop through each rib
for i,rib in enumerate(rib_db):
    if rib_db[i][0] == 'strut' and len(rib_db[i][2]) == 0: ##TODO: remove hardcoded occurence of non-filled strut i == 2            
    #elif rib_db[i][0] == 'strut': #if the rib is strut, but does not yet contain canopy points
        print("i = ",i," is a strut, but does not yet contain canopy points")
        print("ini len: ",len(rib_db[i][2]))

        # 1. find additional canopy points
        # 2. finding the canopy point that's closest to the LE outer point
        # 3. finding the line between the LE canopy point and TE_circle outer point
        # 4. Use proximity to line to filter those canopy points that are part of the new canopy
        #
        # Last step here would be to find the 'mid' points.
        # you can probably take the closest point to the point between TE_circles of accompayning struts, to find a TE point
        # Than do the same trick as above to find the canopy points, but with the change of taking the TE point rather than the outer TE_circle point 

        ## (1) finding additional_canopy_points that are within (tip-to-mid) angle range
        additional_canopy_points = []
        LE_centroid_angle = get_angle_for_point_on_wing(get_centroid(rib_db[i][1]))
        for point in left_vertices_canopy:
            point_angle = get_angle_for_point_on_wing(point) 
            ## if the angle is within the tolerance, append the point to the canopy surface
            if abs(point_angle - LE_centroid_angle) < np.deg2rad(tol_angle_canopy_to_strut):
                additional_canopy_points.append(point)

        ## (2) Get additional_canopy_point_closest_to_LE point 
        # sort LE to have max. distance from mid-point as [-1] index
        def cond_dist_mid(point):
            return euclidian_distance(point,mid_point)
        rib_db[i][1] = sorted(rib_db[i][1],key=cond_dist_mid)
        LE_outer_point = rib_db[i][1][-1] #define the LE outer point

        # sort canopy list to have min. distance from LE outer point as [0] index
        def cond_dist_LE_outer(point):
            return euclidian_distance(point,LE_outer_point)
        additional_canopy_points.sort(key=cond_dist_LE_outer)
        additional_canopy_point_closest_to_LE = additional_canopy_points[0]
        
        ## (3) finding the line between LE canopy point and TE_circle outer point
        # sort TE to have max. distance from mid-point as [-1] index
        rib_db[i][3].sort(key=cond_dist_mid)
        TE_outer_point = rib_db[i][3][-1] #define the TE outer point
        line = [additional_canopy_point_closest_to_LE, TE_outer_point]
        length_rib_canopy = euclidian_distance(line[0],line[1])
        
        ## (4) Use proximity to line to filter those canopy points that are part of the new canopy
        new_canopy = [] #defining the new canopy, first as empty list ofcourse
        ## finding those additional points that are close enough to the line
        for j,point in enumerate(additional_canopy_points):
            ## if close enough to line, then add to new canopy
            if np.linalg.norm((lineseg_dists(additional_canopy_points[j],line[0],line[1])))< tol_off_line*length_rib_canopy: ##TODO: hardcoded
                new_canopy.append(np.array(additional_canopy_points[j])) ## append the point to new_canopy

        ## Append the new canopy to the rib_db
        rib_db[i][2] = new_canopy
        
        new_canopy_lst.append(new_canopy)
        print("FIN len: ",len(rib_db[i][2]))
        print(" ")

plot_strut_LE           = [3,get_rib_db_index_from_rib_db("strut",1,rib_db), "red", 3]
plot_strut_canopy_new   = [3,new_canopy_lst, 'green',5]
plot_strut_TE           = [3,get_rib_db_index_from_rib_db("strut",3,rib_db), 'blue',3]
plot_vertices_canopy    = [2,left_vertices_canopy, 'black',1]
get_3D_plot_N_D_4times(plot_strut_TE,plot_strut_LE,plot_strut_canopy_new,plot_vertices_canopy,"LE: red, new canopy: green (don't need perfect coverage), TE: blue, rest canopy is black")


#%% FILLING THE CANOPYS A BIT MORE, using the vertices_canopy

##################
## STRUT CANOPY ##
##################

## filtered min-distance between line and canopy points

#N_steps = 30 #number of steps to take from TE to LE
tol_angle_canopy_to_strut   = 0.5    #0.5 #tolerance angle in degrees
tol_off_plane               = 0.02 #0.02 #% tolerance off plane (for filled struts)


new_canopy_points_lst = [] #list to store the new canopy points
## loop through each rib
for i,rib in enumerate(rib_db):
   
    if rib_db[i][0] == 'strut' and len(rib_db[i][2]) != 0:#i != 2: #if the rib is a prefilled strut

        print("ini len: ",len(rib_db[i][2]))
        ## sorting on new angle condition from TE to LE (so not tip to mid angle condition)
        rib_db[i][2] = sorted(rib_db[i][2], key=get_angle_TE_to_LE)

        N_canopy_points = len(rib_db[i][2]) #length of the canopy points
        #print(i)
        #p_beg    = rib_db[i][2][0] #first point of the canopy
        #p_middle = rib_db[i][2][int(N_canopy_points/2)] #middle point of the canopy
        #p_end    = rib_db[i][2][-1] #last point of the canopy
        plane = [rib_db[i][2][0],rib_db[i][2][int(N_canopy_points/2)],rib_db[i][2][-1]] #define the plane based on the first, middle and last point of the canopy
        length_rib_canopy = euclidian_distance(rib_db[i][2][0],rib_db[i][2][-1]) #length of existing canopy
        new_canopy_points = [] # The newly found points, that will be added to the last canopy points

        ## finding additional_canopy_points that are within tip-to-mid angle range
        additional_canopy_points = []
        LE_centroid_angle = get_angle_for_point_on_wing(get_centroid(rib_db[i][1]))
        for point in left_vertices_canopy:
            point_angle = get_angle_for_point_on_wing(point) 
            ## if the angle is within the tolerance, append the point to the canopy surface
            if abs(point_angle - LE_centroid_angle) < np.deg2rad(tol_angle_canopy_to_strut):
                additional_canopy_points.append(point)

        ## From the additional canopy points finds that satisfy
        ## By filtering on proximity to plane (defined by existing TE,mid,LE points)
        for j,point in enumerate(additional_canopy_points):
            ## if close enough to plane, then add to new canopy
            if abs(distance_point_to_plane(additional_canopy_points[j],plane))< tol_off_plane*length_rib_canopy: ##TODO: hardcoded
                new_canopy_points.append(np.array(additional_canopy_points[j])) ## append the point to new_canopy

        new_canopy_points_lst.append(new_canopy_points)
        #print(np.shape(rib_db[i][2]))
        #print(np.shape(new_canopy_points))
        ## Append the new canopy to the rib_db
        if len(new_canopy_points) > 0:
            rib_db[i][2] = np.concatenate((rib_db[i][2],new_canopy_points),axis=0)

        print("FIN len: ",len(rib_db[i][2]))
        print(" ")

## plotting struts stuff
rib_strut_LE = get_rib_db_index_from_rib_db("strut",1,rib_db)
rib_strut_canopy =  get_rib_db_index_from_rib_db("strut",2,rib_db)

##TODO: This because some are still 0
rib_strut_canopy_new = []
for i in range(0,len(rib_strut_canopy)):
    if len(rib_strut_canopy[i])>0:
        rib_strut_canopy_new.append(rib_strut_canopy[i])

plot_strut_LE           = [3,get_rib_db_index_from_rib_db("strut",1,rib_db), "red", 3]
plot_strut_canopy       = [3,get_rib_db_index_from_rib_db("strut",2,rib_db), 'green',3]
plot_strut_TE           = [3,get_rib_db_index_from_rib_db("strut",3,rib_db), 'blue',3]
plot_strut_canopy_new   = [3,new_canopy_points_lst, 'black',5]
get_3D_plot_N_D_4times(plot_strut_TE,plot_strut_LE,plot_strut_canopy,plot_strut_canopy_new,"LE: red, covered canopy: green (don't need perfect coverage), TE: blue, rest canopy is black")

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
# The above works, and extracts desired strut: LE, canopy, TE points
# The below will deal with getting the desired points from each strut rib

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################

rib_db_old = copy.deepcopy(rib_db)

#%% Defining new circles from the old circles
## Firstly defining new circles, based on the old ones, as they don't describe a full circle

rib_db = copy.deepcopy(rib_db_old)

def rotate_vector(v, axis, angle):
    """
    Rotate a vector v around an axis by a given angle.
    """
    axis = axis / np.linalg.norm(axis)
    a = np.cos(angle / 2)
    b, c, d = -axis * np.sin(angle / 2)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    rotation_matrix = np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                                [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                                [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])
    return np.dot(rotation_matrix, v)

def find_vector3(vector1, vector2, circle_center, circle_radius, angle):
    # Calculate the angle between vector1 and vector2
    v1 = vector1 - circle_center
    v2 = vector2 - circle_center
    #angle = np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))

    # Rotate vector1 around the circle center by the same angle as vector2
    axis = np.cross(v1, v2)
    vector3 = rotate_vector(v1, axis, angle) + circle_center

    # Ensure that vector3 is on the circle
    vector3 = circle_center + circle_radius * (vector3 - circle_center) / np.linalg.norm(vector3 - circle_center)

    return vector3

def get_new_circle(N, points):
    # N is the number of points on the new circle
    # Find the center and radius of the circle
    center = get_centroid(points)
    radius = np.mean(np.linalg.norm(points - center, axis=1))

    # Define the angles of the new points
    angles = np.linspace(0, 2*np.pi, N)

    new_points = []
    for angle in angles:
        new_points.append(find_vector3(points[0], points[10], center, radius, angle))

    return new_points

for i,rib in enumerate(rib_db):# Print the new points
    if rib_db[i][0] == "strut":
        ##TODO: hardcoded number of points
        #transform LE into new circle with 100 points
        rib_db[i][1] = get_new_circle(100,rib_db[i][1])
        #transform canopy into new circle with 100 points
        rib_db[i][3] = get_new_circle(100,rib_db[i][3])
    
    if rib_db[i][0] == 'mid':
        #transform LE into new circle with 100 points
        rib_db[i][1] = get_new_circle(100,rib_db[i][1])

## plotting struts stuff
plot_strut_LE_old   = [3,get_rib_db_index_from_rib_db("strut",1,rib_db_old), 'red' ,2]
plot_strut_LE       = [3,get_rib_db_index_from_rib_db("strut",1,rib_db), 'black' ,1.5]
plot_strut_TE_old   = [3,get_rib_db_index_from_rib_db("strut",3,rib_db_old), 'blue' ,2]
plot_strut_TE       = [3,get_rib_db_index_from_rib_db("strut",3,rib_db), 'black' ,1.5]

get_3D_plot_N_D_4times(plot_strut_TE,plot_strut_LE,plot_strut_TE_old,plot_strut_LE_old,"Should have complete circles now, LE: red, canopy: green, TE: blue")

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################  
# With all wing circles, lines etc. we can now extract the
# desired points using the code below 
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################


#%% Finding each point on each rib
## Finding each point on each rib

strut_rib_particles = [] # list of the particles

for i in range(len(rib_db)):
    if rib_db[i][0] == "strut":

        ##TODO: this is only added for testing (to not have to run the whole script each-run)
        rib_db[i] = [rib_db[i][0],rib_db[i][1],rib_db[i][2],rib_db[i][3]]

        strut_rib_particles = [] # list of the particles

        ## [0] LE point closest to mid-point (on the lower side)
        # (a) sort LE of current strut based on distance to mid-point and take[0]
        #def cond_dist_mid(point):
        #    return euclidian_distance(point,mid_point)
        #strut_rib_particles.append(sorted(rib_db[i][1],key=cond_dist_mid)[0])
        ##TODO: the point above is removed as there is a bridle attachment point present.

        ## [1] LE point in stagnation point
        # (a) sort on LE to TE angle and take[-1]
        strut_rib_particles.append(sorted(rib_db[i][1], key=get_angle_TE_to_LE)[-1])

        ## [2] LE point closest to most LE canopy point
        # (a) find the canopy_close_to_LE_point
        canopy_close_to_LE_point = sorted(rib_db[i][2], key=get_angle_TE_to_LE)[-1]

        # (b) sort LE on distance to the point found above and take[0]
        def LE_cond_dist_to_canopy(point):
            return euclidian_distance(point,canopy_close_to_LE_point)
        
        strut_rib_particles.append(sorted(rib_db[i][1],key=LE_cond_dist_to_canopy)[0])

        ## [3],[4],[5] 25/50/75% canopy LE to TE
        # (a) sort the canopy point lists from LE to TE
        canopy_points_to_be_discretized = sorted(rib_db[i][2], key=get_angle_TE_to_LE, reverse=True)

        # (b) find the TE point closest to the canopy point
        def TE_cond_dist_to_canopy(point):
            return euclidian_distance(point,canopy_points_to_be_discretized[0])
        
        TE_close_to_canopy =sorted(rib_db[i][3],key=TE_cond_dist_to_canopy)[0]

        # (c) append the TE point, to ensure it goes from LE to TE (assuming LE point is always there ##TODO: fix this assumption)
        canopy_points_to_be_discretized.append(TE_close_to_canopy)

        # (d) use function above to find the points ##TODO: hardcoded..
        interp_values = get_discretized_canopy(np.array(canopy_points_to_be_discretized),[0.25,0.5,0.75]) ##TODO: hardcoded..

        strut_rib_particles.append(interp_values[0])
        strut_rib_particles.append(interp_values[1])
        strut_rib_particles.append(interp_values[2])

        ## [6] TE point closest to canopy
        # (a) this is already found in the last step, so just append
        strut_rib_particles.append(TE_close_to_canopy)

        ## [7] TE point closest to mid-point (on the lower side)
        # (a) sort TE of current strut based on distance to mid-point and take[0]
        strut_rib_particles.append(sorted(rib_db[i][3],key=cond_dist_mid)[0])

        # define the number of strut rib particles
        N_strut_particles = len(strut_rib_particles)

        rib_db[i].append(strut_rib_particles)

## plotting struts stuff
rib_strut_LE = get_rib_db_index_from_rib_db("strut",1,rib_db)
rib_strut_canopy =  get_rib_db_index_from_rib_db("strut",2,rib_db)
rib_strut_TE = get_rib_db_index_from_rib_db("strut",3,rib_db)
rib_strut_particles = get_rib_db_index_from_rib_db("strut",4,rib_db)

plot_strut_LE = [3,rib_strut_LE, 'red',1]
plot_strut_canopy = [3,rib_strut_canopy, 'green',2]
plot_strut_TE = [3,rib_strut_TE, 'blue',1]
plot_strut_particles = [3,rib_strut_particles, 'black',4]

get_3D_plot_N_D_4times(plot_strut_TE,plot_strut_LE,plot_strut_canopy,plot_strut_particles,"black shows strut discretisation")


#%% Defining the mid-bridle special points
## mid-bridle special points

# retrieving the LEs of these points
mid_tip_LE_1 = left_vertices_LE_splitted[3]
mid_tip_LE_2 = left_vertices_LE_splitted[4]
mid_tip_LE_lst = [mid_tip_LE_1,mid_tip_LE_2]

# retrieving the canopy and TE of these points

# 1. find additional canopy points
# 2. finding the canopy point that's closest to the LE outer point
# 3. finding the line between the LE canopy point and TE_circle outer point
# 4. Use proximity to line to filter those canopy points that are part of the new canopy
#
# Last step here would be to find the 'mid' points.
# you can probably take the closest point to the point between TE_circles of accompayning struts, to find a TE point
# Than do the same trick as above to find the canopy points, but with the change of taking the TE point rather than the outer TE_circle point

tol_angle = 1 #tolerance angle in degrees
tol_off_plane = 0.02 #% tolerance off plane (for filled struts)
tol_off_line = 0.01 #% tolerance off line (for non filled struts)

mid_tip_results = [] #initializingg

for mid_tip_LE_i in mid_tip_LE_lst:

    ## (1) finding additional_canopy_points that are within (tip-to-mid) angle range
    additional_canopy_points = []
    LE_centroid_angle = get_angle_for_point_on_wing(get_centroid(mid_tip_LE_i))
    for point in left_vertices_canopy:
        point_angle = get_angle_for_point_on_wing(point) 
        ## if the angle is within the tolerance, append the point to the canopy surface
        if abs(point_angle - LE_centroid_angle) < np.deg2rad(tol_angle):
            additional_canopy_points.append(point)

    ## (2) Get additional_canopy_point_closest_to_LE point 
    # sort LE to have max. distance from mid-point as [-1] index
    def cond_dist_mid(point):
        return euclidian_distance(point,mid_point)
    mid_tip_LE_i = sorted(mid_tip_LE_i,key=cond_dist_mid)
    LE_outer_point = mid_tip_LE_i #define the LE outer point

    # sort canopy list to have min. distance from LE outer point as [0] index
    def cond_dist_LE_outer(point):
        return euclidian_distance(point,LE_outer_point)
    additional_canopy_points.sort(key=cond_dist_LE_outer)
    additional_canopy_point_closest_to_LE = additional_canopy_points[0]
    mid_tip_TE_i = additional_canopy_points[-1]

    ## test plot, blue = TE, purple is closest canopy point to LE.
    #get_3D_plot_N_D_4times([2,additional_canopy_points, 'green',1],[2,mid_tip_LE_i, 'red',1],[2,[mid_tip_TE_i], 'blue',4],[2,[additional_canopy_point_closest_to_LE], 'purple',4])

    ## (3) finding the line between LE canopy point and TE
    # sort TE to have max. distance from mid-point as [-1] index
    line = [additional_canopy_point_closest_to_LE, mid_tip_TE_i]
    length_rib_canopy = euclidian_distance(line[0],line[1])
            
    ## (4) Use proximity to line to filter those canopy points that are part of the new canopy
    new_canopy = [] #defining the new canopy, first as empty list ofcourse
    ## finding those additional points that are close enough to the line
    for j,point in enumerate(additional_canopy_points):
        ## if close enough to line, then add to new canopy
        if np.linalg.norm((lineseg_dists(additional_canopy_points[j],line[0],line[1])))< tol_off_line*length_rib_canopy: ##TODO: hardcoded
            new_canopy.append(additional_canopy_points[j]) ## append the point to new_canopy

    ## Append the new canopy to the rib_db
    mid_tip_canopy_i = copy.deepcopy(new_canopy)

    ## test plot, blue = TE, purple is canopy points
    #get_3D_plot_N_D_4times([2,additional_canopy_points, 'green',1],[2,mid_tip_LE_i, 'red',1],[2,[mid_tip_TE_i], 'blue',4],[2,mid_tip_canopy_i, 'purple',4])

    ##############################################################3
    ## Getting also the particles at once
    mid_tip_particles_i= []
    ## [0] LE point closest to mid-point (on the lower side)
    # (a) sort LE of current strut based on distance to mid-point and take[0]
    #def cond_dist_mid(point):
    #    return euclidian_distance(point,mid_point)
    #strut_rib_particles.append(sorted(rib_db[i][1],key=cond_dist_mid)[0])
    ##TODO: the point above is removed as there is a bridle attachment point present.

    ## [1] LE point in stagnation point
    # (a) sort on LE to TE angle and take[-1]
    mid_tip_particles_i.append(sorted(mid_tip_LE_i, key=get_angle_TE_to_LE)[-1])

    ## [2] LE point closest to most LE canopy point
    # (a) find the canopy_close_to_LE_point
    canopy_close_to_LE_point = sorted(mid_tip_canopy_i, key=get_angle_TE_to_LE)[-1]

    # (b) sort LE on distance to the point found above and take[0]
    def LE_cond_dist_to_canopy(point):
        return euclidian_distance(point,canopy_close_to_LE_point)
    
    mid_tip_particles_i.append(sorted(mid_tip_LE_i,key=LE_cond_dist_to_canopy)[0])

    ## [3],[4],[5] 25/50/75% canopy LE to TE
    # (a) sort the canopy point lists from LE to TE
    canopy_points_to_be_discretized = sorted(mid_tip_canopy_i, key=get_angle_TE_to_LE, reverse=True)

    # (b) find the TE point closest to the canopy point
    TE_close_to_canopy = mid_tip_TE_i

    # (c) append the TE point, to ensure it goes from LE to TE (assuming LE point is always there ##TODO: fix this assumption)
    canopy_points_to_be_discretized.append(TE_close_to_canopy)

    # (d) use function above to find the points ##TODO: hardcoded..
    interp_values = get_discretized_canopy(np.array(canopy_points_to_be_discretized),[0.25,0.5,0.75]) ##TODO: hardcoded..

    mid_tip_particles_i.append(interp_values[0])
    mid_tip_particles_i.append(interp_values[1])
    mid_tip_particles_i.append(interp_values[2])

    ## [6] TE point closest to canopy
    # (a) this is already found in the last step, so just append
    mid_tip_particles_i.append(TE_close_to_canopy)

    # Getting to a whole wing by storing, deepcopying and changing the x location
    mid_tip_left_i = ["mid-tip",mid_tip_LE_i,mid_tip_canopy_i,[mid_tip_TE_i],mid_tip_particles_i]
    mid_tip_right_i = copy.deepcopy(mid_tip_left_i)

    mid_tip_right_i[1] = [np.array([-point[0],point[1],point[2]]) for point in mid_tip_right_i[1]]
    mid_tip_right_i[2] = [np.array([-point[0],point[1],point[2]]) for point in mid_tip_right_i[2]]
    mid_tip_right_i[3] = [np.array([-point[0],point[1],point[2]]) for point in mid_tip_right_i[3]]
    mid_tip_right_i[4] = [np.array([-point[0],point[1],point[2]]) for point in mid_tip_right_i[4]]

    mid_tip_results.append(mid_tip_left_i)
    mid_tip_results.append(mid_tip_right_i)

# weird order is to ensure sorting for the whole wing from left-to-right
rib_db_mid_tip1_left  = mid_tip_results[0] #most outer
rib_db_mid_tip2_left  = mid_tip_results[2] #most inner
rib_db_mid_tip1_right = mid_tip_results[3] #most inner
rib_db_mid_tip2_right = mid_tip_results[1] #most outer

## Verification plots
plot_mid_tip_LE = [2, np.concatenate((rib_db_mid_tip1_left[1],rib_db_mid_tip2_left[1],rib_db_mid_tip1_right[1],rib_db_mid_tip2_right[1])), 'red',1]
plot_mid_tip_particles = [2, np.concatenate((rib_db_mid_tip1_left[4],rib_db_mid_tip2_left[4],rib_db_mid_tip1_right[4],rib_db_mid_tip2_right[4])), 'black',3]
get_3D_plot_N_D_double(plot_mid_tip_LE,plot_mid_tip_particles)
#get_3D_plot_N_D_triple([2,rib_db_mid_tip1_left[1], 'red',1],[2,rib_db_mid_tip1_left[3], 'red',4],[2,rib_db_mid_tip1_left[4], 'black',3])
#get_3D_plot_N_D_triple([2,rib_db_mid_tip2_left[1], 'red',1],[2,rib_db_mid_tip2_left[3], 'red',4],[2,rib_db_mid_tip2_left[4], 'black',3])



#%% Get the whole wing, by copying and mirroring over x-axis

## Defining the left_side (forcing a deepcopy!!)
rib_db_left = copy.deepcopy(rib_db)

## Defining the right-side, to get to a whole wing
rib_db_right = copy.deepcopy(rib_db) #making a DEEPcopy
rib_db_right = rib_db_right[::-1] # reverse the list

for i in range(len(rib_db_right)):
    
    ##TODO: should also do this for the MID ofcourse
    if rib_db_right[i][0] == "strut":
        # mirror over x, using list comprehension
        rib_db_right[i][1] = [np.array([-point[0],point[1],point[2]]) for point in rib_db_right[i][1]]
        rib_db_right[i][2] = [np.array([-point[0],point[1],point[2]]) for point in rib_db_right[i][2]]
        rib_db_right[i][3] = [np.array([-point[0],point[1],point[2]]) for point in rib_db_right[i][3]]
        rib_db_right[i][4] = [np.array([-point[0],point[1],point[2]]) for point in rib_db_right[i][4]]

## Concatenating the left and right side
rib_db_whole = np.concatenate((rib_db_left,rib_db_right),axis=0)
import copy ##TODO: move this import to the top
rib_db_whole_old = copy.deepcopy(rib_db_whole) #making a DEEPcopy

##TODO: once you are using mid-points, you must add the one in the middle plane as well!!!
## Deleting the mid-points for now, since they are not needed
rib_db_whole = []
for i in range(len(rib_db_whole_old)):
    if rib_db_whole_old[i][0]=='strut':
        rib_db_whole.append(rib_db_whole_old[i])

###########################################
###########################################
##TODO: temporary fix for the mid-tip points
# appending to the back of the list, in order LEFT-TO-RIGHT over the whole wing
rib_db_whole.append(rib_db_mid_tip1_left)
rib_db_whole.append(rib_db_mid_tip2_left)
rib_db_whole.append(rib_db_mid_tip1_right)
rib_db_whole.append(rib_db_mid_tip2_right)

### making the circles full again
for i,rib in enumerate(rib_db_whole):# Print the new points
    if rib_db_whole[i][0] == 'mid-tip':
        rib_db_whole[i][1] = get_new_circle(100,rib_db_whole[i][1])

## plotting struts stuff
rib_strut_LE = get_rib_db_index_from_rib_db("strut",1,rib_db_whole)
rib_strut_canopy =  get_rib_db_index_from_rib_db("strut",2,rib_db_whole)
rib_strut_TE = get_rib_db_index_from_rib_db("strut",3,rib_db_whole)
rib_strut_particles = get_rib_db_index_from_rib_db("strut",4,rib_db_whole)

plot_mid_tip_LE = [3, get_rib_db_index_from_rib_db("mid-tip",1,rib_db_whole), 'red',1]
plot_strut_LE = [3,rib_strut_LE, 'red',1]
plot_strut_canopy = [3,rib_strut_canopy, 'green',1]
plot_strut_TE = [3,rib_strut_TE, 'blue',1]
plot_strut_particles = [3,rib_strut_particles, 'black',4]

#get_3D_plot_N_D_4times(plot_strut_TE,plot_strut_LE,plot_strut_canopy,plot_strut_particles)
get_3D_plot_N_D_5times(plot_strut_LE,plot_strut_particles,plot_mid_tip_LE,plot_mid_tip_particles,plot_strut_TE)

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
# With the whole wing, all particles defined
# We can start with connectivity
# Because we need the bridles, to define the whole rib, we'll do that first
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################



#%% Extracting bridles

obj_raw = read_obj_file_bridles_faces(filepath)
vertices_bridle = np.array(obj_raw[0])
faces_raw = obj_raw[2] # extracting the vertices that form a face
#get_3D_plot(vertices, 'red',1)

## (2) Splitting the vertices on proximity
# tol_variables
tol_proximity_bridle_vertices        = 0.005 #m
tol_min_points_bridle_knot_or_pulley = 7    #[-]
tol_proximity_bridle_vertices_pulley = 0.07 #m
tol_min_points_pulley                = 60   #[-]

# Defining both KNOTS & PULLEYS
proximity_list          = split_for_proximity(vertices_bridle,tol_proximity_bridle_vertices,tol_min_points_bridle_knot_or_pulley)
proximity_list_pulley   = split_for_proximity(vertices_bridle,tol_proximity_bridle_vertices_pulley,tol_min_points_pulley)

plot_vertices_pulleys   = [3,proximity_list_pulley, 'red',1]
plot_vertices_bridle    = [3,proximity_list, 'blue',6]
get_3D_plot_N_D_double(plot_vertices_bridle,plot_vertices_pulleys)

##TODO: write a function that removes the proximity_list points that have more than 8 neighbours

#%% (TAKES LONG) Getting the indices of the proximity list
##  (TAKES LONG) Getting the indices of the proximity list
proximity_list_index = get_indices_from_proximity_list(proximity_list,vertices_bridle)

#%% CONNECTIVITY  
## (4) Getting the connectivity of the bridle
bridle_connectivity, bridle_connectivity_i, bridle_connectivity_j = get_bridle_connectivity(faces_raw,proximity_list_index)

## (5) Getting the midpoints of each bridle vertex
midpoints_surfplan = get_midpoints_of_proximity_list(proximity_list)
midpoints_surfplan_pulley = get_midpoints_of_proximity_list(proximity_list_pulley)

plot_knots              = [2,midpoints_surfplan, 'blue',4]
plot_pulleys            = [2,midpoints_surfplan_pulley, 'red',7]
plot_bridle_lines_conn  = [midpoints_surfplan,bridle_connectivity_i,bridle_connectivity_j,'black',1]
get_3D_plot_N_D_4times_PLUS_lines(plot_pulleys,plot_knots,plot_knots,plot_knots,plot_bridle_lines_conn)

#%% Defining the indices that correspond to pulleys

tol_proximity_bridle_vertices_is_pulley = 0.005 #m

pulley_idx_list = [] #initializing
# loop through each vertex point
for i in range(len(midpoints_surfplan)):
    for j in range(len(midpoints_surfplan_pulley)):
        if euclidian_distance(midpoints_surfplan[i],midpoints_surfplan_pulley[j]) < tol_proximity_bridle_vertices_is_pulley:
            pulley_idx_list.append(i)

## extracting the pulley points from the pulley_idx_list to check if it worked
midpoints_surfplan_pulley_from_idx = []
for i in pulley_idx_list:
    midpoints_surfplan_pulley_from_idx.append(midpoints_surfplan[i])

## plotting
plot_knots              = [2,midpoints_surfplan, 'blue',4]
plot_pulleys            = [2,midpoints_surfplan_pulley_from_idx, 'red',7]
plot_bridle_lines_conn  = [midpoints_surfplan,bridle_connectivity_i,bridle_connectivity_j,'black',1]
get_3D_plot_N_D_4times_PLUS_lines(plot_pulleys,plot_knots,plot_knots,plot_knots,plot_bridle_lines_conn)



#%% COMMENTED SECTION, should be used later on
##TODO: uncomment the below, should actually get rid of the knots that are not connected to anything


# saving the old configuration
bridle_connectivity_i_old = copy.deepcopy(bridle_connectivity_i)
bridle_connectivity_j_old = copy.deepcopy(bridle_connectivity_j)
midpoints_surfplan_old    = copy.deepcopy(midpoints_surfplan)
#%% Getting rid of knots, were no lines are splitted and only two lines are connected
###   
###   
###   ## find min y-coordinate
###   def sorting_condition_y(point):
###           return point[1]
###   vertices_bridle_min_y = min(vertices_bridle,key=sorting_condition_y)[1]
###   tol_vertices_bridle_min_y_to_be_KCU = 0.1 #m
###   condition_vertices_bridle_min_y_to_be_KCU = vertices_bridle_min_y + tol_vertices_bridle_min_y_to_be_KCU
### 
# 
# 
# 

##TODO: new attempt below here, also decided to not actually be usefull... 
# Reason being the non straight bridles near the KCU, there you can't remove the knots
# Just leave it as it, physically the knots shouldn't be a problem
#  
##########     # defining the np.arrays
##########     bridle_connectivity_i = np.array(copy.deepcopy(bridle_connectivity_i_old))
##########     bridle_connectivity_j = np.array(copy.deepcopy(bridle_connectivity_j_old))
##########     midpoints_surfplan    = copy.deepcopy(midpoints_surfplan_old)
##########     
##########     counter = 0
##########     counter_deleted_points = []
##########     #while len(midpoints_surfplan)-counter > 0:
##########     idx_false_knots_list = []
##########     
##########     for i in range(0,len(midpoints_surfplan)):
##########         # get the number of occurences of the current index
##########     #    index_occ_i = np.count_nonzero(bridle_connectivity_i == bridle_connectivity_i[counter])
##########     #    index_occ_j = np.count_nonzero(bridle_connectivity_j == bridle_connectivity_i[counter])
##########     #    index_occ    = index_occ_i + index_occ_j
##########     #    print(index_occ)
##########     
##########         index_occ_i = np.where(bridle_connectivity_i == i)[0]
##########         index_occ_j = np.where(bridle_connectivity_j == i)[0]
##########         index_occ   = np.hstack((index_occ_i,index_occ_j))
##########         
##########         if len(index_occ) == 2:
##########             idx_false_knots_list.append(index_occ)
##########     
##########             ##TODO: Add a comment to delete indices from bridle_conn_i and j
##########     
##########     
##########     
##########     ## plotting to check
##########     plot_pulleys            = [3,proximity_list_pulley, 'red',1]
##########     plot_knots              = [2,midpoints_surfplan, 'blue',3]
##########     
##########     ## making this lists into one
##########     idx_false_knots_list = np.concatenate(idx_false_knots_list,axis=0)
##########     
##########     ##TODO: remove duplicates?
##########     
##########     ## making indices into a list
##########     false_knot_list = []
##########     for i in idx_false_knots_list:
##########         false_knot_list.append(midpoints_surfplan[i])
##########     
##########     plot_false_knots = [2,false_knot_list, 'green',6]
##########     
##########     plot_bridle_lines_conn  = [midpoints_surfplan,bridle_connectivity_i,bridle_connectivity_j,'black',1]
##########     get_3D_plot_N_D_4times_PLUS_lines(plot_pulleys,plot_knots,plot_false_knots,plot_knots,plot_bridle_lines_conn)
##########     
##########     


###   
###       # first look for neighbouring points
###       i_idx_list = np.where(bridle_connectivity_i == counter)[0]
###       j_idx_list = np.where(bridle_connectivity_j == counter)[0]
###       index_occ = len(i_idx_list) + len(j_idx_list)
###       
###       # if wrong point
###       # (if the index_occ is not 2 and if it's above the KCU attached points)
###       y_coordinate = midpoints_surfplan[counter][1]
###       if index_occ == 2 and y_coordinate > condition_vertices_bridle_min_y_to_be_KCU:
###           
###           ## if both occurences are in the i_idx_list
###           if len(i_idx_list) == 2:
###               # then we must append both other occurences, the neighbours
###               ##TODO: something still goes wrong in the appending part
###   
###               np.append(bridle_connectivity_i,bridle_connectivity_j[i_idx_list[0]])
###               np.append(bridle_connectivity_j,bridle_connectivity_j[i_idx_list[1]])
###                           
###               # removing the points
###               bridle_connectivity_i   = np.delete(bridle_connectivity_i,i_idx_list[0])
###               bridle_connectivity_j   = np.delete(bridle_connectivity_j,i_idx_list[0])
###               bridle_connectivity_i   = np.delete(bridle_connectivity_i,i_idx_list[1]-1)
###               bridle_connectivity_j   = np.delete(bridle_connectivity_j,i_idx_list[1]-1)
###   
###           ## if both occurences are in the j_idx_list    
###           elif len(j_idx_list) == 2:
###               # then we must append both other occurences, the neighbours
###               np.append(bridle_connectivity_i,bridle_connectivity_i[j_idx_list[0]])
###               np.append(bridle_connectivity_j,bridle_connectivity_i[j_idx_list[1]])
###   
###               # removing the points
###               bridle_connectivity_i   = np.delete(bridle_connectivity_i,j_idx_list)
###               bridle_connectivity_j   = np.delete(bridle_connectivity_j,j_idx_list)
###   
###           ## if one occurence is in the i_idx_list and one in the j_idx_list
###           else:
###               # then we must append the opposite occurence, the neighbour
###               np.append(bridle_connectivity_i,bridle_connectivity_i[j_idx_list[0]])
###               np.append(bridle_connectivity_j,bridle_connectivity_j[i_idx_list[0]])
###               
###               # removing the points
###               #idx_removal_list = np.concatenate((i_idx_list,j_idx_list))
###               bridle_connectivity_i   = np.delete(bridle_connectivity_i,[i_idx_list[0],j_idx_list[0]])
###               bridle_connectivity_j   = np.delete(bridle_connectivity_j,[i_idx_list[0],j_idx_list[0]])
###   
###           # once appended we can delete the points
###           # Alternative method is needed to keep midpoints_surfplan a list!
###           #del midpoints_surfplan[counter]
###           #midpoints_surfplan      = np.delete(midpoints_surfplan,counter)
###   
###           # Storing the counter
###           counter_deleted_points.append(counter)
###   
###       # Updating the counter
###       counter += 1
###   

###   
###   
###   #%%
###   #get_3D_plot_with_lines(midpoints_surfplan,bridle_connectivity_i,bridle_connectivity_j,show_legend=False)
###   
###   #plot_removed_knots = [3,proximity_list, 'blue',1]
###   plot_knots         =    [3,proximity_list, 'blue',2]
###   plot_pulleys       = [3,proximity_list_pulley, 'red',2]
###   plot_bridle_lines_conn = [midpoints_surfplan,bridle_connectivity_i,bridle_connectivity_j,'black',1]
###   get_3D_plot_N_D_4times_PLUS_lines(plot_pulleys,plot_knots,plot_knots,plot_knots,plot_bridle_lines_conn)
###   #get_3D_plot(midpoints_surfplan, 'blue',1)

#%% Defining the supervector and the particles as indices
## The wing_supervector
# Ordered left-to-right and in the rib-count direction ([0] to [7])
# Thus 12 particles, 8 inherently wing particles and 4~ bridle line attachment particles

supervector_wing = [] # the supervector
idx_counter = 0 # the counter
indices_list = [] # the list of indices

for i in range(len(rib_db_whole)):
    ##TODO: should also do this for the MID ofcourse
    if rib_db_whole[i][0] == "strut" or "mid-tip":
        # creating supervector, by appending particles
        supervector_wing.append(rib_db_whole[i][4])

        idx_counter += len(rib_db_whole[i][4]) # updating the counter
        new_indices = np.arange(idx_counter-len(rib_db_whole[i][4]),idx_counter) # getting the new indices
        
        # ensuring there are only 5 lists ##TODO: only here to run this cell multiple times
        rib_db_whole[i] = rib_db_whole[i][:5]

        # appending a 6th element, which contains the indices of the particles
        rib_db_whole[i].append(new_indices)

        indices_list.append(new_indices) ##TODO: used for testing, remove later (test checks out)

# supervector wing is now a list of lists, we need to convert it to a numpy array
supervector_wing = np.concatenate(supervector_wing,axis=0)
N_wing_particles = int(len(supervector_wing))

plot_strut_LE = [3,rib_strut_LE, 'red',1]
plot_strut_canopy = [3,rib_strut_canopy, 'green',1]
plot_strut_TE = [3,rib_strut_TE, 'blue',1]
plot_strut_particles = [2,supervector_wing, 'black',3]

# same plot as last cell, but now with particles called as supervector
get_3D_plot_N_D_5times(plot_strut_TE,plot_strut_LE,plot_strut_canopy,plot_strut_particles,plot_mid_tip_LE)

rib_db_whole_old = copy.deepcopy(rib_db_whole) #making a copy in case I mess up

#%% Find those bridle points that are part of the ribs

##TODO: left-off here, need to fix the mid_tip_bridle_connection_idx "problem"
supervector_bridles = copy.deepcopy(midpoints_surfplan) #making a copy in case I mess up
rib_db_whole = copy.deepcopy(rib_db_whole_old) #making a copy in case I mess up

## For each rib strut we need to:
    ## 1. Define a line between particle [0] and particle 7, with idx [6]
    ## 2. Find the indices of the points that lie within a defined proximity of that line
    ##     - Must correct for the supervector ordering WING -> BRIDLES
    ## 3. Order the indices, based on the points their proximity to TE (particle 7)
    ## 4. Add the indices to the rib_db "indices" list

## DEFINING TOLERANCE
tol_bridle_points_off_line = .2 #m
tol_mid_tip_bridle_point_from_LE_particle = 0.1 #m

rib_bridle_indices_list = [] # storing the indices

for i in range(len(rib_db_whole)):
    if rib_db_whole[i][0] == "strut": ##TODO: should also do this for the MID ofcourse

        ## 1. Define a line between particle [0] and the last particle 7 with idx [6]
        rib_lower_surface_line = [rib_db_whole[i][4][0],rib_db_whole[i][4][N_strut_particles-1]]

        ## 2. Find the indices of the points that lie within a defined proximity of that line
        rib_bridle_indices = [] # storing the indices of the bridle_attachtment points

        for j,bridle_point in enumerate(supervector_bridles): ##TODO: below is hardcoded with 100 distance
            if np.linalg.norm(lineseg_dists(supervector_bridles[j],rib_db_whole[i][4][0],rib_db_whole[i][4][N_strut_particles-1])) < tol_bridle_points_off_line:
                ## Supervector is ordered WING -> Bridles. So we need to correct the indices for this
                rib_bridle_indices.append(int(j + N_wing_particles))

        ## 3. Order the indices, based on the points their proximity to TE (particle 7)
        def cond_dist_to_particle7(rib_bridle_idx):
            point = supervector_bridles[rib_bridle_idx-N_wing_particles] # again correcting for supervector ordering
            return euclidian_distance(point,rib_db_whole[i][4][N_strut_particles-1])

        rib_bridle_indices = sorted(rib_bridle_indices,key=cond_dist_to_particle7)

        ## 4. Add them to the rib_db "indices" list
        rib_db_whole[i][5] = np.concatenate((rib_db_whole[i][5],rib_bridle_indices),axis=0)

        rib_bridle_indices_list.append(rib_bridle_indices) ##TODO: used for test

    ## DEALING WITH THE MID-TIPS
    if rib_db_whole[i][0] == "mid-tip":
        
        mid_tip_bridle_connection_idx = 0 ##TODO: remove this, just here to make sure it's empty

        for j,bridle_point in enumerate(supervector_bridles):
            # if the bridle point is within proximity (tol_) of the first particle on the mid-tip LE
            if euclidian_distance(supervector_bridles[j],rib_db_whole[i][4][0]) < tol_mid_tip_bridle_point_from_LE_particle: 
                # then, we store the idx
                mid_tip_bridle_connection_idx = int(j + N_wing_particles)

        # appending to the rib_db list
        rib_db_whole[i][5] = np.concatenate((rib_db_whole[i][5],[mid_tip_bridle_connection_idx]),axis=0)
        
        # used for test-plotting
        rib_bridle_indices_list.append([mid_tip_bridle_connection_idx]) ##TODO: used for test

## testing if the above worked
rib_bridle_indices_list = np.concatenate(rib_bridle_indices_list,axis=0) #going to 1D array

rib_bridle_points_list = []
for i in np.arange(len(rib_bridle_indices_list)):
    idx = int(rib_bridle_indices_list[i]-N_wing_particles)
    point = supervector_bridles[idx]
    rib_bridle_points_list.append(point)

plot_strut_LE = [3,rib_strut_LE, 'red',1]
plot_strut_canopy = [3,rib_strut_canopy, 'green',1]
plot_strut_TE = [3,rib_strut_TE, 'blue',1]
plot_strut_particles = [2,supervector_wing, 'black',3]
plot_strut_particles_bridle = [2,rib_bridle_points_list, 'orange',5]

#get_3D_plot_N_D_4times(plot_strut_particles_bridle,plot_strut_LE,plot_strut_canopy,plot_strut_particles)
get_3D_plot_N_D_triple(plot_strut_particles_bridle,[2,supervector_bridles, 'blue',2],plot_strut_particles)

#%% WING CONNECTIVITY
# wing connectivity
## Defining the supervector
supervector = np.concatenate((supervector_wing,supervector_bridles),axis=0)

###### CONNECTIVITY ######
# Use an connectivity / adjacency matrix to define the connectivity between the particles
# This is a symmetric matrix of size (n_particles x n_particles)
# All the connections are 1's and the rest is 0's
# Which makes a sparse matrix definition attractive
# Choice is the scipy.sparse.coo_matrix (COOrdinate format)
# The input here is 3fold, all being 1D np.arrays:
    # 1) the row indices (from supervector, index of particle a, connected to b)
    # 2) the column indices (from supervector, index of particle b, connected to a)
    # 3) the values (1's in this case)
# The output can be written as sparse matrix or dense matrix

row_ind, col_ind, conn_indicator = [],[],[] #defining the empty lists
counter_mid_tip = 0 #counter for the mid-tip rib

for i in range(len(rib_db_whole)):
    if rib_db_whole[i][0] == "strut": ##TODO: should also do this for the MID ofcourse
        
        ########################################################
        ## inter-rib connectivity - I, the surface of the rib

        #adding the indices of the particles
        row_ind.append(rib_db_whole[i][5]) 
        # indices are ordered, so rows are just connected to cols in order
        col_ind.append(np.append(rib_db_whole[i][5][1:],rib_db_whole[i][5][0])) #last point is connected to first point
        # the values are all 1's
        conn_indicator.append(np.ones(len(rib_db_whole[i][5])))

        ########################################################
        ## inter-rib connectivity - II, the inner connections of the rib itself
        
        # connecting point 0 to 2 and point 4 to 6
        row_ind.append([rib_db_whole[i][5][0],rib_db_whole[i][5][4]])
        col_ind.append([rib_db_whole[i][5][2],rib_db_whole[i][5][6]])
        conn_indicator.append(np.ones(2))

        # Defining the other lines, coming from the bridle attachment points
        # Each bridle line point is connected to the 2 closest rib points that are non-neighbours
        #   - To remove the neighbours, [0] and [6] are remove
        #   - The 2 closest points are found by sorting, based on proximity
        # j loops over each bridle line attachment point (from [7] to [end])
        for j in range(N_strut_particles,len(rib_db_whole[i][5])):
            
            # defining the current index of the bridle line attachment point
            index_j = rib_db_whole[i][5][j]

            # Define a sorting condition, based on the distance to the current particle
            def cond_dist_to_current_particle(index):
                return euclidian_distance(supervector[index],supervector[index_j])
            
            # find possible connections, which should include all the points in the rib
            # except for the first [0] and point [6]
            list_of_possible_connections = rib_db_whole[i][5][1:int(N_strut_particles-1)]

            # sort the list based on proximity
            proximity_list = sorted(list_of_possible_connections,key=cond_dist_to_current_particle)

            row_ind.append([index_j,index_j])
            col_ind.append([proximity_list[0],proximity_list[1]])
            conn_indicator.append(np.ones(2))


        ########################################################
        ## inter-rib connectivity - III, the canopy surface

        # treat treat the rest {1,..,-2} normal (last-non-tip is also skipped on purpose)
        ##TODO: -4 is added to skip the last 4 tip-mid ribs
        if i != 0 and i != (len(rib_db_whole)-2 -4) and i < (len(rib_db_whole)-1-4):
            # a list is created that has the indices of the current rib
            # the list describes the points [last (11 in this case),2,3,4,5]
            row_ind.append(np.concatenate(([rib_db_whole[i][5][-1]],rib_db_whole[i][5][0:(N_strut_particles-1)]),axis=0))
            # the same points in the next STRUT rib (i+1) are connected
            col_ind.append(np.concatenate(([rib_db_whole[i+1][5][-1]],rib_db_whole[i+1][5][0:(N_strut_particles-1)]),axis=0))
            # the values are all 1's
            conn_indicator.append(np.ones(N_strut_particles)) ##TODO: this does not always work, but works for now


    if rib_db_whole[i][0] == "mid-tip": ##TODO: should also do this for ALL the MIDs ofcourse
        

        ########################################################
        ## inter-rib connectivity - I, the surface of the rib

        ##TODO: could add a 4th point on the LE to better resemble the full circular structure

        # Walking from the front to the back
        #rows go like; [ bridle attachment point, LE1, LE2, canopy1 ... , TE]    
        row_ind.append(np.concatenate(([rib_db_whole[i][5][-1]],rib_db_whole[i][5][:-2]),axis=0) ) #(itself) adding the indices of the particles, except for the last one
        #cols go like; [ LE1, LE2, canopy1 ... , TE, bridle attachment point]
        col_ind.append(rib_db_whole[i][5][:-1]) #(itself)
        conn_indicator.append(np.ones(len(rib_db_whole[i][5][:-1])))

        ########################################################
        ## inter-rib connectivity - II, the inner connections of the rib itself
        
        # connecting bridle point to LE2
        row_ind.append([rib_db_whole[i][5][-1]])
        col_ind.append([rib_db_whole[i][5][1]])
        conn_indicator.append(np.ones(1))


        ########################################################
        ## inter-rib connectivity - III, the canopy surface

        if counter_mid_tip == 0: #most left mid-tip
            
            #from most left-tip to most left mid-tip
            row_ind.append(rib_db_whole[0][5][:6]) ##TODO: hardcoded
            col_ind.append(rib_db_whole[i][5][:6])
            conn_indicator.append(np.ones(6))

        elif counter_mid_tip == 1: #2nd most left mid-tip

            #first connection, from 1st most left to 2nd most left ##TODO: only works because they are in order
            row_ind.append(np.concatenate(([rib_db_whole[i-1][5][-1]],rib_db_whole[i-1][5][:6]),axis=0)) ##TODO: hardcoded
            col_ind.append(np.concatenate(([rib_db_whole[i][5][-1]],rib_db_whole[i][5][:6]),axis=0))
            conn_indicator.append(np.ones(7))

            #second connection, from 2nd most left to first (non-tip) strut
            row_ind.append(np.concatenate(([rib_db_whole[i][5][-1]],rib_db_whole[i][5][:6]),axis=0)) ##TODO: hardcoded
            col_ind.append(np.concatenate(([rib_db_whole[1][5][-1]],rib_db_whole[1][5][:6]),axis=0))
            conn_indicator.append(np.ones(7))

        elif counter_mid_tip == 2: #3rd most left mid-tip

            #first connection, from  first (non-tip) strut to 2nd most right
            row_ind.append(np.concatenate(([rib_db_whole[12][5][-1]],rib_db_whole[12][5][:6]),axis=0)) ##TODO: hardcoded
            col_ind.append(np.concatenate(([rib_db_whole[i][5][-1]],rib_db_whole[i][5][:6]),axis=0))
            conn_indicator.append(np.ones(7))

            #second connection, from 2nd most right to most right mid-tip ##TODO: only works because they are in order
            row_ind.append(np.concatenate(([rib_db_whole[i][5][-1]],rib_db_whole[i][5][:6]),axis=0)) ##TODO: hardcoded
            col_ind.append(np.concatenate(([rib_db_whole[i+1][5][-1]],rib_db_whole[i+1][5][:6]),axis=0))
            conn_indicator.append(np.ones(7))

        elif counter_mid_tip == 3: #most right mid-tip
            
            # from most right mid-tip to right tip
            row_ind.append(rib_db_whole[i][5][:6]) ##TODO: hardcoded
            col_ind.append(rib_db_whole[13][5][:6])
            conn_indicator.append(np.ones(6))
        
        counter_mid_tip += 1


# converting to 1D arrays
row_ind = np.concatenate(row_ind,axis=0).tolist()
col_ind = np.concatenate(col_ind,axis=0).tolist()
conn_indicator = np.concatenate(conn_indicator,axis=0).tolist()

# plotting the wing connectivity
plot_line_varibles = [supervector,row_ind,col_ind, 'blue',2]
get_3D_plot_N_D_4times_PLUS_lines(plot_strut_particles_bridle,plot_strut_LE,plot_mid_tip_LE,plot_strut_particles,plot_line_varibles)

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
# With the supervector and connectivities defined
# We can rotate the axis and generate the desired output
# points_ini, ci,cj and ci_bridle, cj_bridle
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################


#%% Generating desired output
## Generating the desired output

## Rotating axis
points_surfplan = copy.deepcopy(supervector)
points_ini_model = transform_coordinate_system_Surfplan_to_Model(points_surfplan)

## Redefining the terms, to match naming convention
ci_wing = copy.deepcopy(row_ind)
cj_wing = copy.deepcopy(col_ind)
ci_bridle = copy.deepcopy(bridle_connectivity_i)
cj_bridle = copy.deepcopy(bridle_connectivity_j)

## Correcting bridle_connectivity as its placed aft of the wing_connectivity points
for i,ci_bridle_index in enumerate(ci_bridle):
    ci_bridle[i] = ci_bridle_index + len(supervector_wing)
for i,cj_bridle_index in enumerate(cj_bridle):
    cj_bridle[i] = cj_bridle_index + len(supervector_wing)

## Checking output
plot_supervector = [2,supervector,'black',3]
plot_connectivity = [supervector,np.concatenate((ci_wing,ci_bridle),axis=0),np.concatenate((cj_wing,cj_bridle),axis=0), 'blue',2]
#get_3D_plot_N_D_4times_PLUS_lines(plot_supervector,plot_strut_LE,plot_mid_tip_LE,plot_strut_particles_bridle,plot_line_varibles)
#get_3D_plot_with_lines(points_ini_model,plot_connectivity[1],plot_connectivity[2],show_legend=False)

## wing
# LE circles
plot_mid_tip_LE = [3, plot_mid_tip_LE[1], 'red',.1]
plot_strut_LE = [3,rib_strut_LE, 'red',.1]

# canopy
plot_strut_particles = [2,supervector_wing, 'blue',5]

# bridle line connection points
plot_strut_particles_bridle = [2,rib_bridle_points_list, 'orange',5]

# connections
plot_line_varibles_wing = [supervector,ci_wing,cj_wing, 'blue',2]

## bridles
plot_bridle_points = [2,supervector_bridles, 'black',3]
plot_line_varibles_bridle = [supervector,ci_bridle,cj_bridle, 'black',1]

#pulleys
supervector_bridles_pulleys = [] ## extracting the pulley points from the pulley_idx_list to check if it worked
for i in pulley_idx_list:
    supervector_bridles_pulleys.append(supervector_bridles[i])
plot_pulley_points = [2,supervector_bridles_pulleys, 'red',5]

get_3D_plot_N_D_6times_PLUS_2lines(plot_strut_particles_bridle,plot_strut_LE,plot_mid_tip_LE,plot_strut_particles,plot_bridle_points,plot_pulley_points,plot_line_varibles_wing,plot_line_varibles_bridle)

##TODO: left-off here
##TODO: left-off here

#------------- Realize that you are in v4 and in v5 you have some improvements
#------------- Improvements towards finding the mid-surface mainly

##TODO: left-off here
##TODO: left-off here

#%% ROTATING ALL THE OUTPUT DATA
# rotating all output data

# all particles
points_surfplan = copy.deepcopy(supervector)
points_ini_model = transform_coordinate_system_Surfplan_to_Model(points_surfplan)

## WING
# wing vertices
points_surfplan_wing_vertices = copy.deepcopy(vertices_wing)
points_ini_model_wing_vertices = transform_coordinate_system_Surfplan_to_Model(points_surfplan_wing_vertices)

# wing particles
points_surfplan_wing = copy.deepcopy(supervector_wing)
points_ini_model_wing = transform_coordinate_system_Surfplan_to_Model(points_surfplan_wing)

# strut LE particles
points_surfplan_LE_per_strut = copy.deepcopy(rib_strut_LE)
points_ini_model_LE_per_strut = [transform_coordinate_system_Surfplan_to_Model(points_surfplan_LE_per_strut[i]) for i in range(len(points_surfplan_LE_per_strut))]

# mid-tip LE particles
points_surfplan_LE_per_mid = copy.deepcopy(plot_mid_tip_LE[1])
points_ini_model_LE_per_mid_tip = [transform_coordinate_system_Surfplan_to_Model(points_surfplan_LE_per_mid[i]) for i in range(len(points_surfplan_LE_per_mid))]

# canopy particles
points_surfplan_canopy_per_strut = copy.deepcopy(rib_strut_canopy)
points_ini_model_canopy_per_strut = [transform_coordinate_system_Surfplan_to_Model(points_surfplan_canopy_per_strut[i]) for i in range(len(points_surfplan_canopy_per_strut))]

# TE particles
points_surfplan_TE_per_strut = copy.deepcopy(rib_strut_TE)
points_ini_model_TE_per_strut = [transform_coordinate_system_Surfplan_to_Model(points_surfplan_TE_per_strut[i]) for i in range(len(points_surfplan_TE_per_strut))]

# bridle line attachment points
points_surfplan_bridle_attachment_points = copy.deepcopy(rib_bridle_points_list)
points_ini_model_bridle_attachment_points = transform_coordinate_system_Surfplan_to_Model(points_surfplan_bridle_attachment_points)


## BRIDLES
# bridle vertices
points_surfplan_bridle_vertices = copy.deepcopy(vertices_bridle)
points_ini_model_bridle_vertices = transform_coordinate_system_Surfplan_to_Model(points_surfplan_bridle_vertices)

# bridles particles
points_surfplan_bridle = copy.deepcopy(supervector_bridles)
points_ini_model_bridle = transform_coordinate_system_Surfplan_to_Model(points_surfplan_bridle)

# pulley particles
supervector_bridles_pulleys = [] ## extracting the pulley points from the pulley_idx_list to check if it worked
for i in pulley_idx_list:
    supervector_bridles_pulleys.append(supervector_bridles[i])

points_surfplan_pulleys = copy.deepcopy(supervector_bridles_pulleys)
points_ini_model_pulleys = transform_coordinate_system_Surfplan_to_Model(points_surfplan_pulleys)

## rib_db_whole
rib_db_whole_surfplan = copy.deepcopy(rib_db_whole)
rib_db_whole_model = rib_db_whole_surfplan
for i in range(len(rib_db_whole_surfplan)):
    rib_db_whole_model[i][1] = transform_coordinate_system_Surfplan_to_Model(rib_db_whole_surfplan[i][1])
    rib_db_whole_model[i][2] = transform_coordinate_system_Surfplan_to_Model(rib_db_whole_surfplan[i][2])
    rib_db_whole_model[i][3] = transform_coordinate_system_Surfplan_to_Model(rib_db_whole_surfplan[i][3])
    rib_db_whole_model[i][4] = transform_coordinate_system_Surfplan_to_Model(rib_db_whole_surfplan[i][4])

#%% Preparing the connectivity output
##% CONNECTIVITY

## Redefining the terms, to match naming convention
ci_bridle_indexed_on_bridle = copy.deepcopy(bridle_connectivity_i)
cj_bridle_indexed_on_bridle = copy.deepcopy(bridle_connectivity_j)

## Redefining the terms, to match naming convention
ci_wing = copy.deepcopy(row_ind)
cj_wing = copy.deepcopy(col_ind)
ci_bridle = copy.deepcopy(bridle_connectivity_i)
cj_bridle = copy.deepcopy(bridle_connectivity_j)

## Correcting bridle_connectivity as its placed aft of the wing_connectivity points
for i,ci_bridle_index in enumerate(ci_bridle):
    ci_bridle[i] = ci_bridle_index + len(supervector_wing)
for i,cj_bridle_index in enumerate(cj_bridle):
    cj_bridle[i] = cj_bridle_index + len(supervector_wing)


#%% CHECKING THE OUTPUT - VERTICES VS PARTICLES
## CHECKING THE OUTPUT - VERTICES VS PARTICLES

plot_points_ini_model = [2,points_ini_model,'blue',2]
plot_conn_model       = [points_ini_model,np.concatenate((ci_wing,ci_bridle),axis=0),np.concatenate((cj_wing,cj_bridle),axis=0), 'blue',2]

plot_points_ini_model_wing_vertices = [2,points_ini_model_wing_vertices,'black',1]
plot_points_ini_model_bridle_vertices = [2,points_ini_model_bridle_vertices,'black',1]

plot_points_ini_model_pulleys = [2,points_ini_model_pulleys,'red',3]

get_3D_plot_N_D_4times_PLUS_lines(plot_points_ini_model,plot_points_ini_model_wing_vertices,plot_points_ini_model_bridle_vertices,plot_points_ini_model_pulleys,plot_conn_model)


#%% CHECKING THE OUTPUT - WING DETAIL AND BRIDLE
## CHECKING THE OUTPUT - WING DETAIL AND BRIDLE DETAIL

plot_points_ini_model_wing = [2,points_ini_model_wing,'blue',2]
##TODO: realize that the WING CONNECTIVITY, is made using the bridle line attachment points 
## Therefore it requires all points to plot the lines, as it also needs the bridle line attachmentpoint points
plot_conn_wing       = [points_ini_model,ci_wing,cj_wing, 'blue',1]

plot_point_ini_model_strut_LE = [3,points_ini_model_LE_per_strut,'red',2]
plot_point_ini_model_mid_tip_LE = [3,points_ini_model_LE_per_mid_tip,'red',6]
plot_point_ini_model_strut_canopy = [3,points_ini_model_canopy_per_strut,'green',3]
plot_point_ini_model_strut_TE = [3,points_ini_model_TE_per_strut,'purple',2]
plot_points_ini_model_bridle_attachment_points = [2,points_ini_model_bridle_attachment_points,'orange',4.5]

plot_points_ini_model_bridle = [2,points_ini_model_bridle,'black',1]
plot_conn_bridle = [points_ini_model_bridle,ci_bridle_indexed_on_bridle,cj_bridle_indexed_on_bridle, 'black',1]

get_3D_plot_N_D_6times_PLUS_2lines(plot_points_ini_model_bridle,plot_point_ini_model_strut_LE,plot_point_ini_model_mid_tip_LE,plot_point_ini_model_strut_canopy,plot_point_ini_model_strut_TE,plot_points_ini_model_bridle_attachment_points,plot_conn_model,plot_conn_bridle)


##TODO: left-off here
#------------- Realize that you are in v4 and in v5 you have some improvements
#------------- Improvements towards finding the mid-surface mainly

#%% Storing 
### STORING 
#filename = '../AerostructuralModel_v'+date+'/run_results/pos_up_'+str(int(100*u_p))+'_'+str(delta_ld_used)+'_billowing.csv'

filepath = "/home/jellepoland/surfdrive/phd/code/phd_Jelle_Poland/Simulations/initialisation/V9_60C/results/"

## points
# wing
np.save(filepath+'vertices_wing_model.npy',transform_coordinate_system_Surfplan_to_Model(vertices_wing))
np.save(filepath+'points_ini_model.npy',points_ini_model)
np.save(filepath+'points_ini_model_wing.npy',points_ini_model_wing)
np.save(filepath+'points_ini_model_wing_vertices.npy',points_ini_model_wing_vertices)
np.save(filepath+'points_ini_model_LE_per_strut.npy',points_ini_model_LE_per_strut)
np.save(filepath+'points_ini_model_LE_per_mid_tip.npy',points_ini_model_LE_per_mid_tip)
np.save(filepath+'points_ini_model_canopy_per_strut.npy',points_ini_model_canopy_per_strut)
np.save(filepath+'points_ini_model_TE_per_strut.npy',points_ini_model_TE_per_strut)
np.save(filepath+'points_ini_model_bridle_attachment_points.npy',points_ini_model_bridle_attachment_points)

# bridle line system
np.save(filepath+'points_ini_model_bridle_vertices.npy',points_ini_model_bridle_vertices)
np.save(filepath+'points_ini_model_bridle.npy',points_ini_model_bridle)
np.save(filepath+'points_ini_model_pulleys.npy',points_ini_model_pulleys)
np.save(filepath+'pulley_point_indices.npy',pulley_idx_list)

# rib_db_whole
np.save(filepath+'rib_db_whole_model.npy',rib_db_whole_model)

## connectivity
np.save(filepath+'ci_bridle_indexed_on_bridle.npy',ci_bridle_indexed_on_bridle)
np.save(filepath+'cj_bridle_indexed_on_bridle.npy',cj_bridle_indexed_on_bridle)
np.save(filepath+'ci_wing.npy',ci_wing)
np.save(filepath+'cj_wing.npy',cj_wing)
np.save(filepath+'ci_bridle.npy',ci_bridle)
np.save(filepath+'cj_bridle.npy',cj_bridle)



#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################





