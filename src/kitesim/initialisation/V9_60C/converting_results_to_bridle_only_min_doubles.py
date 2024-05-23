#%% IMPORTING
import numpy as np
import copy

from kitesim.post_processing import functions_plot

## making things autoreload
%load_ext autoreload
%autoreload 2

filename = 'V9_60C_4_53_M_lines'
filepath = 'src/initialisation/V9_60C/results/'


## points
# wing
points_ini_model                            =np.load(filepath+'points_ini_model.npy',allow_pickle=True)
points_ini_model_wing                       =np.load(filepath+'points_ini_model_wing.npy',allow_pickle=True) 
points_ini_model_wing_vertices              =np.load(filepath+'points_ini_model_wing_vertices.npy',allow_pickle=True)
points_ini_model_LE_per_strut               =np.load(filepath+'points_ini_model_LE_per_strut.npy',allow_pickle=True)
points_ini_model_LE_per_mid_tip             =np.load(filepath+'points_ini_model_LE_per_mid_tip.npy',allow_pickle=True)
points_ini_model_canopy_per_strut           =np.load(filepath+'points_ini_model_canopy_per_strut.npy',allow_pickle=True)
points_ini_model_TE_per_strut               =np.load(filepath+'points_ini_model_TE_per_strut.npy',allow_pickle=True)
points_ini_model_bridle_attachment_points   =np.load(filepath+'points_ini_model_bridle_attachment_points.npy',allow_pickle=True)

# bridles
points_ini_model_bridle_vertices    =np.load(filepath+'points_ini_model_bridle_vertices.npy',allow_pickle=True)
points_ini_model_bridle             =np.load(filepath+'points_ini_model_bridle.npy',allow_pickle=True)
points_ini_model_pulleys            =np.load(filepath+'points_ini_model_pulleys.npy',allow_pickle=True)
pulley_point_indices                =np.load(filepath+'pulley_point_indices.npy',allow_pickle=True)

# rib_db_whole
rib_db_whole_model  = np.load(filepath+'rib_db_whole_model.npy',allow_pickle=True)

## connectivity
ci_bridle_indexed_on_bridle     =np.load(filepath+'ci_bridle_indexed_on_bridle.npy',allow_pickle=True)
cj_bridle_indexed_on_bridle     =np.load(filepath+'cj_bridle_indexed_on_bridle.npy',allow_pickle=True)
ci_wing                         =np.load(filepath+'ci_wing.npy',allow_pickle=True)
cj_wing                         =np.load(filepath+'cj_wing.npy',allow_pickle=True)
ci_bridle                       =np.load(filepath+'ci_bridle.npy',allow_pickle=True)
cj_bridle                       =np.load(filepath+'cj_bridle.npy',allow_pickle=True)

# #%% CHECKING THE INPUT - VERTICES VS PARTICLES
## CHECKING THE INPUT - VERTICES VS PARTICLES

plot_points_ini_model = [2,points_ini_model,'blue',2]
plot_conn_model       = [points_ini_model,np.concatenate((ci_wing,ci_bridle),axis=0),np.concatenate((cj_wing,cj_bridle),axis=0), 'blue',2]

plot_points_ini_model_wing_vertices = [2,points_ini_model_wing_vertices,'black',1]
plot_points_ini_model_bridle_vertices = [2,points_ini_model_bridle_vertices,'red',1]

plot_points_ini_model_pulleys = [2,points_ini_model_pulleys,'red',3]

# get_3D_plot_N_D_4times_PLUS_lines(plot_points_ini_model,plot_points_ini_model_wing_vertices,plot_points_ini_model_bridle_vertices,plot_points_ini_model_pulleys,plot_conn_model)


# #%% CHECKING THE OUTPUT - WING DETAIL AND BRIDLE
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

# get_3D_plot_N_D_6times_PLUS_2lines(plot_points_ini_model_bridle,plot_point_ini_model_strut_LE,plot_point_ini_model_mid_tip_LE,plot_point_ini_model_strut_canopy,plot_point_ini_model_strut_TE,plot_points_ini_model_bridle_attachment_points,plot_conn_model,plot_conn_bridle)


# #%% MAKING A SIMPLER VERSION
## making a simpler version
###########################################
###########################################
###########################################
###########################################
###########################################
###########################################

# loop through each bridle point list
# check if there is a ci_wing and cj_wing involved
# store those

## find those wing connections that connect bridle points to one another
# as it turns out this quitenicely makes up the tubular frame, so we name it after that
conn_tubular_frame_i,conn_tubular_frame_j = [],[]

# first one loops through the connectivity of the wing
for i in np.arange(0,len(ci_wing)):

    # then check if the connection is a bridle attachment point 
    # and directly also if it connects to another
    if points_ini_model[ci_wing[i]] in points_ini_model_bridle_attachment_points and points_ini_model[cj_wing[i]] in points_ini_model_bridle_attachment_points:
        
        # if it does, store both the i and j connectivities
        conn_tubular_frame_i.append(ci_wing[i])
        conn_tubular_frame_j.append(cj_wing[i])

## plotting to check if it works
plot_conn_tubular_frame = [points_ini_model,conn_tubular_frame_i,conn_tubular_frame_j, 'grey',10]
plot_outer_points = [2,points_ini_model_bridle_attachment_points[:5], 'black',8]
# get_3D_plot_N_D_4times_PLUS_lines(plot_outer_points,plot_points_ini_model_bridle_attachment_points,plot_points_ini_model_bridle,plot_points_ini_model_pulleys,plot_conn_tubular_frame)

# #%% Order is not correct, so we must order it, angle wise. 
## left should be red, and right should be blue

# copying to not overwrite the original
points_ini_model_bridle_attachment_points_sorted = copy.deepcopy(points_ini_model_bridle_attachment_points)

## Order is found to not be correct, so we must reorder it from left-to-right
# first must define the sorting condition
def condition_sorting_left_to_right(point): #defining the sorting condition ##TODO: this can also be written as a lambda function
    #min_z = min(point_i[2] for point_i in points_ini_model_wing)
    min_z = 0
    angle = np.arctan( (point[1]) / (point[2]-min_z) )
    return angle

## sorting on the set angle condition from tip to mid-span
#points_ini_model_bridle_attachment_points.sort(key=condition_sorting_left_to_right)
points_ini_model_bridle_attachment_points_sorted = sorted(points_ini_model_bridle_attachment_points_sorted,key=condition_sorting_left_to_right)

plot_outer_points_first3 = [2,points_ini_model_bridle_attachment_points_sorted[:3], 'red',8]
plot_outer_points_last3 = [2,points_ini_model_bridle_attachment_points_sorted[-3:], 'blue',8]

# get_3D_plot_N_D_4times_PLUS_lines(plot_outer_points_first3,plot_outer_points_last3,plot_points_ini_model_bridle_attachment_points,plot_points_ini_model_bridle,plot_conn_tubular_frame)

# #%% Now we can append the left and right lines to the connectivity
# append the left line to connectivity
conn_tubular_frame_i = np.append(conn_tubular_frame_i,np.where(points_ini_model==points_ini_model_bridle_attachment_points_sorted[1])[0])
conn_tubular_frame_j = np.append(conn_tubular_frame_j,np.where(points_ini_model==points_ini_model_bridle_attachment_points_sorted[2])[0])

# append the right line to connectivity
conn_tubular_frame_i = np.append(conn_tubular_frame_i,np.where(points_ini_model==points_ini_model_bridle_attachment_points_sorted[-3])[0])
conn_tubular_frame_j = np.append(conn_tubular_frame_j,np.where(points_ini_model==points_ini_model_bridle_attachment_points_sorted[-2])[0])

plot_conn_tubular_frame = [points_ini_model,conn_tubular_frame_i,conn_tubular_frame_j, 'grey',10]
# get_3D_plot_N_D_4times_PLUS_lines(plot_outer_points_first3,plot_outer_points_last3,plot_points_ini_model_bridle_attachment_points,plot_points_ini_model_bridle,plot_conn_tubular_frame)


# #%% Finding the membrane connectivity
## red shows the canopy springs, grey the tubular springs

## first define a list that only contains struts
rib_db_whole_model_struts = []
for i in np.arange(0,len(rib_db_whole_model)):
    # if not middle but strut
    if rib_db_whole_model[i][0]=="strut":
        rib_db_whole_model_struts.append(rib_db_whole_model[i])

# loop through ribs, if not first or last
# then make a connection line between the things

conn_canopy_i,conn_canopy_j = [],[]
# loop through ribs (arrange left-to-right as [rib,mid,rib,mid,..] )
for i in np.arange(0,len(rib_db_whole_model_struts)-1):

    # The rear (TE) line first
    conn_canopy_i.append(rib_db_whole_model_struts[i]  [5][7])
    conn_canopy_j.append(rib_db_whole_model_struts[i+1][5][7])

    # the next-up line
    conn_canopy_i.append(rib_db_whole_model_struts[i]  [5][8])
    conn_canopy_j.append(rib_db_whole_model_struts[i+1][5][8])

    if i > 0 and i < (len(rib_db_whole_model_struts)-2):
        # the third-up line
        conn_canopy_i.append(rib_db_whole_model_struts[i]  [5][9])
        conn_canopy_j.append(rib_db_whole_model_struts[i+1][5][9])

        # the fourth-up line
        conn_canopy_i.append(rib_db_whole_model_struts[i]  [5][10])
        conn_canopy_j.append(rib_db_whole_model_struts[i+1][5][10])

        # constraining the additional points present on these struts
        if i < 3 or i > 9:
            # 11
            conn_canopy_i.append(rib_db_whole_model_struts[i]  [5][11])
            conn_canopy_j.append(rib_db_whole_model_struts[i+1][5][11])


## still need to do something with the most outer strut attached points

## append the left two lines to connectivity
# more rear one
##TODO: the -1 is a bit of a hack, because the np.where finds multiple solutions, it works however.. 
conn_canopy_i.append(np.where(points_ini_model==points_ini_model_bridle_attachment_points_sorted[2])[0][-1])
conn_canopy_j.append(rib_db_whole_model_struts[1][5][9])
# more front one
conn_canopy_i.append(np.where(points_ini_model==points_ini_model_bridle_attachment_points_sorted[3])[0][-1])
conn_canopy_j.append(rib_db_whole_model_struts[1][5][10])


## append the right two lines to connectivity
# more rear one
conn_canopy_i.append(np.where(points_ini_model==points_ini_model_bridle_attachment_points_sorted[-3])[0][0])
conn_canopy_j.append(rib_db_whole_model_struts[-2][5][9])
# more front one
conn_canopy_i.append(np.where(points_ini_model==points_ini_model_bridle_attachment_points_sorted[-4])[0][0])
conn_canopy_j.append(rib_db_whole_model_struts[-2][5][10])

## plotting results
plot_points_ini_model_bridle_attachment_points = [2,points_ini_model_bridle_attachment_points, 'orange',4]
plot_points_ini_model_canopy_per_strut = [3,points_ini_model_canopy_per_strut, 'blue',1]
plot_points_ini_model_TE_per_strut = [3,points_ini_model_TE_per_strut, 'blue',1]
plot_points_ini_model_LE_per_strut = [3,points_ini_model_LE_per_strut, 'blue',1]
plot_points_ini_model_wing = [2,points_ini_model_wing, 'black',3]

plot_conn_canopy = [points_ini_model,conn_canopy_i,conn_canopy_j, 'red',5]
plot_conn_tubular_frame = [points_ini_model,conn_tubular_frame_i,conn_tubular_frame_j, 'grey',4]

#  get_3D_plot_N_D_6times_PLUS_2lines(plot_points_ini_model_bridle_attachment_points,plot_points_ini_model_wing,plot_points_ini_model_LE_per_strut,plot_points_ini_model_canopy_per_strut,plot_points_ini_model_wing,plot_points_ini_model_TE_per_strut,plot_conn_canopy,plot_conn_tubular_frame)


# #%% Combining the connectivity

conn_wing_i = np.concatenate((conn_tubular_frame_i,conn_canopy_i),axis=0)
conn_wing_j = np.concatenate((conn_tubular_frame_j,conn_canopy_j),axis=0)

plot_line_new_wing_conn = [points_ini_model,conn_wing_i,conn_wing_j, 'green',5]

# get_3D_plot_N_D_6times_PLUS_2lines(plot_points_ini_model_bridle_attachment_points,plot_points_ini_model_wing,plot_points_ini_model_LE_per_strut,plot_points_ini_model_LE_per_strut,plot_points_ini_model_LE_per_strut,plot_points_ini_model_TE_per_strut,plot_line_new_wing_conn,plot_conn_wing)


# #%% Defining the special tip-plates
## defining the special tip-plates

# first defining the left-plate
tip_plate_left_LE_left  = np.where(points_ini_model==points_ini_model_bridle_attachment_points_sorted[2])[0][-1]
tip_plate_left_LE_right = rib_db_whole_model_struts[1][5][11]
tip_plate_left_TE_right = rib_db_whole_model_struts[1][5][7]
tip_plate_left_TE_left  = rib_db_whole_model_struts[0][5][7]
tip_plate_left = [tip_plate_left_LE_left,tip_plate_left_LE_right,tip_plate_left_TE_right,tip_plate_left_TE_left]

# second, defining the right-plate
tip_plate_right_LE_left  = rib_db_whole_model_struts[-2][5][11]
plate_rigth_LE_right = np.where(points_ini_model==points_ini_model_bridle_attachment_points_sorted[-3])[0][0]
tip_plate_right_TE_right = rib_db_whole_model_struts[-1][5][7]
tip_plate_right_TE_left  = rib_db_whole_model_struts[-2][5][7]
tip_plate_right = [tip_plate_right_LE_left,plate_rigth_LE_right,tip_plate_right_TE_right,tip_plate_right_TE_left]

tip_plate_left_points = [points_ini_model[idx] for idx in tip_plate_left] #[points_ini[tip_plate_left[0]],points_ini[tip_plate_left[1]],points_ini[tip_plate_left[2]],points_ini[tip_plate_left[3]]]
tip_plate_right_points = [points_ini_model[idx] for idx in tip_plate_right] # [points_ini[tip_plate_right[0]],points_ini[tip_plate_right[1]],points_ini[tip_plate_right[2]],points_ini[tip_plate_right[3]]]
plot_tip_plate_left  = [2,tip_plate_left_points, 'red',10]
plot_tip_plate_right = [2,tip_plate_right_points, 'blue',10]

## plotting to check if this worked out
plot_points_ini_model_bridle_attachment_points = [2,points_ini_model_bridle_attachment_points, 'orange',4]
#  get_3D_plot_N_D_6times_PLUS_2lines(plot_tip_plate_left,plot_tip_plate_right,plot_points_ini_model_bridle_attachment_points,plot_points_ini_model_bridle_attachment_points,plot_points_ini_model_bridle_attachment_points,plot_points_ini_model_TE_per_strut,plot_line_new_wing_conn,plot_conn_wing)


# #%% Defining the rest of the plates

plate_point_indices = [tip_plate_left]
for i in np.arange(1,len(rib_db_whole_model_struts)-2):
    
    TE_idx = 7
    LE_idx_i = 11
    LE_idx_i1 = 11

    # change the LE_idx, for two special struts at i = 2,11
    # the notation is a bit weird because we start at i = 1
    # and we also have to consider i+1
    if i == 2 or i == 11:
        LE_idx_i = 12
    elif (i+1) == 2 or (i+1) == 11:
        LE_idx_i1  = 12

    plate_LE_left_i     = rib_db_whole_model_struts[i]  [5][LE_idx_i]
    plate_LE_right_i    = rib_db_whole_model_struts[i+1][5][LE_idx_i1]
    plate_TE_right_i    = rib_db_whole_model_struts[i+1][5][TE_idx]
    plate_TE_left_i     = rib_db_whole_model_struts[i]  [5][TE_idx]

    plate_point_indices.append([plate_LE_left_i,plate_LE_right_i,plate_TE_right_i,plate_TE_left_i])

plate_point_indices.append(tip_plate_right)

## plotting to check if this worked out
plate_points = [points_ini_model[idx] for idx in plate_point_indices]
plot_plate_points = [3,plate_points, 'black',7]
# get_3D_plot_N_D_6times_PLUS_2lines(plot_tip_plate_left,plot_tip_plate_right,plot_points_ini_model_bridle_attachment_points,plot_plate_points,plot_points_ini_model_wing,plot_points_ini_model_TE_per_strut,plot_line_new_wing_conn,plot_conn_wing)


# #%% extracting the right set of starting points

## this also means that the indices must be adjusted
# format: "for idx,value in enumerate(list)"
# essentially what's happening as that the normal points are wing,bridle
# and the new points are only the bridle points, so those indices must be adjusted

for i in enumerate(conn_tubular_frame_i):
    conn_tubular_frame_i[i[0]] = conn_tubular_frame_i[i[0]] - len(points_ini_model_wing)
    conn_tubular_frame_j[i[0]] = conn_tubular_frame_j[i[0]] - len(points_ini_model_wing)

for i in enumerate(conn_canopy_i):
    conn_canopy_i[i[0]] = conn_canopy_i[i[0]] - len(points_ini_model_wing)
    conn_canopy_j[i[0]] = conn_canopy_j[i[0]] - len(points_ini_model_wing)

for i in enumerate(conn_wing_i):
    conn_wing_i[i[0]] = conn_wing_i[i[0]] - len(points_ini_model_wing)
    conn_wing_j[i[0]] = conn_wing_j[i[0]] - len(points_ini_model_wing)

for i in enumerate(plate_point_indices):
    for j in enumerate(plate_point_indices[i[0]]):
        plate_point_indices[i[0]][j[0]] = plate_point_indices[i[0]][j[0]] - len(points_ini_model_wing)

## this must also be done for the indices of the rib_db_whole_model_struts
for i,rib in enumerate(rib_db_whole_model_struts):
    for j,rib in enumerate(rib_db_whole_model_struts[i][5]):
        rib_db_whole_model_struts[i][5][j] = rib_db_whole_model_struts[i][5][j] - len(points_ini_model_wing)

## plotting results
plot_points_ini_model_bridle_attachment_points = [2,points_ini_model_bridle_attachment_points, 'orange',4]
plot_points_ini_model_canopy_per_strut = [3,points_ini_model_canopy_per_strut, 'blue',1]
plot_points_ini_model_TE_per_strut = [3,points_ini_model_TE_per_strut, 'blue',1]
plot_points_ini_model_LE_per_strut = [3,points_ini_model_LE_per_strut, 'blue',1]
plot_points_ini_model_wing = [2,points_ini_model_wing, 'black',3]

plot_conn_canopy = [points_ini_model_bridle,conn_canopy_i,conn_canopy_j, 'red',5]
plot_conn_tubular_frame = [points_ini_model_bridle,conn_tubular_frame_i,conn_tubular_frame_j, 'grey',3]

# get_3D_plot_N_D_6times_PLUS_2lines(plot_points_ini_model_bridle_attachment_points,plot_points_ini_model_wing,plot_points_ini_model_LE_per_strut,plot_points_ini_model_canopy_per_strut,plot_points_ini_model_wing,plot_point_ini_model_strut_TE,plot_conn_canopy,plot_conn_tubular_frame)


# #%% Checking if tube_point_indices & plate_point_indices is correct AND if the connectivities are written as functions of the bridle points only
## Checking if tube_point_indices is correct AND if the connectivities are written as functions of the bridle points only

tube_point_indices = []
for point in points_ini_model_bridle_attachment_points:
    indices = np.where(points_ini_model_bridle==point)[0]
    for idx in indices:
        if idx not in tube_point_indices:
            tube_point_indices.append(idx)

points_ini_model_tube = [points_ini_model_bridle[idx] for idx in tube_point_indices]
## plotting results

plot_dummy = [2,[[0,0,0]], 'black',1]
plot_points_ini_model_bridle = [2,points_ini_model_bridle, 'black',2]
plot_points_tubular_frame = [2,points_ini_model_tube, 'purple',6]
plot_points_plate_point_indices = [3,[points_ini_model_bridle[idx] for idx in plate_point_indices], 'orange',4]

#plot_conn_canopy = [points_ini_model_bridle,conn_canopy_i,conn_canopy_j, 'grey',4]
#plot_conn_tubular_frame = [points_ini_model_bridle,conn_tubular_frame_i,conn_tubular_frame_j, 'grey',8]
plot_conn_wing = [points_ini_model_bridle,conn_wing_i,conn_wing_j, 'blue',2]
plot_conn_bridle = [points_ini_model_bridle,ci_bridle_indexed_on_bridle,cj_bridle_indexed_on_bridle, 'black',1]

# get_3D_plot_N_D_6times_PLUS_2lines(plot_points_ini_model_bridle,plot_points_tubular_frame,plot_points_ini_model_pulleys,plot_points_plate_point_indices,plot_dummy,plot_dummy,plot_conn_bridle,plot_conn_wing)


# #%% further processing
## Add a bridle point 
points_bridle_copy = copy.deepcopy(points_ini_model_bridle)
bridlepoint = np.array([-0.55,0,-9.5])
points_bridle = np.vstack((points_bridle_copy,bridlepoint))
bridlepoint_index = len(points_bridle)-1

## Find the closest points to the bridle point
def find_closest_points(point_list, target_index, num_closest=4):
    target_point = point_list[target_index]
    distances = np.linalg.norm(point_list - target_point, axis=1)
    sorted_indices = np.argsort(distances)
    closest_indices = sorted_indices[1:num_closest+1]  # Include the target point
    return closest_indices

# the points closest to the bridlepoint are the kcu points
kcu_point_indices = find_closest_points(points_bridle, bridlepoint_index,num_closest=4)
kcu_plate_indices = [[bridlepoint_index,kcu_point_indices[0],kcu_point_indices[1],bridlepoint_index],[bridlepoint_index,kcu_point_indices[1],kcu_point_indices[2]],[bridlepoint_index,kcu_point_indices[2],kcu_point_indices[3]],[bridlepoint_index,kcu_point_indices[3],kcu_point_indices[0]],kcu_point_indices]

## Adding the KCU lines
kcu_conn_i = [bridlepoint_index,bridlepoint_index,bridlepoint_index,bridlepoint_index,            kcu_point_indices[0],kcu_point_indices[1],kcu_point_indices[2],kcu_point_indices[3],kcu_point_indices[1],kcu_point_indices[2]] 
kcu_conn_j = [kcu_point_indices[0],kcu_point_indices[1],kcu_point_indices[2],kcu_point_indices[3],kcu_point_indices[1],kcu_point_indices[2],kcu_point_indices[3],kcu_point_indices[0],kcu_point_indices[3],kcu_point_indices[0]]
bridle_ci = np.hstack((ci_bridle_indexed_on_bridle,kcu_conn_i))
bridle_cj = np.hstack((cj_bridle_indexed_on_bridle,kcu_conn_j))
kcu_line_indices = np.arange(len(bridle_ci)-len(kcu_conn_i),len(bridle_ci))


# #%% acquiring the tube_line_indices

tube_line_indices = []
for idx, (ci, cj) in enumerate(zip(conn_wing_i, conn_wing_j)):
    # Check if the line (ci, cj) is present in conn_tubular_frame_i and conn_tubular_frame_j
    if (ci, cj) in zip(conn_tubular_frame_i, conn_tubular_frame_j):
        tube_line_indices.append(idx)

# #%% CORRECTING ALL POINTS 

## Correcting all points such that bridle point becomes the origin
points_bridle = points_bridle - bridlepoint
##TODO: this should be done through indices, not through points
## you should only have 1 point list!
for rib in rib_db_whole_model_struts:
    rib[1] = rib[1] - bridlepoint
    rib[2] = rib[2] - bridlepoint
    rib[3] = rib[3] - bridlepoint
    rib[4] = rib[4] - bridlepoint




## PLOTTING TO CHECK
# plot the results
plot_kcu_points = [2,points_bridle[kcu_point_indices], 'red',2]
plot_points_ini_model_pulleys = [2,points_bridle[pulley_point_indices], 'orange',4]
plot_points_bridle = [2,points_bridle, 'black',2]
plot_points_tubular_frame = [2,points_bridle[tube_point_indices], 'purple',3]

plot_conn_kcu = [points_bridle,bridle_ci[kcu_line_indices],bridle_cj[kcu_line_indices],'red',2,[]]
plot_conn_wing = [points_bridle,conn_wing_i,conn_wing_j, 'blue',2,[]]
plot_conn_bridle = [points_bridle,ci_bridle_indexed_on_bridle,cj_bridle_indexed_on_bridle, 'black',1,[]]
# plot_conn_tubular_frame = [points_bridle,conn_tubular_frame_i,conn_tubular_frame_j, 'grey',8,[]]
plot_conn_tubular_frame = [points_bridle,conn_wing_i[tube_line_indices],conn_wing_j[tube_line_indices], 'grey',8,[]]

plot_surface_kcu = [True,points_bridle,kcu_plate_indices,'black',1]
functions_plot.plot_kite([plot_kcu_points,plot_points_bridle,plot_points_tubular_frame,plot_points_ini_model_pulleys],[plot_conn_kcu,plot_conn_wing,plot_conn_bridle,plot_conn_tubular_frame],title='Bridle point added',surface_list=[plot_surface_kcu])


#%% Saving results
## saving to result_bridle_only

# filepath = "initialisation/V9_60C/results_bridle_only/"
filepath = "src/initialisation/V9_60C/results_bridle_only_min_doubles/"
### saving the newly found connectivities - as func of bridle points

# tubular frame
np.save(filepath+'conn_tubular_frame_i.npy',conn_tubular_frame_i )
np.save(filepath+'conn_tubular_frame_j.npy',conn_tubular_frame_j )
np.save(filepath+'tube_point_indices.npy',tube_point_indices )
np.save(filepath+'tube_line_indices.npy',tube_line_indices)

# canopy
np.save(filepath+'conn_canopy_i.npy',conn_canopy_i )
np.save(filepath+'conn_canopy_j.npy',conn_canopy_j )

# wing
np.save(filepath+'conn_wing_i.npy',conn_wing_i )
np.save(filepath+'conn_wing_j.npy',conn_wing_j )
np.save(filepath+'plate_point_indices.npy',plate_point_indices )

### saving the newly found points 
np.save(filepath+'rib_db_whole_model_struts.npy',rib_db_whole_model_struts ) #list with only the struts
np.save(filepath+'points_ini_model_bridle.npy',points_bridle )

## saving old stuff, in the new folder
np.save(filepath+'bridle_ci.npy',bridle_ci )
np.save(filepath+'bridle_cj.npy',bridle_cj ) 

## Bridle line system
np.save(filepath+'pulley_point_indices.npy',pulley_point_indices)
np.save(filepath+'bridlepoint_index.npy',bridlepoint_index)

# KCU THINGS
np.save(filepath+'kcu_point_indices.npy',kcu_point_indices )
np.save(filepath+'kcu_line_indices.npy',kcu_line_indices )
np.save(filepath+'kcu_plate_indices.npy',kcu_plate_indices )

