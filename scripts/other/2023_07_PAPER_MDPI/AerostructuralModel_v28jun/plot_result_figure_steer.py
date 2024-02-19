#%% Printing solution

### DEFINING INPUT
date = '28jun'
delta_ld_used = 8#[%] (depower-tape extension, in percentage)
u_p = 0
billowing_boolean = False

### Define depower-tape extenstion relation to the power-setting u_p
def up_to_ld(u_p,delta_ld_used):
    '''
    input:  u_p         = [0-1]   (power-setting)
            delta_ld    = 8 or 13 [%] (depower-tape extension, in percentage)
    output: depower-tape extension [mm]
    '''
    if delta_ld_used == 8:
        depower_tape_max = 1.482
    elif delta_ld_used == 13:
        depower_tape_max = 1.722
    else:
        raise Exception('delta_ld wrong value, should be: 8 or 13')

    depower_tape_0 = 1.098 #minimum depower tape length

    return 1e3*0.5*(1-u_p)*(depower_tape_max-depower_tape_0)

import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.integrate import solve_ivp
import sys
#sys.path.insert(0, '../functions/') ##this opens a folder one directory higher
sys.path.insert(0, '../AerostructuralModel_v'+date+'/functions/') 
import functions_VSM_LLT as VSM
import functions_PSM as PSM
import functions_plot as functions_plot
import functions_plate_aero as plate_aero
import random

### Loading the red kite file
filename = '../AerostructuralModel_v'+date+'/run_results/pos_up_100_13_False_12.0_v2.csv'

pos = np.loadtxt(filename,delimiter = ',')

## Getting the initial position for LENGTHS
## GEOMETRY AND CONNECTIVITY
## geometry and connectivity
import sys
sys.path.insert(0, '../AerostructuralModel_v'+date+'/geometry_and_connectivity/') 
import functions_connectivity as conn

## Getting the initial position for LENGTHS
#CAD shape
points_CAD = np.loadtxt('../AerostructuralModel_v'+date+'/geometry_and_connectivity/Geometry_modified_kcu.csv',delimiter = ',')

### Initial Guess of POSITION
filename = '../AerostructuralModel_v'+date+'/run_results/pos_up_'+str(int(100*(1)))+'_'+str(delta_ld_used)+'.csv'
pos_initial_guess = np.loadtxt(filename,delimiter = ',')
#pos_initial_guess = points_CAD

## Connectivity grids
ci, cj = conn.get_bridle_line_system_connectivity_KCU() #imports the defined connectivity from the geometry file
ci_kite,cj_kite, plates = conn.get_kite_plate_connectivity()
tube_idx = conn.inflatable_tubes_idx()

## IMPORT BILLOWING using photogrammetry data
sys.path.insert(0, '../AerostructuralModel_v'+date+'/billowing/') 
from Photogrammetry import ballooning_up 

## Moving points to make things SYMMETRICAL
if ballooning_up == True:  #billowing applied to the central points 
    points_CAD[14,1] = points_CAD[14,1]*(ballooning_up(u_p)[0]) # Move points too to avoid asymmetry
    points_CAD[15,1] = points_CAD[15,1]*(ballooning_up(u_p)[0]) # Move points too to avoid asymmetry   

## DYNAMIC SOLVER
## solving the system

## SETTING SPRING LENGTHS
springL = PSM.get_springL_bridle(ci,cj,points_CAD) 
# KITE = spring lengths + billowing applied
ballooning_boolean = False
springL_kite = PSM.get_springL_kite(ci_kite,cj_kite,points_CAD,ballooning_boolean,u_p,ballooning_up)

## ACTUATION INPUT
delta_ld = up_to_ld(u_p,delta_ld_used) # [mm] depower-tape extension
springL[1] += delta_ld #springL[1] corresponds to the depower tape

##TODO: start with slightly shorter springs, for numerical reasons
springL = np.array(springL)*0.999
springL_kite = np.array(springL_kite)*0.999


def plot_kite_pretty(pos,ci,cj,ci_kite,cj_kite,ax,line_boolean,col_kite,elev,plot_params):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
	is is one possible solution to Matplotlib's
	ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

	Input
	  ax: a matplotlib axis, e.g., as output from plt.gca().
      cubes as cubes, etc..  Th
	'''
    width,msize_pulley,msize_lines,col_kite,color_front_lines, color_back_lines, col_surface, alfa_surface,TE_width,tube_width = plot_params


    ### Plot tapes
    tape_colour = col_kite #'#0D23C2'
    tape_width = 3*width

    id1 = 21
    id2 = 27
    ax.plot([pos[id1, 0],pos[id2, 0]],[pos[id1, 1],pos[id2, 1]],[pos[id1, 2],pos[id2, 2]],color = tape_colour,linewidth = tape_width) 
    id2 = 22    
    ax.plot([pos[id1, 0],pos[id2, 0]],[pos[id1, 1],pos[id2, 1]],[pos[id1, 2],pos[id2, 2]],color = tape_colour,linewidth = tape_width) 
    id2 = 23
    ax.plot([pos[id1, 0],pos[id2, 0]],[pos[id1, 1],pos[id2, 1]],[pos[id1, 2],pos[id2, 2]],color = tape_colour,linewidth = tape_width) 
    ##Plot KCU
    ax.plot(pos[21, 0],pos[21, 1],pos[21, 2],color='#000000',linewidth = width,marker = 'v',markersize = msize_pulley*3)


    #plotting the plates    
    for plate in plates:
        p1 = tuple(pos[plate[0], :].tolist())
        p2 = tuple(pos[plate[1], :].tolist())
        p3 = tuple(pos[plate[2], :].tolist())
        p4 = tuple(pos[plate[3], :].tolist())
        verts= [p1,p2,p3,p4]

        # Create the polygon
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
        poly = Poly3DCollection([verts], alpha=alfa_surface, facecolor=col_surface)

        # attempt at shading the surface.. failed
        #from matplotlib import cm
        #from matplotlib.colors import LightSource
        #ls = LightSource(270, 45)
        ## in the rgb colors of the shaded surface calculated from "shade".
        #rgb = np.tile(col_surface, (4,3, 1))
        ##rgb = ls.shade(p1, cmap=cm.gist_earth, vert_exag=0.9, blend_mode='soft')
        #illuminated_surface = ls.shade_rgb(rgb, poly)

        # add the polygon to the plot
        ax.add_collection3d(poly)

        # plot some additional lines on top
        ax.add_collection3d(Line3DCollection([[p1,p2]], colors=col_kite, linewidths=tube_width*2.2, linestyles='-')) #linewidth 3 is for the TE


    #Plot the kite
    for i in range(0, len(ci_kite)):
        
        #LE
        if i in np.arange(0,len(ci_kite),4):  
            ax.plot([pos[ci_kite[i], 0], pos[cj_kite[i], 0]], [pos[ci_kite[i], 1], pos[cj_kite[i], 1]], [pos[ci_kite[i], 2], pos[cj_kite[i], 2]],color=col_kite,linewidth = tube_width)

        #Right-strut
        elif i in np.arange(1,len(ci_kite),4):  
            ax.plot([pos[ci_kite[i], 0], pos[cj_kite[i], 0]], [pos[ci_kite[i], 1], pos[cj_kite[i], 1]], [pos[ci_kite[i], 2], pos[cj_kite[i], 2]],color=col_kite,linewidth = tube_width)
        
        #TE
        elif i in np.arange(2,len(ci_kite),4):  
            ax.plot([pos[ci_kite[i], 0], pos[cj_kite[i], 0]], [pos[ci_kite[i], 1], pos[cj_kite[i], 1]], [pos[ci_kite[i], 2], pos[cj_kite[i], 2]],color=col_kite,linewidth = TE_width)

        #Left-strut
        elif i in np.arange(3,len(ci_kite),4):  
            ax.plot([pos[ci_kite[i], 0], pos[cj_kite[i], 0]], [pos[ci_kite[i], 1], pos[cj_kite[i], 1]], [pos[ci_kite[i], 2], pos[cj_kite[i], 2]],color=col_kite,linewidth = tube_width)

        
 
    # plotting the bridle lines
    for i in range(0, len(ci)):
        if i > 20:
            g = 3
            ax.plot([pos[ci[i], 0], pos[cj[i], 0]], [pos[ci[i], 1], pos[cj[i], 1]], [pos[ci[i], 2], pos[cj[i], 2]],color=color_front_lines,linewidth = width,marker = '.',markersize = msize_lines)  
        else:
            ax.plot([pos[ci[i], 0], pos[cj[i], 0]], [pos[ci[i], 1], pos[cj[i], 1]], [pos[ci[i], 2], pos[cj[i], 2]],color= color_back_lines,linewidth = width,marker = '.',markersize = msize_lines)  
    

    #Plot KCU
    #ax.plot(pos[21, 0],pos[21, 1],pos[21, 2],color='#000000',linewidth = width,marker = 'v',markersize = msize*3)
    #Plot pulleys
    pulley_1 = 23
    pulley_2 = 27
    ax.plot(pos[pulley_1, 0],pos[pulley_1, 1],pos[pulley_1, 2],color='red',linewidth = width,marker = "X",markersize = msize_pulley) 
    ax.plot(pos[pulley_2, 0],pos[pulley_2, 1],pos[pulley_2, 2],color='red',linewidth = width,marker = "X",markersize = msize_pulley)

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()
    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)
    
    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5 * max([x_range, y_range, z_range])
    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])
    plt.tight_layout()

    return

import matplotlib.pyplot as plt

if u_p == 0:
    points_CAD = pos_initial_guess

#PLOTTING
plt.rcParams.update({'font.size': 10})
fig = plt.figure(figsize= (10,10)) # set up the axes for the first plot
ax = fig.add_subplot(1, 1, 1, projection='3d')
elev                = 1 #20
azim                = -90#-130

width               = 0.4 #0.7 #BRIDLES
TE_width            = 1. #0.7 is for full-system
tube_width          = 1. #2.3 is for full-system
msize_pulley        = 0#12
msize_lines         = 0 #9
col_kite            = 'black'
color_front_lines   = 'black' #23B521'
color_back_lines    = 'black'   #0D23C2'
col_surface         = 'grey'  #'whitesmoke'
col_surface         = 'grey'
alfa_surface        = 0.15 # 0.6 is for full-system
#tube_idx_with_diags = [0,1,3,6,7,9,12,13,15,19,21,24,25,27,30,31,33,36,37,39,42,43,45,48,49,51]
plot_params = [width,msize_pulley,msize_lines,col_kite,color_front_lines,color_back_lines,col_surface,alfa_surface,TE_width,tube_width]
plot_kite_pretty(points_CAD/1000,ci,cj,ci_kite,cj_kite,ax,True,'black',elev,plot_params)
#elev                = 0
#azim                = -90
width               = 0 #0.7 #BRIDLES
#TE_width            = 1.5 #0.7 is for full-system
#tube_width          = 1.5 #2.3 is for full-system
msize_pulley        = 0#12
msize_lines         = 0 #9
col_kite            = 'red'
color_front_lines   = 'red' #23B521'
color_back_lines    = 'red'   #0D23C2'
col_surface         = 'red'  #'whitesmoke'
col_surface         = 'red'
alfa_surface        = 0.15 # 0.6 is for full-system
#tube_idx_with_diags = [0,1,3,6,7,9,12,13,15,19,21,24,25,27,30,31,33,36,37,39,42,43,45,48,49,51]
plot_params = [width,msize_pulley,msize_lines,col_kite,color_front_lines,color_back_lines,col_surface,alfa_surface,TE_width,tube_width]
plot_kite_pretty(pos/1000,ci,cj,ci_kite,cj_kite,ax,True,'red',elev,plot_params)
ax.grid(False)
ax.axis('off')
ax.view_init(elev = elev, azim = azim)
bbox = fig.bbox_inches.from_bounds(2, 2, 8, 6)
ax.plot([0,0],[0,0],[0,0],color='black',linewidth = 0.5,marker = '.',markersize = 0.1)
## SAVING FIGURE
filename='PSM_full_particles_side'+str(100*u_p)+'_steer.svg'
print(filename)
plt.savefig(filename,bbox_inches = bbox)

#PLOTTING
plt.rcParams.update({'font.size': 10})
fig = plt.figure(figsize= (10,10)) # set up the axes for the first plot
ax = fig.add_subplot(1, 1, 1, projection='3d')
elev                = 90
azim                = 180
width               = 0 #0.7 #BRIDLES
#TE_width            = 1.5 #0.7 is for full-system
#tube_width          = 1.5 #2.3 is for full-system
msize_pulley        = 0#12
msize_lines         = 0 #9
col_kite            = 'black'
color_front_lines   = 'black' #23B521'
color_back_lines    = 'black'   #0D23C2'
col_surface         = 'grey'  #'whitesmoke'
col_surface         = 'grey'
alfa_surface        = 0.15 # 0.6 is for full-system
#tube_idx_with_diags = [0,1,3,6,7,9,12,13,15,19,21,24,25,27,30,31,33,36,37,39,42,43,45,48,49,51]
plot_params = [width,msize_pulley,msize_lines,col_kite,color_front_lines,color_back_lines,col_surface,alfa_surface,TE_width,tube_width]
plot_kite_pretty(points_CAD/1000,ci,cj,ci_kite,cj_kite,ax,True,'black',elev,plot_params)
width               = 0 #0.7 #BRIDLES
#TE_width            = 1.5 #0.7 is for full-system
#tube_width          = 1.5 #2.3 is for full-system
msize_pulley        = 0#12
msize_lines         = 0 #9
col_kite            = 'red'
color_front_lines   = 'red' #23B521'
color_back_lines    = 'red'   #0D23C2'
col_surface         = 'red'  #'whitesmoke'
col_surface         = 'red'
alfa_surface        = 0.15 # 0.6 is for full-system
#tube_idx_with_diags = [0,1,3,6,7,9,12,13,15,19,21,24,25,27,30,31,33,36,37,39,42,43,45,48,49,51]
plot_params = [width,msize_pulley,msize_lines,col_kite,color_front_lines,color_back_lines,col_surface,alfa_surface,TE_width,tube_width]
plot_kite_pretty(pos/1000,ci,cj,ci_kite,cj_kite,ax,True,'red',elev,plot_params)
ax.grid(False)
ax.axis('off')
ax.view_init(elev = elev, azim = azim)
bbox = fig.bbox_inches.from_bounds(2, 2, 8, 6)
ax.plot([0,0],[0,0],[0,0],color='black',linewidth = 0.5,marker = '.',markersize = 0.1)
## SAVING FIGURE
filename='PSM_full_particles_top'+str(100*u_p)+'_steer.svg'
print(filename)
plt.savefig(filename,bbox_inches = bbox)

#PLOTTING
plt.rcParams.update({'font.size': 10})
fig = plt.figure(figsize= (10,10)) # set up the axes for the first plot
ax = fig.add_subplot(1, 1, 1, projection='3d')
elev                = 15
azim                = 230

width               = 0 #0.7 #BRIDLES
#TE_width            = 1.5 #0.7 is for full-system
#tube_width          = 1.5 #2.3 is for full-system
msize_pulley        = 0#12
msize_lines         = 0 #9
col_kite            = 'black'
color_front_lines   = 'black' #23B521'
color_back_lines    = 'black'   #0D23C2'
col_surface         = 'grey'  #'whitesmoke'
col_surface         = 'grey'
alfa_surface        = 0.15 # 0.6 is for full-system
#tube_idx_with_diags = [0,1,3,6,7,9,12,13,15,19,21,24,25,27,30,31,33,36,37,39,42,43,45,48,49,51]
plot_params = [width,msize_pulley,msize_lines,col_kite,color_front_lines,color_back_lines,col_surface,alfa_surface,TE_width,tube_width]
plot_kite_pretty(points_CAD/1000,ci,cj,ci_kite,cj_kite,ax,True,'black',elev,plot_params)
width               = 0 #0.7 #BRIDLES
#TE_width            = 1.5 #0.7 is for full-system
#tube_width          = 1.5 #2.3 is for full-system
msize_pulley        = 0#12
msize_lines         = 0 #9
col_kite            = 'red'
color_front_lines   = 'red' #23B521'
color_back_lines    = 'red'   #0D23C2'
col_surface         = 'red'  #'whitesmoke'
col_surface         = 'red'
alfa_surface        = 0.15 # 0.6 is for full-system
#tube_idx_with_diags = [0,1,3,6,7,9,12,13,15,19,21,24,25,27,30,31,33,36,37,39,42,43,45,48,49,51]
plot_params = [width,msize_pulley,msize_lines,col_kite,color_front_lines,color_back_lines,col_surface,alfa_surface,TE_width,tube_width]
plot_kite_pretty(pos/1000,ci,cj,ci_kite,cj_kite,ax,True,'red',elev,plot_params)
ax.grid(False)
ax.axis('off')
ax.view_init(elev = elev, azim = azim)
bbox = fig.bbox_inches.from_bounds(2, 2, 8, 6)
ax.plot([0,0],[0,0],[0,0],color='black',linewidth = 0.5,marker = '.',markersize = 0.1)
## SAVING FIGURE
filename='PSM_full_particles_isometric'+str(100*u_p)+'_steer.svg'
print(filename)
plt.savefig(filename,bbox_inches = bbox)


#PLOTTING
plt.rcParams.update({'font.size': 10})
fig = plt.figure(figsize= (10,10)) # set up the axes for the first plot
ax = fig.add_subplot(1, 1, 1, projection='3d')
elev                = 0
azim                = 180

width               = 0 #0.7 #BRIDLES
#TE_width            = 1.5 #0.7 is for full-system
#tube_width          = 1.5 #2.3 is for full-system
msize_pulley        = 0#12
msize_lines         = 0 #9
col_kite            = 'black'
color_front_lines   = 'black' #23B521'
color_back_lines    = 'black'   #0D23C2'
col_surface         = 'grey'  #'whitesmoke'
col_surface         = 'grey'
alfa_surface        = 0.15 # 0.6 is for full-system
#tube_idx_with_diags = [0,1,3,6,7,9,12,13,15,19,21,24,25,27,30,31,33,36,37,39,42,43,45,48,49,51]
plot_params = [width,msize_pulley,msize_lines,col_kite,color_front_lines,color_back_lines,col_surface,alfa_surface,TE_width,tube_width]
plot_kite_pretty(points_CAD/1000,ci,cj,ci_kite,cj_kite,ax,True,'black',elev,plot_params)
width               = 0 #0.7 #BRIDLES
#TE_width            = 1.5 #0.7 is for full-system
#tube_width          = 1.5 #2.3 is for full-system
msize_pulley        = 0#12
msize_lines         = 0 #9
col_kite            = 'red'
color_front_lines   = 'red' #23B521'
color_back_lines    = 'red'   #0D23C2'
col_surface         = 'red'  #'whitesmoke'
col_surface         = 'red'
alfa_surface        = 0.15 # 0.6 is for full-system
#tube_idx_with_diags = [0,1,3,6,7,9,12,13,15,19,21,24,25,27,30,31,33,36,37,39,42,43,45,48,49,51]
plot_params = [width,msize_pulley,msize_lines,col_kite,color_front_lines,color_back_lines,col_surface,alfa_surface,TE_width,tube_width]
plot_kite_pretty(pos/1000,ci,cj,ci_kite,cj_kite,ax,True,'red',elev,plot_params)
ax.grid(False)
ax.axis('off')
ax.view_init(elev = elev, azim = azim)
bbox = fig.bbox_inches.from_bounds(2, 2, 8, 6)
ax.plot([0,0],[0,0],[0,0],color='black',linewidth = 0.5,marker = '.',markersize = 0.1)
## SAVING FIGURE
filename='PSM_full_particles_front'+str(100*u_p)+'_steer.svg'
print(filename)
plt.savefig(filename,bbox_inches = bbox)