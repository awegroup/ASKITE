# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 10:18:54 2023

@author: ocayon

Functions to plot the VSM-PSM in various styles
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.proj3d import proj_transform
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib.lines import Line2D
import functions_PSM as PSM
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
#%% Useful functions
def _arrow3D(ax, x, y, z, dx, dy, dz, *args, **kwargs):
    '''Add an 3d arrow to an `Axes3D` instance.'''

    arrow = Arrow3D(x, y, z, dx, dy, dz, *args, **kwargs)
    ax.add_artist(arrow)
    
class Arrow3D(FancyArrowPatch):

    def __init__(self, x, y, z, dx, dy, dz, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._xyz = (x, y, z)
        self._dxdydz = (dx, dy, dz)

    def draw(self, renderer):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        super().draw(renderer)
        
    def do_3d_projection(self, renderer=None):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))

        return np.min(zs) 


#%% Plotting functions


def plot_VSM_PSM(wingpanels,controlpoints,rings,F,coord_L,pos,ax,col_kite,elev):
    # Aerodynamic part
    width = 1
    mksize = 1
    N_struct = 9
    N_split = int(len(wingpanels)/N_struct)
    secp = 10
    for panel in wingpanels:
        coord = np.array([panel['p1'],panel['p2'],panel['p3']])
        ax.plot(coord[:,0], coord[:,1], coord[:,2], '#000000', linewidth = width, linestyle='--')
        coord = np.array([panel['p3'],panel['p4']])
        ax.plot(coord[:,0], coord[:,1], coord[:,2], '#000000', linewidth = width)
    for i in range(len(wingpanels)):
        sec = (N_struct-1)-int((i+1)/N_split-0.01)
        if sec != secp:
            
            coord = np.array([wingpanels[i]['p1'],wingpanels[i]['p4']])
            ax.plot(coord[:,0], coord[:,1], coord[:,2], '#000000', linewidth = width*5)
        coord = np.array([wingpanels[i]['p1'],wingpanels[i]['p2']])
        ax.plot(coord[:,0], coord[:,1], coord[:,2], '#000000', linewidth = width*5)                  
        secp = sec
    coord = np.array([wingpanels[i]['p2'],wingpanels[i]['p3']])
    ax.plot(coord[:,0], coord[:,1], coord[:,2], '#000000', linewidth = width*5)   
    # for cp in controlpoints:       
        # ax.plot(cp['coordinates'][0],cp['coordinates'][1],cp['coordinates'][2],'orange', marker = '.', markersize = mksize)
        # ax.plot(cp['coordinates_aoa'][0],cp['coordinates_aoa'][1],cp['coordinates_aoa'][2],'#00008B', marker = '.', markersize = mksize)
    for ring in rings:
        for filament in ring:
            
            if filament['id'] == 'trailing_inf1' or filament['id'] == 'trailing_inf2':
                coord = np.array([filament['x1'],filament['x1']+filament['dir']*4])
                # ax.plot(coord[:,0], coord[:,1], coord[:,2], '#0D23C2',linewidth = width,linestyle = '--', alpha = 0.3)
            elif filament['id'] == 'bound':
                coord = np.array([filament['x1'],filament['x2']])
                ax.plot(coord[:,0], coord[:,1], coord[:,2],'#1E90FF',linewidth = width, alpha=0.6)
            else:
                coord = np.array([filament['x1'],filament['x2']])
                # ax.plot(coord[:,0], coord[:,1], coord[:,2],'#0D23C2',linewidth = width,linestyle = '--', alpha=0.6)
                    
    setattr(Axes3D, 'arrow3D', _arrow3D)
    for i in range(len(F)):
        a = coord_L[i]
        b = (F[i][0])/np.linalg.norm(F[int(len(F)/2)][0])*3   
        c = (F[i][1])/np.linalg.norm(F[int(len(F)/2)][1])*1 
        ax.arrow3D(a[0],a[1],a[2],
                    b[0],b[1],b[2],
                    mutation_scale=5,
                    linewidth = width,
                    arrowstyle="-|>",
                    fc = '#2E8B57',
                    ec = '#2E8B57'
                    )
        ax.arrow3D(a[0],a[1],a[2],
                    c[0],c[1],c[2],
                    mutation_scale=5,
                    linewidth = width,
                    arrowstyle="-|>",
                    fc = '#800000',
                    ec = '#800000'
                    )
        
    # Structural part
    ci_kite,cj_kite,plates = PSM.get_kite_plate_connectivity()
    ci, cj = PSM.get_bridle_line_system_connectivity_KCU()
    TE = [2,8,14,20,26,32,38,44,50]
    width = 0.7
    msize = 5
    tube_idx = PSM.inflatable_tubes_idx()
    
    
    # Fill canopy
    for i in range(len(plates)):
        x = np.array([pos[plates[i][0],0],pos[plates[i][1],0],pos[plates[i][2],0],pos[plates[i][3],0]])
        y = np.array([pos[plates[i][0],1],pos[plates[i][1],1],pos[plates[i][2],1],pos[plates[i][3],1]])
        z = np.array([pos[plates[i][0],2],pos[plates[i][1],2],pos[plates[i][2],2],pos[plates[i][3],2]])
       

        # Create a polygonal surface between the lines
        vertices = [(xi, yi, zi) for xi, yi, zi in zip(x, y, z)]
        # vertices += [(xi, yi, zi) for xi, yi, zi in zip(x[::-1], z[::-1], np.zeros_like(x))]
        poly = Poly3DCollection([vertices], alpha=0.1, facecolor='black')
        ax.add_collection3d(poly)
    
        
    # ip = 1    
    # for i in range(0, len(ci_kite)):
    #     if i in tube_idx:
    #         ax.plot([pos[ci_kite[i], 0], pos[cj_kite[i], 0]], [pos[ci_kite[i], 1], pos[cj_kite[i], 1]], [pos[ci_kite[i], 2], pos[cj_kite[i], 2]],color=col_kite,linewidth = width*6)
    #     elif i in TE:  
    #         ip += 1
    #         ax.plot([pos[ci_kite[i], 0], pos[cj_kite[i], 0]], [pos[ci_kite[i], 1], pos[cj_kite[i], 1]], [pos[ci_kite[i], 2], pos[cj_kite[i], 2]],color=col_kite,linewidth = width)
            
    for i in range(0, len(ci)):
        if i > 21:
            ax.plot([pos[ci[i], 0], pos[cj[i], 0]], [pos[ci[i], 1], pos[cj[i], 1]], [pos[ci[i], 2], pos[cj[i], 2]],color='#000000',linewidth = width,marker = '.',markersize = msize)  
        else:
            ax.plot([pos[ci[i], 0], pos[cj[i], 0]], [pos[ci[i], 1], pos[cj[i], 1]], [pos[ci[i], 2], pos[cj[i], 2]],color='#000000',linewidth = width,marker = '.',markersize = msize)  
        
    
   
    


        
    #Plot KCU
    # ax.plot(pos[21, 0],pos[21, 1],pos[21, 2],color='#000000',linewidth = width,marker = 'v',markersize = msize*3)
    #Plot pulleys
    # ax.plot(pos[24, 0],pos[24, 1],pos[24, 2],color='orange',linewidth = width,marker = 'P',markersize = msize*1.5) 
    # ax.plot(pos[28, 0],pos[28, 1],pos[28, 2],color='orange',linewidth = width,marker = 'P',markersize = msize*1.5)
    #Plot tapes
    # id1 = 21
    # id2 = 27
    # ax.plot([pos[id1, 0],pos[id2, 0]],[pos[id1, 1],pos[id2, 1]],[pos[id1, 2],pos[id2, 2]],color = '#333333',linewidth = width*3) 
    # id2 = 22    
    # ax.plot([pos[id1, 0],pos[id2, 0]],[pos[id1, 1],pos[id2, 1]],[pos[id1, 2],pos[id2, 2]],color = '#333333',linewidth = width*3) 
    # id2 = 23
    # ax.plot([pos[id1, 0],pos[id2, 0]],[pos[id1, 1],pos[id2, 1]],[pos[id1, 2],pos[id2, 2]],color = '#333333',linewidth = width*3) 
     #Plot KCU
    # ax.plot(pos[21, 0],pos[21, 1],pos[21, 2],color='#000000',linewidth = width,marker = 'v',markersize = msize*3)
    # Plot point mass particles
    for i in range(0, len(ci)):
        ax.plot(pos[ci[i], 0], pos[ci[i], 1], pos[ci[i], 2],color='#A605CD',linewidth = width,marker = '.',markersize = msize)  
        ax.plot(pos[cj[i], 0], pos[cj[i], 1], pos[cj[i], 2],color='#A605CD',linewidth = width,marker = '.',markersize = msize)
    
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
    
    legend_elements = [Line2D([0], [0], color=col_kite, lw=width*5),
                       Line2D([0], [0], color=col_kite, lw=width, ls='--'),
                       Line2D([0], [0], color='#2E8B57', lw=width),
                       Line2D([0], [0], color='#800000', lw=width),
                       Line2D([0], [0], color='#1E90FF', lw=width),
                       Line2D([0], [0], marker='.', color='w',
                              markerfacecolor='#A605CD',markersize=10)
            
                   ]    
        
    ax.legend(legend_elements,['Inflatable tubes','VSM discretization','Lift force', 'Drag force','Lifting line','Point mass particles'], frameon = False,loc = 'center right')
            