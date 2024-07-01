import numpy as np
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.proj3d import proj_transform
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

##TODO: used
def vec_norm(v):
	return np.linalg.norm(v)
#    return np.sqrt(v[0]**2+v[1]**2+v[2]**2)

##TODO: used
def distance(A, B):
	return np.linalg.norm(B - A)
	#return np.sqrt( ((B[0] - A[0]) ** 2) + ((B[1] - A[1]) ** 2) + ((B[2] - A[2]) ** 2) )

##TODO: used
def get_springL_bridle(ci,cj,pos):
	springL = []
	for i in range(0, len(ci)):
		springL.append(distance(pos[ci[i], :], pos[cj[i], :]))
	return springL

##TODO: used
def get_springL_kite(ci_kite,cj_kite,pos,ballooning_boolean,u_p,ballooning_up):
    # kite spring rest lengths
    springL_kite = []
    for i in range(0, len(ci_kite)):
        springL_kite.append(distance(pos[ci_kite[i], :], pos[cj_kite[i], :]))

    if ballooning_boolean == True:
		
        #for i in np.arange(2,len(ci_kite),4): #for the TE occurences
        springL_kite[8]  = springL_kite[8]*ballooning_up(u_p)[3] # 17 -> 18
        # diagonals
        #springL_kite[10] = springL_kite[8]*ballooning_up(u_p)[3]
        #springL_kite[11] = springL_kite[8]*ballooning_up(u_p)[3]

        springL_kite[14] = springL_kite[14]*ballooning_up(u_p)[2] # 16->17
        # diagonals
        #springL_kite[16] = springL_kite[16]*ballooning_up(u_p)[2]
        #springL_kite[17] = springL_kite[17]*ballooning_up(u_p)[2]

        springL_kite[20] = springL_kite[20]*ballooning_up(u_p)[1] # 15- >16
        # diagonals
        #springL_kite[22] = springL_kite[22]*ballooning_up(u_p)[1]
        #springL_kite[23] = springL_kite[23]*ballooning_up(u_p)[1]

        springL_kite[26] = springL_kite[26]*ballooning_up(u_p)[0] # 14- >15           
        # diagonals
        #springL_kite[28] = springL_kite[28]*ballooning_up(u_p)[0]
        #springL_kite[29] = springL_kite[29]*ballooning_up(u_p)[0]

        springL_kite[32] = springL_kite[32]*ballooning_up(u_p)[1] # 13 to 14
        # diagonals
        #springL_kite[34] = springL_kite[34]*ballooning_up(u_p)[1]
        #springL_kite[35] = springL_kite[35]*ballooning_up(u_p)[1]

        springL_kite[38] = springL_kite[38]*ballooning_up(u_p)[2] # 12 to 13
        # diagonals
        #springL_kite[40] = springL_kite[40]*ballooning_up(u_p)[2]
        #springL_kite[41] = springL_kite[41]*ballooning_up(u_p)[2]

        springL_kite[44] = springL_kite[44]*ballooning_up(u_p)[3] # 11 to 12
        #diagonals
        #springL_kite[46] = springL_kite[46]*ballooning_up(u_p)[3]
        #springL_kite[47] = springL_kite[47]*ballooning_up(u_p)[3]

    return springL_kite

# apply box boundary conditions, reverse velocity if outside box
def applyBoundary(pos, vel):
	# kind of the only real boundary is the ground
	is_out = np.where(pos[:, 2] < 0)
	pos[is_out, 2] *= -1
	vel[is_out, 2] *= -1
	return (pos, vel)

def get_factor_area_lst(pos,plates):
	factor_area_lst = []
	for i in np.arange(0,len(plates)): #looping through each panel
		### Alternative factor, based on the area and orientation of the plate wrt x,y plane
		side_1 = distance(pos[plates[i][0]], pos[plates[i][1]])
		side_2 = distance(pos[plates[i][1]], pos[plates[i][2]])
		side_3 = distance(pos[plates[i][2]], pos[plates[i][3]])
		side_4 = distance(pos[plates[i][3]], pos[plates[i][0]])

		# Calculating the semi-perimeter (s)
		# of the given quadilateral
		s = (side_1 + side_2 + side_3 + side_4) / 2

		# Applying Brahmagupta's formula to #https://en.wikipedia.org/wiki/Brahmagupta%27s_formula
		# get maximum area of quadrilateral
		import math
		area = math.sqrt((s - side_1) * (s - side_2) * (s - side_3) * (s - side_4))
		factor_area_lst.append(area) #2919545

	factor_area_lst = np.array(factor_area_lst) / np.max(factor_area_lst) #Non-dimensionalizing it

	return factor_area_lst

##TODO: used
def get_spring_force_kite_no_damp(ci,cj,pos,springL,K_kite,K_diag,K_TE,spring_force,tube_idx):

    for i in range(0,len(ci)): #i loops over each (bridle or kite) line

        sep_vec= pos[ci[i]] - pos[cj[i]] 	# Vector(ci --> cj) separating the points, indicating direction
        sep = vec_norm(sep_vec)		# Absolute magnitude of the vector (works both ways), indicating strength
        unit_vector = sep_vec/sep 			# Define the unit_vector (ci --> cj)
        dL = (sep - springL[i])/1e3 			# SpringL is defined on a range(0,len(ci)) loop (works both ways)
        dL_ratio = sep/springL[i] 			# SpringL is defined on a range(0,len(ci)) loop (works both ways)
        
        K_factor = K_kite
        TE_idx = [2,8,14,20,26,32,38,44,50] ##TODO: weirdly hardcoded?
        #sec = 0

        #IF DIAGONAL         
        if i not in tube_idx and i not in TE_idx: #IF DIAGONAL
            #if i%3 == 0:  
            #    sec +=1            
            K_factor = K_diag
            if dL<0: 
                K_factor = 0
        #IF TE                        
        elif i in TE_idx:
            K_factor = 0
            
            ##TODO: could put this on only give stretch resistance after 5% stretch
            if  dL_ratio > 1 : K_factor = K_TE 
             

        #Apply spring force to ci[i]
        spring_force[ci[i], 0] += -K_factor * dL * unit_vector[0]  # fill x_coord of point in spring_force with: K*dL*unit_vectir
        spring_force[ci[i], 1] += -K_factor * dL * unit_vector[1] 
        spring_force[ci[i], 2] += -K_factor * dL * unit_vector[2] 

        #Apply spring force in opposite direction to cj[i]
        spring_force[cj[i], 0] += -K_factor * dL * -unit_vector[0]   # fill x_coord of point in spring_force with: K*dL*unit_vectir
        spring_force[cj[i], 1] += -K_factor * dL * -unit_vector[1]  
        spring_force[cj[i], 2] += -K_factor * dL * -unit_vector[2] 



    return spring_force


##TODO: used
def plot_kiteforces(pos,ci,cj,ci_kite,cj_kite,lift_force,ax,elev):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
	is is one possible solution to Matplotlib's
	ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

	Input
	  ax: a matplotlib axis, e.g., as output from plt.gca().
      cubes as cubes, etc..  Th
	'''
    setattr(Axes3D, 'arrow3D', _arrow3D)
    width = 0.5
    pos = pos/10
    for i in range(0, len(ci_kite)):
            ax.plot([pos[ci_kite[i], 0], pos[cj_kite[i], 0]], [pos[ci_kite[i], 1], pos[cj_kite[i], 1]], [pos[ci_kite[i], 2], pos[cj_kite[i], 2]],color='black',linewidth = width)#,legend='kite')
    width = 1
    # for i in range(0,len(lift_force)):
    #     ax.plot([pos[i, 0], pos[i, 0]+lift_force[i,0]], [pos[i, 1], pos[i, 1]+lift_force[i,1]], [pos[i, 2], pos[i, 2]+lift_force[i,2]],color='green',linewidth = width)#,legend='kite')
    
    for i in range(0,len(lift_force)):
            a = pos[i]
            b = lift_force[i]
            ax.arrow3D(a[0],a[1],a[2],
                        b[0],b[1],b[2],
                        mutation_scale=5,
                        linewidth = 0.5,
                        arrowstyle="-|>",
                        fc = 'green',
                        ec = 'green'
                        )
            
    
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()
    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)
    ax.view_init(elev = elev, azim = 45)
    
    	# The plot bounding box is a sphere in the sense of the infinity
    	# norm, hence I call half the max range the plot radius.
    plot_radius = 0.5 * max([x_range, y_range, z_range])
    
    # ax.view_init(elev=2000, azim=None)
    # if dist is not None:
    #     ax.dist = dist
    # if elev is not None:
    #     ax.elev = elev
    
    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])
    return

##TODO: used
def plot_kite_deformation(pos,ci,cj,ci_kite,cj_kite,plates,col_kite,ax,line_boolean,elev,tube_idx):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
	is is one possible solution to Matplotlib's
	ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

	Input
	  ax: a matplotlib axis, e.g., as output from plt.gca().
      cubes as cubes, etc..  Th
	'''
    TE = [2,8,14,20,26,32,38,44,50]
    width = 0.7
    msize = 5

    
    # Fill canopy
    for i in range(len(plates)):
        x = np.array([pos[plates[i][0],0],pos[plates[i][1],0],pos[plates[i][2],0],pos[plates[i][3],0]])
        y = np.array([pos[plates[i][0],1],pos[plates[i][1],1],pos[plates[i][2],1],pos[plates[i][3],1]])
        z = np.array([pos[plates[i][0],2],pos[plates[i][1],2],pos[plates[i][2],2],pos[plates[i][3],2]])
       

        # Create a polygonal surface between the lines
        vertices = [(xi, yi, zi) for xi, yi, zi in zip(x, y, z)]
        # vertices += [(xi, yi, zi) for xi, yi, zi in zip(x[::-1], z[::-1], np.zeros_like(x))]
        poly = Poly3DCollection([vertices], alpha=0.1, facecolor=col_kite)
        ax.add_collection3d(poly)
        
    # ax.scatter(pos[:,0], pos[:,1], pos[:,2], color=color)
    
    for i in range(0, len(ci_kite)):
        if i in tube_idx:
            ax.plot([pos[ci_kite[i], 0], pos[cj_kite[i], 0]], [pos[ci_kite[i], 1], pos[cj_kite[i], 1]], [pos[ci_kite[i], 2], pos[cj_kite[i], 2]],color=col_kite,linewidth = width*5,solid_capstyle="round")
        elif i in TE:  
            ax.plot([pos[ci_kite[i], 0], pos[cj_kite[i], 0]], [pos[ci_kite[i], 1],
                                                               pos[cj_kite[i], 1]], [pos[ci_kite[i], 2], pos[cj_kite[i], 2]],color=col_kite,linewidth = width, solid_capstyle="round")
            
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
    plot_radius =  max([0.5*x_range, 0.5*y_range, 0.3*z_range])
    
    # ax.view_init(elev=2000, azim=None)
    # if dist is not None:
    #     ax.dist = dist
    # if elev is not None:
    #     ax.elev = elev
    
    
    
    
    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])
    return

##TODO: used
def plot_kite(pos,ci,cj,ci_kite,cj_kite,plates,slack_index,col_kite,ax,line_boolean,elev,tube_idx):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
	is is one possible solution to Matplotlib's
	ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

	Input
	  ax: a matplotlib axis, e.g., as output from plt.gca().
      cubes as cubes, etc..  Th
	'''
    TE = [2,8,14,20,26,32,38,44,50]
    width = 0.7
    msize = 5
    
    
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
    plt.rcParams.update({'font.size': 11})
    plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"]})
    # ax.scatter(pos[:,0], pos[:,1], pos[:,2], color=color)
    if line_boolean == True:
        for i in range(0, len(ci)):
            color = slack_index[i] 
            ax.plot([pos[ci[i], 0], pos[cj[i], 0]], [pos[ci[i], 1], pos[cj[i], 1]], [pos[ci[i], 2], pos[cj[i], 2]],color=color,linewidth = width)
            
    for i in range(0, len(ci_kite)):
        if i in tube_idx:
            ax.plot([pos[ci_kite[i], 0], pos[cj_kite[i], 0]], [pos[ci_kite[i], 1], pos[cj_kite[i], 1]], [pos[ci_kite[i], 2], pos[cj_kite[i], 2]],color=col_kite,linewidth = width*5,solid_capstyle="round")
        elif i in TE:   
            ax.plot([pos[ci_kite[i], 0], pos[cj_kite[i], 0]], [pos[ci_kite[i], 1], pos[cj_kite[i], 1]], [pos[ci_kite[i], 2], pos[cj_kite[i], 2]],color=col_kite,linewidth = width,solid_capstyle="round")
            
    legend_elements = [Line2D([0], [0], color='#00000080', lw=0.5),
                       Line2D([0], [0], color='#1CEFCC', lw=0.5),
                        Line2D([0], [0], color='#FFBD19', lw=0.5)
                       # Line2D([0], [0], color='#EF1CEC', lw=0.5)
                  ]  
    ax.legend(legend_elements,['no slack','slack $< 2.5\%$','$2.5\% <$ slack $< 5\%$','slack $> 5\%$'], frameon = False,loc = 'center right')
    
    #Plot tapes
    id1 = 21
    id2 = 27
    ax.plot([pos[id1, 0],pos[id2, 0]],[pos[id1, 1],pos[id2, 1]],[pos[id1, 2],pos[id2, 2]],color = '#000000',linewidth = width*3) 
    id2 = 22    
    ax.plot([pos[id1, 0],pos[id2, 0]],[pos[id1, 1],pos[id2, 1]],[pos[id1, 2],pos[id2, 2]],color = '#000000',linewidth = width*3) 
    id2 = 23
    ax.plot([pos[id1, 0],pos[id2, 0]],[pos[id1, 1],pos[id2, 1]],[pos[id1, 2],pos[id2, 2]],color = '#000000',linewidth = width*3) 
    
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
    
    # ax.view_init(elev=2000, azim=None)
    # if dist is not None:
    #     ax.dist = dist
    # if elev is not None:
    #     ax.elev = elev
    
    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])
    return
    

##TODO: used
def plot_kite_pretty(pos,ci,cj,ax,line_boolean,col_kite,elev,tube_idx,ci_kite,cj_kite,plates):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
	is is one possible solution to Matplotlib's
	ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

	Input
	  ax: a matplotlib axis, e.g., as output from plt.gca().
      cubes as cubes, etc..  Th
	'''

    TE = [2,8,14,20,26,32,38,44,50]
    width = 0.7
    msize = 5    
    
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
      
    ip = 1    
    for i in range(0, len(ci_kite)):
        if i in tube_idx:
            ax.plot([pos[ci_kite[i], 0], pos[cj_kite[i], 0]], [pos[ci_kite[i], 1], pos[cj_kite[i], 1]], [pos[ci_kite[i], 2], pos[cj_kite[i], 2]],color=col_kite,linewidth = width*6)
        elif i in TE:  
            ip += 1
            ax.plot([pos[ci_kite[i], 0], pos[cj_kite[i], 0]], [pos[ci_kite[i], 1], pos[cj_kite[i], 1]], [pos[ci_kite[i], 2], pos[cj_kite[i], 2]],color=col_kite,linewidth = width)
            
    for i in range(0, len(ci)):
        if i > 21:
            ax.plot([pos[ci[i], 0], pos[cj[i], 0]], [pos[ci[i], 1], pos[cj[i], 1]], [pos[ci[i], 2], pos[cj[i], 2]],color='#23B521',linewidth = width,marker = '.',markersize = msize)  
        else:
            ax.plot([pos[ci[i], 0], pos[cj[i], 0]], [pos[ci[i], 1], pos[cj[i], 1]], [pos[ci[i], 2], pos[cj[i], 2]],color='#0D23C2',linewidth = width,marker = '.',markersize = msize)  
        
    
   
    


        
    #Plot KCU
    ax.plot(pos[21, 0],pos[21, 1],pos[21, 2],color='#000000',linewidth = width,marker = 'v',markersize = msize*3)
    #Plot pulleys
    ax.plot(pos[24, 0],pos[24, 1],pos[24, 2],color='orange',linewidth = width,marker = 'P',markersize = msize*1.5) 
    ax.plot(pos[28, 0],pos[28, 1],pos[28, 2],color='orange',linewidth = width,marker = 'P',markersize = msize*1.5)
    #Plot tapes
    id1 = 21
    id2 = 27
    ax.plot([pos[id1, 0],pos[id2, 0]],[pos[id1, 1],pos[id2, 1]],[pos[id1, 2],pos[id2, 2]],color = '#0D23C2',linewidth = width*3) 
    id2 = 22    
    ax.plot([pos[id1, 0],pos[id2, 0]],[pos[id1, 1],pos[id2, 1]],[pos[id1, 2],pos[id2, 2]],color = '#0D23C2',linewidth = width*3) 
    id2 = 23
    ax.plot([pos[id1, 0],pos[id2, 0]],[pos[id1, 1],pos[id2, 1]],[pos[id1, 2],pos[id2, 2]],color = '#0D23C2',linewidth = width*3) 
     #Plot KCU
    ax.plot(pos[21, 0],pos[21, 1],pos[21, 2],color='#000000',linewidth = width,marker = 'v',markersize = msize*3)
    
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
    from matplotlib.lines import Line2D
    legend_elements = [Line2D([0], [0], color=col_kite, lw=width*5),
                       Line2D([0], [0], color=col_kite, lw=width),
                        Line2D([0], [0], color='#23B521', lw=width),
                        Line2D([0], [0], color='#0D23C2', lw=width),
                        Line2D([0], [0], color='#0D23C2', lw=width*3),
                        
                       Line2D([0], [0], marker='v', color='w',
                               markerfacecolor='black',markersize=15
                               ),
                       Line2D([0], [0], marker='.', color='w',
                              markerfacecolor='#A605CD',markersize=10),
                       Line2D([0], [0], marker='P', color='w',
                               markerfacecolor='orange',markersize=10
                               )
                   ]    
        
    ax.legend(legend_elements,['Inflatable tubes','Trailing edge','Power lines','Steering lines','Steering tapes','KCU','Point mass particles','Pulleys'], frameon = False,loc = 'center right')
            
    plt.tight_layout()
    return

    
##TODO: used
def plot_kite_pretty_steering(pos,ci,cj,ax,line_boolean,col_kite,elev,tube_idx,ci_kite,cj_kite,plates):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
	is is one possible solution to Matplotlib's
	ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

	Input
	  ax: a matplotlib axis, e.g., as output from plt.gca().
      cubes as cubes, etc..  Th
	'''

    TE = [2,8,14,20,26,32,38,44,50]
    width = 0.7
    msize = 5    
    
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
      
    ip = 1    
    for i in range(0, len(ci_kite)):
        if i in tube_idx:
            ax.plot([pos[ci_kite[i], 0], pos[cj_kite[i], 0]], [pos[ci_kite[i], 1], pos[cj_kite[i], 1]], [pos[ci_kite[i], 2], pos[cj_kite[i], 2]],color=col_kite,linewidth = width*6)
        elif i in TE:  
            ip += 1
            ax.plot([pos[ci_kite[i], 0], pos[cj_kite[i], 0]], [pos[ci_kite[i], 1], pos[cj_kite[i], 1]], [pos[ci_kite[i], 2], pos[cj_kite[i], 2]],color=col_kite,linewidth = width)
            
    for i in range(0, len(ci)):
        if i > 21:
            ax.plot([pos[ci[i], 0], pos[cj[i], 0]], [pos[ci[i], 1], pos[cj[i], 1]], [pos[ci[i], 2], pos[cj[i], 2]],color='#23B521',linewidth = width,marker = '.',markersize = msize)  
        else:
            ax.plot([pos[ci[i], 0], pos[cj[i], 0]], [pos[ci[i], 1], pos[cj[i], 1]], [pos[ci[i], 2], pos[cj[i], 2]],color='#0D23C2',linewidth = width,marker = '.',markersize = msize)  
        
    
   
    


        
    #Plot KCU
    ax.plot(pos[21, 0],pos[21, 1],pos[21, 2],color='#000000',linewidth = width,marker = 'v',markersize = msize*3)
    #Plot pulleys
    ax.plot(pos[24, 0],pos[24, 1],pos[24, 2],color='orange',linewidth = width,marker = 'P',markersize = msize*1.5) 
    ax.plot(pos[28, 0],pos[28, 1],pos[28, 2],color='orange',linewidth = width,marker = 'P',markersize = msize*1.5)
    #Plot tapes
    id1 = 21
    id2 = 27
    ax.plot([pos[id1, 0],pos[id2, 0]],[pos[id1, 1],pos[id2, 1]],[pos[id1, 2],pos[id2, 2]],color = '#0D23C2',linewidth = width*3) 
    id2 = 22    
    ax.plot([pos[id1, 0],pos[id2, 0]],[pos[id1, 1],pos[id2, 1]],[pos[id1, 2],pos[id2, 2]],color = '#0D23C2',linewidth = width*3) 
    id2 = 23
    ax.plot([pos[id1, 0],pos[id2, 0]],[pos[id1, 1],pos[id2, 1]],[pos[id1, 2],pos[id2, 2]],color = '#0D23C2',linewidth = width*3) 
     #Plot KCU
    ax.plot(pos[21, 0],pos[21, 1],pos[21, 2],color='#000000',linewidth = width,marker = 'v',markersize = msize*3)
    #Plot Steerin anchorg
    ax.plot(pos[37, 0],pos[37, 1],pos[37, 2],color='red',linewidth = 3*width,marker = '*',markersize = msize*3)

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
    from matplotlib.lines import Line2D
    legend_elements = [Line2D([0], [0], color=col_kite, lw=width*5),
                       Line2D([0], [0], color=col_kite, lw=width),
                        Line2D([0], [0], color='#23B521', lw=width),
                        Line2D([0], [0], color='#0D23C2', lw=width),
                        Line2D([0], [0], color='#0D23C2', lw=width*3),
                        
                       Line2D([0], [0], marker='v', color='w',
                               markerfacecolor='black',markersize=15
                               ),
                       Line2D([0], [0], marker='.', color='w',
                              markerfacecolor='#A605CD',markersize=10),
                       Line2D([0], [0], marker='P', color='w',
                               markerfacecolor='orange',markersize=10
                               )
                   ]    
        
    ax.legend(legend_elements,['Inflatable tubes','Trailing edge','Power lines','Steering lines','Steering tapes','KCU','Point mass particles','Pulleys'], frameon = False,loc = 'center right')
            
    plt.tight_layout()
    return

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

##TODO: used
def get_spring_force_no_damp(ci,cj,pos,springL,K,spring_force,tube_idx):
	for i in range(0,len(ci)): #i loops over each (bridle or kite) line

		if i == 4 or i == 5: #Hardcode the pulley occurences
			sep_vec_1 = pos[ci[i]] - pos[cj[i]]  # Vector(ci --> cj) separating the points, indicating direction
			sep_1 = vec_norm(sep_vec_1)  # Absolute magnitude of the vector (works both ways), indicating strength
			unit_vector_1 = sep_vec_1 / sep_1  # Define the unit_vector (ci --> cj)

			i_n = i*2 -2 #4 --> 6, 5 --> 8
			sep_vec_2 = pos[ci[i_n]] - pos[cj[i_n]]  # Vector(ci --> cj) separating the points, indicating direction
			sep_2 = vec_norm(sep_vec_2)  # Absolute magnitude of the vector (works both ways), indicating strength
			unit_vector_2 = sep_vec_2 / sep_2  # Define the unit_vector (ci --> cj)



			dL =((sep_1+sep_2) - (springL[i]+springL[i_n]))/1e3 # SpringL is defined on a range(0,len(ci)) loop (works both ways)
            

			### Making a smooth decrease in spring_force, for stability
			if dL >= 0 : K_factor = K*(springL[i]+springL[i_n])*dL/1e3
			elif dL <0	: K_factor = 0  # Think mm's
                   
			# Apply spring force to ci[i]
			spring_force[ci[i], 0] += -K_factor * unit_vector_1[0]  # fill x_coord of point in spring_force with: K*dL*unit_vectir
			spring_force[ci[i], 1] += -K_factor * unit_vector_1[1]  
			spring_force[ci[i], 2] += -K_factor * unit_vector_1[2]  

			# Apply spring force in opposite direction to cj[i]
			spring_force[cj[i], 0] += -K_factor * -unit_vector_1[0]  # fill x_coord of point in spring_force with: K*dL*unit_vectir
			spring_force[cj[i], 1] += -K_factor * -unit_vector_1[1] 
			spring_force[cj[i], 2] += -K_factor * -unit_vector_1[2] 

			# Apply spring force to ci[i_n]
			spring_force[ci[i_n], 0] += -K_factor * unit_vector_2[0]   # fill x_coord of point in spring_force with: K*dL*unit_vectir
			spring_force[ci[i_n], 1] += -K_factor * unit_vector_2[1]  
			spring_force[ci[i_n], 2] += -K_factor * unit_vector_2[2] 

			# Apply spring force in opposite direction to cj[i_n]
			spring_force[cj[i_n], 0] += -K_factor * -unit_vector_2[0]  # fill x_coord of point in spring_force with: K*dL*unit_vectir
			spring_force[cj[i_n], 1] += -K_factor * -unit_vector_2[1] 
			spring_force[cj[i_n], 2] += -K_factor * -unit_vector_2[2] 

		elif i != 6 and i != 8: #Knot bridle points
			sep_vec= pos[ci[i]] - pos[cj[i]] 	# Vector(ci --> cj) separating the points, indicating direction
			sep = vec_norm(sep_vec)		# Absolute magnitude of the vector (works both ways), indicating strength
			unit_vector = sep_vec/sep 			# Define the unit_vector (ci --> cj)
			dL = (sep - springL[i])/1e3  			# SpringL is defined on a range(0,len(ci)) loop (works both ways)

			### Making a smooth decrease in spring_force, for stability
			if dL >= 0 : K_factor = K*springL[i]*dL/1e3
			elif dL <0	: K_factor = 0  # Think mm's
            
			#Apply spring force to ci[i]
            # fill x_coord of point in spring_force with: K*dL*unit_vectir
			spring_force[ci[i], 0] += -K_factor * unit_vector[0] 
			spring_force[ci[i], 1] += -K_factor * unit_vector[1]
			spring_force[ci[i], 2] += -K_factor * unit_vector[2] 

			#Apply spring force in opposite direction to cj[i]
			spring_force[cj[i], 0] += -K_factor * -unit_vector[0]  # fill x_coord of point in spring_force with: K*dL*unit_vectir
			spring_force[cj[i], 1] += -K_factor * -unit_vector[1] 
			spring_force[cj[i], 2] += -K_factor * -unit_vector[2] 

	return spring_force

##TODO: used
def get_spring_force_no_damp_steering(ci,cj,pos,springL,K,spring_force,tube_idx,ratio_K):
    for i in range(0,len(ci)): #i loops over each (bridle or kite) line

        if i == 4 or i == 5: #Hardcode the pulley occurences
            sep_vec_1 = pos[ci[i]] - pos[cj[i]]  # Vector(ci --> cj) separating the points, indicating direction
            sep_1 = vec_norm(sep_vec_1)  # Absolute magnitude of the vector (works both ways), indicating strength
            unit_vector_1 = sep_vec_1 / sep_1  # Define the unit_vector (ci --> cj)

            i_n = i*2 -2 #4 --> 6, 5 --> 8
            sep_vec_2 = pos[ci[i_n]] - pos[cj[i_n]]  # Vector(ci --> cj) separating the points, indicating direction
            sep_2 = vec_norm(sep_vec_2)  # Absolute magnitude of the vector (works both ways), indicating strength
            unit_vector_2 = sep_vec_2 / sep_2  # Define the unit_vector (ci --> cj)



            dL =((sep_1+sep_2) - (springL[i]+springL[i_n]))/1e3 # SpringL is defined on a range(0,len(ci)) loop (works both ways)
            

            ### Making a smooth decrease in spring_force, for stability
            if dL >= 0 : K_factor = K*(springL[i]+springL[i_n])*dL/1e3
            elif dL <0	: K_factor = 0  # Think mm's
                    
            # Apply spring force to ci[i]
            spring_force[ci[i], 0] += -K_factor * unit_vector_1[0]  # fill x_coord of point in spring_force with: K*dL*unit_vectir
            spring_force[ci[i], 1] += -K_factor * unit_vector_1[1]  
            spring_force[ci[i], 2] += -K_factor * unit_vector_1[2]  

            # Apply spring force in opposite direction to cj[i]
            spring_force[cj[i], 0] += -K_factor * -unit_vector_1[0]  # fill x_coord of point in spring_force with: K*dL*unit_vectir
            spring_force[cj[i], 1] += -K_factor * -unit_vector_1[1] 
            spring_force[cj[i], 2] += -K_factor * -unit_vector_1[2] 

            # Apply spring force to ci[i_n]
            spring_force[ci[i_n], 0] += -K_factor * unit_vector_2[0]   # fill x_coord of point in spring_force with: K*dL*unit_vectir
            spring_force[ci[i_n], 1] += -K_factor * unit_vector_2[1]  
            spring_force[ci[i_n], 2] += -K_factor * unit_vector_2[2] 

            # Apply spring force in opposite direction to cj[i_n]
            spring_force[cj[i_n], 0] += -K_factor * -unit_vector_2[0]  # fill x_coord of point in spring_force with: K*dL*unit_vectir
            spring_force[cj[i_n], 1] += -K_factor * -unit_vector_2[1] 
            spring_force[cj[i_n], 2] += -K_factor * -unit_vector_2[2] 

        elif i != 6 and i != 8 and ci[i] != 37: #Knot bridle points
            sep_vec= pos[ci[i]] - pos[cj[i]] 	# Vector(ci --> cj) separating the points, indicating direction
            sep = vec_norm(sep_vec)		# Absolute magnitude of the vector (works both ways), indicating strength
            unit_vector = sep_vec/sep 			# Define the unit_vector (ci --> cj)
            dL = (sep - springL[i])/1e3  			# SpringL is defined on a range(0,len(ci)) loop (works both ways)

            ### Making a smooth decrease in spring_force, for stability
            if dL >= 0 : K_factor = K*springL[i]*dL/1e3
            elif dL <0	: K_factor = 0  # Think mm's
            
            #Apply spring force to ci[i]
            # fill x_coord of point in spring_force with: K*dL*unit_vectir
            spring_force[ci[i], 0] += -K_factor * unit_vector[0] 
            spring_force[ci[i], 1] += -K_factor * unit_vector[1]
            spring_force[ci[i], 2] += -K_factor * unit_vector[2] 

            #Apply spring force in opposite direction to cj[i]
            spring_force[cj[i], 0] += -K_factor * -unit_vector[0]  # fill x_coord of point in spring_force with: K*dL*unit_vectir
            spring_force[cj[i], 1] += -K_factor * -unit_vector[1] 
            spring_force[cj[i], 2] += -K_factor * -unit_vector[2] 

        ### IF SPECIAL STEERINGPOINT
        elif ci[i] == 37:

            sep_vec= pos[ci[i]] - pos[cj[i]] 	# Vector(ci --> cj) separating the points, indicating direction
            sep = vec_norm(sep_vec)		# Absolute magnitude of the vector (works both ways), indicating strength
            unit_vector = sep_vec/sep 			# Define the unit_vector (ci --> cj)
            dL = (sep - springL[i])/1e3  			# SpringL is defined on a range(0,len(ci)) loop (works both ways)

            ### Making a smooth decrease in spring_force, for stability            
            if dL >= 0 : K_factor  = ratio_K*K*springL[i]*dL/1e3
            elif dL <0	: K_factor = 0
            
            #Apply spring force to ci[i]
            # fill x_coord of point in spring_force with: K*dL*unit_vectir
            spring_force[ci[i], 0] += -K_factor * unit_vector[0] 
            spring_force[ci[i], 1] += -K_factor * unit_vector[1]
            spring_force[ci[i], 2] += -K_factor * unit_vector[2] 

            #Apply spring force in opposite direction to cj[i]
            spring_force[cj[i], 0] += -K_factor * -unit_vector[0]  # fill x_coord of point in spring_force with: K*dL*unit_vectir
            spring_force[cj[i], 1] += -K_factor * -unit_vector[1] 
            spring_force[cj[i], 2] += -K_factor * -unit_vector[2] 
             
    return spring_force

##TODO: not used
def get_kite_plate_connectivity_without_diag():

	### Plate connectivity
	plate_1 = [19, 2, 18, 1]
	plate_2 = [2, 3, 17, 18]
	plate_3 = [3, 4, 16, 17]
	plate_4 = [4, 5, 15, 16]
	plate_5 = [5, 6, 14, 15]
	plate_6 = [6, 7, 13, 14]
	plate_7 = [7, 8, 12, 13]
	plate_8 = [8, 9, 11, 12]
	plate_9 = [9, 20, 10, 11]
	plates = [plate_1, plate_2, plate_3, plate_4, plate_5, plate_6, plate_7, plate_8, plate_9]

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
		#ci_kite.append(plates[i][0])
		#cj_kite.append(plates[i][2])

		#ci_kite.append(plates[i][1])
		#cj_kite.append(plates[i][3])

	ci_kite = np.reshape(ci_kite,len(ci_kite))
	cj_kite = np.reshape(cj_kite, len(cj_kite))

	return ci_kite, cj_kite, plates

##TODO: not used, but contains colouring possibilities
def get_solution_print_no_damp_old(pos,points_ini,ci,cj,ci_kite,cj_kite,springL,springL_kite,TE_extension,u_p):
							   
    slack_index = []
    for i in range(0,len(ci)):
        if i > 21:
            # LE bridles
            slack_index.append('#00000080')  
        else:      
            #TE bridles
            slack_index.append('#00000080')  
    
    ### KITE
    dL_kite_stetch_lst,dL_kite_slack_lst = np.zeros(len(ci_kite)),np.zeros(len(ci_kite)) #Initializing, with additional entry for when it remains 0
    for i in range(0, len(ci_kite)):  # i loops over each (bridle or kite) line
        sep_vec = pos[ci_kite[i]] - pos[cj_kite[i]]  # Vector(ci --> cj) separating the points, indicating direction
        sep = vec_norm(sep_vec)  # Absolute magnitude of the vector (works both ways), indicating
        dL_perc = 100*((sep - springL_kite[i]) / sep )
        if dL_perc > 0:  dL_kite_stetch_lst[i] = (dL_perc)  # Append when bridle is longer
        else:       dL_kite_slack_lst[i] = (dL_perc)
    
    ### BRIDLE
    dL_bridle_slack_lst,dL_bridle_stretch_lst = np.zeros(len(ci_kite)),np.zeros(len(ci_kite)) #Initializing, with additional entry for when it remains 0
    for i in range(0, len(ci)):  # i loops over each (bridle or kite) line
        if i == 4 or i == 5: #An attempt at making the pulley line equal colour
            sep_vec_1 = pos[ci[i]] - pos[cj[i]]  # Vector(ci --> cj) separating the points, indicating direction
            sep_1 = vec_norm(sep_vec_1)  # Absolute magnitude of the vector (works both ways), indicating strength
            unit_vector_1 = sep_vec_1 / sep_1  # Define the unit_vector (ci --> cj)

            i_n = i * 2 - 2  # 4 --> 6, 6 --> 8
            sep_vec_2 = pos[ci[i_n]] - pos[cj[i_n]]  # Vector(ci --> cj) separating the points, indicating direction
            sep_2 = vec_norm(sep_vec_2)  # Absolute magnitude of the vector (works both ways), indicating strength
            unit_vector_2 = sep_vec_2 / sep_2  # Define the unit_vector (ci --> cj)

            dL = ((sep_1 + sep_2) - (springL[i] + springL[i_n])) / ( springL[i] + springL[i_n])  # SpringL is defined on a range(0,len(ci)) loop (works both ways)
            dL_perc = 100*dL
            


        elif i!=6 and i!=8:
            sep_vec = pos[ci[i]] - pos[cj[i]]  # Vector(ci --> cj) separating the points, indicating direction
            sep = vec_norm(sep_vec)  # Absolute magnitude of the vector (works both ways), indicating strength
            dL_perc = 100* ((sep - springL[i]) / springL[i])
        elif i == 6:    dL_perc  = dL_bridle_stretch_lst[3]
        elif i == 8:    dL_perc  = dL_bridle_stretch_lst[5]
        
          
        if dL_perc > 0:  
            dL_bridle_stretch_lst[i]=(dL_perc)  # Append when bridle is longer          
        else:       
            dL_bridle_slack_lst[i] = (dL_perc)  
            if dL_perc > -2.5:
                slack_index[i] = '#1CEFCC'
            elif dL_perc <= -2.5 and dL_perc > -5:
                slack_index[i] = '#FFBD19'
            elif dL_perc < -5:
                slack_index[i] = '#EF1CEC'


    print('------------------ Slack & Stretch --------------')
    print('Ave slack kite frame         :', np.round(np.average(dL_kite_slack_lst),2),'%')
    print('Max slack kite  frame        :', np.round(np.min(dL_kite_slack_lst),2),'%')
    print('Ave stretch kite frame       :', np.round(np.average(dL_kite_stetch_lst),2),'%')
    print('Max stretch kite frame       :', np.round(np.max(dL_kite_stetch_lst),2),'%')
    print(' ')
    print('Ave slack bridle system      :', np.round(np.average(dL_bridle_slack_lst),2),'%')
    print('Max slack bridle  system     :', np.round(np.min(dL_bridle_slack_lst),2),'%')
    print('Ave stretch bridle system    :', np.round(np.average(dL_bridle_stretch_lst),2),'%')
    print('Max stretch bridle system    :', np.round(np.max(dL_bridle_stretch_lst),2),'%')
    print(' ')
    return pos,slack_index,dL_bridle_stretch_lst,dL_bridle_slack_lst,dL_kite_stetch_lst,dL_kite_slack_lst

##TODO: used
def get_solution_print_no_damp(pos,points_ini,ci,cj,ci_kite,cj_kite,springL,springL_kite,TE_extension,u_p,tube_idx):
							   
    ### KITE
    #initializing
    dL_tubular_frame_stretch = [0]
    dL_tubular_frame_slack   = [0]
    dL_TE_stretch            = [0]
    dL_TE_slack              = [0]
    dL_diagonal_stretch      = [0]
    dL_diagonal_slack        = [0] 
    #dL_kite_stetch_lst,dL_kite_slack_lst = np.zeros(len(ci_kite)),np.zeros(len(ci_kite)) #Initializing, with additional entry for when it remains 0
    
    
    for i in range(0, len(ci_kite)):  # i loops over each (bridle or kite) line
        sep_vec = pos[ci_kite[i]] - pos[cj_kite[i]]  # Vector(ci --> cj) separating the points, indicating direction
        sep = vec_norm(sep_vec)  # Absolute magnitude of the vector (works both ways), indicating
        #dL_perc = 100*((sep - springL_kite[i]) / sep )
        #dL_perc = sep-springL_kite[i]

        sep_vec_CAD = points_ini[ci_kite[i]] - points_ini[cj_kite[i]]  # Vector(ci --> cj) separating the points, indicating direction
        sep_CAD = vec_norm(sep_vec_CAD)  # Absolute magnitude of the vector (works both ways), indicating
        dL_perc = 100*((sep - sep_CAD) / sep_CAD )
        dL_perc = sep-sep_CAD
        dL_perc = sep/sep_CAD

        TE_idx = [2,8,14,20,26,32,38,44,50]

        # Tubular Frame
        if i in tube_idx: 
            if dL_perc > 0:  # Append when kite element is longer
                dL_tubular_frame_stretch.append(dL_perc) 
            else:       
                dL_tubular_frame_slack.append(dL_perc)
        
        # Trailinge edge
        elif i in TE_idx:
            if dL_perc > 0:  # Append when kite element is longer
                dL_TE_stretch.append(dL_perc)
            else:       
                dL_TE_slack.append(dL_perc)
        
        # Diagonal canopy
        else:
            if dL_perc > 0:  # Append when kite element is longer
                dL_diagonal_stretch.append(dL_perc)  
            else:       
                dL_diagonal_slack.append(dL_perc)
             
    
    ### BRIDLE
    dL_bridle_slack_lst,dL_bridle_stretch_lst = np.zeros(len(ci_kite)),np.zeros(len(ci_kite)) #Initializing, with additional entry for when it remains 0
    for i in range(0, len(ci)):  # i loops over each (bridle or kite) line
        sep_vec = pos[ci[i]] - pos[cj[i]]  # Vector(ci --> cj) separating the points, indicating direction
        sep = vec_norm(sep_vec)  # Absolute magnitude of the vector (works both ways), indicating
        #dL_perc = 100*((sep - springL[i]) / sep )
        #dL_perc = sep-springL_kite[i]

        #sep_vec_CAD = points_ini[ci[i]] - points_ini[cj[i]]  # Vector(ci --> cj) separating the points, indicating direction
        #sep_CAD = vec_norm(sep_vec_CAD)  # Absolute magnitude of the vector (works both ways), indicating
        sep_CAD = springL[i]
        dL_perc = 100*((sep - sep_CAD) / sep_CAD )
        dL_perc = sep-sep_CAD
        dL_perc = sep/sep_CAD


        if dL_perc > 0:  # Append when bridle is longer
            dL_bridle_stretch_lst[i] = (dL_perc)  
        else:       
            dL_bridle_slack_lst[i] = (dL_perc)
        
    print('------------------ Slack & Stretch --------------')
    print('Ave slack bridle system      :', np.round(np.average(dL_bridle_slack_lst),2),'mm')
    print('Max slack bridle  system     :', np.round(np.min(dL_bridle_slack_lst),2),'mm')
    print('Ave stretch bridle system    : ', np.round(np.average(dL_bridle_stretch_lst),2),'mm')
    print('Max stretch bridle system    : ', np.round(np.max(dL_bridle_stretch_lst),2),'mm')
    print(' ')
    print('Ave slack tubular frame      :', np.round(np.average(dL_tubular_frame_slack),2),'mm')
    print('Max slack tubular  frame     :', np.round(np.min(dL_tubular_frame_slack),2),'mm')
    print('Ave stretch tubular frame    : ', np.round(np.average(dL_tubular_frame_stretch),2),'mm')
    print('Max stretch tubular frame    : ', np.round(np.max(dL_tubular_frame_stretch),2),'mm')
    print(' ')
    print('Ave slack TE                 :', np.round(np.average(dL_TE_slack),2),'mm')
    print('Max slack TE                 :', np.round(np.min(dL_TE_slack),2),'mm')
    print('Ave stretch TE               : ', np.round(np.average(dL_TE_stretch),2),'mm')
    print('Max stretch TE               : ', np.round(np.max(dL_TE_stretch),2),'mm')
    print(' ')
    print('Ave slack diagonals          :', np.round(np.average(dL_diagonal_slack),2),'mm')
    print('Max slack diagonals          :', np.round(np.min(dL_diagonal_slack),2),'mm')
    print('Ave stretch diagonals        : ', np.round(np.average(dL_diagonal_stretch),2),'mm')
    print('Max stretch diagonals        : ', np.round(np.max(dL_diagonal_stretch),2),'mm')
    print(' ')

    
    return dL_bridle_stretch_lst,dL_bridle_slack_lst,dL_tubular_frame_stretch,dL_tubular_frame_slack,dL_TE_stretch,dL_TE_slack,dL_diagonal_stretch,dL_diagonal_slack

##TODO: used
def get_aero_forces2nodes(pos,ci,cj,plates,F, M,ringvec,cp,lift_force):
    
    N_struct = len(plates)
    N_split = int(len(cp)/N_struct)
    for i in np.arange(0,len(cp)): #looping through each panel
        sec = (N_struct-1)-int((i+1)/N_split-0.01)
        Fi = (F[i][0] +F[i][1])*vec_norm(ringvec[i]['r0'])
        Mi = M[i]*vec_norm(ringvec[i]['r0'])
        Mi = Mi*cp[i]['airf_coord'][:,2]  
        
        if sec>4:
            Pnodes = np.array([pos[plates[sec][0],:],
                           pos[plates[sec][1],:],
                           pos[plates[sec][2],:],
                           pos[plates[sec][3],:]])/1000
        else:
            Pnodes = np.array([pos[plates[sec][1],:],
                           pos[plates[sec][0],:],
                           pos[plates[sec][3],:],
                           pos[plates[sec][2],:]])/1000
        Fnode = force2nodes(Fi, cp[i]['coordinates_aoa'], Pnodes, cp[i]['tangential'])
        # print(sum(Fnode)-Fi)
        # M1 = cross_product(Pnodes[0]-cp[i]['coordinates_aoa'], Fnode[0,:])
        # M2 = cross_product(Pnodes[1]-cp[i]['coordinates_aoa'], Fnode[1,:])
        # M3 = cross_product(Pnodes[2]-cp[i]['coordinates_aoa'], Fnode[2,:])
        # M4 = cross_product(Pnodes[3]-cp[i]['coordinates_aoa'], Fnode[3,:])
        # MT = M1+M2+M3+M4
        # print(MT)
        Fnode += moment2nodes(Mi, cp[i]['coordinates_aoa'], Pnodes, cp[i]['tangential'])
        
        
        
        if sec>4:
            
            lift_force[plates[sec][0],:] += Fnode[0]
            lift_force[plates[sec][1],:] += Fnode[1]
            lift_force[plates[sec][2],:] += Fnode[2]
            lift_force[plates[sec][3],:] += Fnode[3]
        else: 
            lift_force[plates[sec][1],:] += Fnode[0]
            lift_force[plates[sec][0],:] += Fnode[1]
            lift_force[plates[sec][3],:] += Fnode[2]
            lift_force[plates[sec][2],:] += Fnode[3]
            
    return lift_force
                                                                            			
##TODO: not used?            
def force2nodes(F,Fpoint,nodes,tangential): 
    
    P1 = line_intersect(nodes[0,:],nodes[1,:],Fpoint, Fpoint+tangential)
    d1 = Fpoint-P1
    
    M1 = cross_product(d1, F)
    
    P2 = line_intersect(nodes[3,:],nodes[2,:],Fpoint, Fpoint+tangential)
    
    d2 = P2-P1
    Fp2 = cross_product(M1,d2)/vec_norm(d2)**2
    
    Fp1 = F-Fp2
    
    M3 = cross_product(P1-nodes[0,:], Fp1)
    d3 = nodes[1,:]-nodes[0,:]
    F3 = cross_product(M3,d3)/vec_norm(d3)**2
    
    node1 = Fp1-F3
    node2 = F3
    
    M4 = cross_product(P2-nodes[2,:], Fp2)
    d4 = nodes[3,:]-nodes[2,:]
    F4 = cross_product(M4,d4)/vec_norm(d4)**2
    
    node4 = F4
    node3 = Fp2-F4
    
    Fnode = np.array([node1,
             node2,
             node3,
             node4])
    
    return Fnode


def cross_product(r1,r2):
    
    return np.array([r1[1]*r2[2]-r1[2]*r2[1], r1[2]*r2[0]-r1[0]*r2[2], r1[0]*r2[1]-r1[1]*r2[0]])

def vec_norm(v):
    
    return np.sqrt(v[0]**2+v[1]**2+v[2]**2)

def dot_product(r1,r2):
    return r1[0]*r2[0]+r1[1]*r2[1]+r1[2]*r2[2]

def line_intersect(p1,p2,p3,p4):
    
    p13 = np.empty(3)
    p43 = np.empty(3)
    p21 = np.empty(3)
    pa = np.empty(3)
    pb = np.empty(3)
    
    p13[0] = p1[0] - p3[0];
    p13[1] = p1[1] - p3[1];
    p13[2] = p1[2] - p3[2];
    
    p43[0] = p4[0] - p3[0];
    p43[1] = p4[1] - p3[1];
    p43[2] = p4[2] - p3[2];
    
    
    p21[0] = p2[0] - p1[0];
    p21[1] = p2[1] - p1[1];
    p21[2] = p2[2] - p1[2];
    
    
    d1343 = p13[0] * p43[0] + p13[1] * p43[1] + p13[2] * p43[2];
    d4321 = p43[0] * p21[0] + p43[1] * p21[1] + p43[2] * p21[2];
    d1321 = p13[0] * p21[0] + p13[1] * p21[1] + p13[2] * p21[2];
    d4343 = p43[0] * p43[0] + p43[1] * p43[1] + p43[2] * p43[2];
    d2121 = p21[0] * p21[0] + p21[1] * p21[1] + p21[2] * p21[2];
    
    denom = d2121 * d4343 - d4321 * d4321;
    
    
    
    numer = d1343 * d4321 - d1321 * d4343;
    
    mua = numer / denom;
    mub = (d1343 + d4321 * mua) / d4343;
    
    pa[0] = p1[0] + mua * p21[0];
    pa[1] = p1[1] + mua * p21[1];
    pa[2] = p1[2] + mua * p21[2];
    pb[0] = p3[0] + mub * p43[0];
    pb[1] = p3[1] + mub * p43[1];
    pb[2] = p3[2] + mub * p43[2];

    return pa

def moment2nodes(M,Mpoint,nodes,tangential):  
    d = tangential*0.05
    
    dF = cross_product(M,d)
    dF = dF/vec_norm(dF)
    Fmag = vec_norm(M)/vec_norm(cross_product(dF,d))
    F = dF*Fmag
    
    P1 = Mpoint+d
    
    Fnode1 = force2nodes(F, P1, nodes,tangential)
    Fnode2 = force2nodes(-F, Mpoint, nodes,tangential)

    Fnode = np.array(Fnode1)+np.array(Fnode2)
    
    return Fnode

def get_symmetrical(pos):
        
	pos[22] = [pos[22][0],0,pos[22][2]]
	left_kite  = [5,4,3,2,19,15,16,17,18,1 ]
	right_kite = [6,7,8,9,20,14,13,12,11,10]

	left_bridle_LE  = [24,31,32,33]
	right_bridle_LE = [28,34,35,36]

	left_bridle_TE = [23,24,25,26]
	right_bridle_TE= [27,28,29,30]
    
	def mirroring(pos,left,right):
		for i in range(0,len(left)):
			pos[right[i]] = [pos[left[i]][0],-pos[left[i]][1],pos[left[i]][2]]
		return pos
	
	pos = mirroring(pos,left_kite,right_kite)
	pos = mirroring(pos,left_bridle_LE,right_bridle_LE)
	pos = mirroring(pos,left_bridle_TE,right_bridle_TE)
	
	return pos

def get_F_gravity(pos):
    m_kite = 22.4
    m_KCU = 0.38# 8.4
    m_wing_bridles = m_kite - m_KCU
    m_particle = m_wing_bridles/(len(pos)-1)
    
    F_gravity = np.zeros((len(pos),3))
    for i in range(0,len(pos)):
        #if i = KCU
        if i == 21:
            F_gravity[i] = [0,0,-m_KCU*9.81]
        # else i is some particle
        else:
            F_gravity[i] = [0,0,-m_particle*9.81]

    return F_gravity