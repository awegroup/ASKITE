#%% Defining mega function
## Defining a mega function that runs the analysis for all the cases

def run_full_analysis(delta_ld_used,billowing_boolean):
    print("Begin with: delta_ld_used = "+str(delta_ld_used)+" and billowing_boolean = "+str(billowing_boolean))

    date = '28jun'

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

    def fun_dynamic(t,w):
        """
        Defines the differential equations for the coupled spring-mass system.

        Arguments:
            w :  vector of the state variables:
                    w = [vel1_x,pos1_x,vel1_y,pos1_y,...,veln_z,posn_z]
            t :  time
            p :  vector of the parameters:
                    p = [K, c, L, ci, cj, ci_kite, cj_kite, springL, springL_kite, plates]
        """
        ### Changing up the input vector back to a usable structure
        q_lst,u_lst,i = [],[], 0 #Initializing
        while i < 222: #converting supervector back to coordinate vector
            q_lst.append([w[i  ],w[i+2],w[i+4]]) #positions
            u_lst.append([w[i+1],w[i+3],w[i+5]]) #velocities
            i += 6 #Incrementing it per 3

        pos,vel = np.array(q_lst),np.array(u_lst) #Renaming for consistency

        # making it perfectly symmetrical
        pos = PSM.get_symmetrical(pos) ##TODO: connectivity is hardcoded

        ## forcing a convergence at stoptime
        c_vel = 0.1+0.9*(t/stoptime)
        c_acc = 0.1 *( 1 -t/stoptime)

        # ### Getting the lift_force 
        # lift_force = np.zeros(pos.shape) #Defining a matrix, filled with zeros of similar shape as the pos matrix
        # lift_force = PSM.get_aero_forces2nodes(pos,ci,cj,plates,F_rel, Fmag[:,2],ringvec,controlpoints,lift_force)  # Put in the new positions of the points
        
        ### Getting the lift_force 
        lift_force = np.zeros(pos.shape) #Defining a matrix, filled with zeros of similar shape as the pos matrix
        lift_force = plate_aero.get_lift_force_orientation(plates,pos,lift_force,alpha_0,Vw,equal_boolean)  # Put in the new positions of the points

        ### Getting the spring_force
        spring_force = np.zeros(pos.shape) #Initialising with zero matrix in same shape as points
        spring_force = PSM.get_spring_force_no_damp(ci,cj,pos,springL,K,spring_force,tube_idx) # Spring works only one way
        spring_force = PSM.get_spring_force_kite_no_damp(ci_kite, cj_kite, pos, springL_kite, K_tube,K_diag,K_TE, spring_force,tube_idx) # Spring now works both ways

        ### Filling the vector f = Function that fills equation representation of vector w = [vel1_x,pos1_x,vel1_y,pos1_y,...,veln_z,posn_z]
        f = [] # Initializing
        for i in np.arange(0,len(pos)): #i loops over each point, including the bridle line point
            ### Time to append the values to f
            f.append((1-c_vel)*vel[i][0])                                         # vel in x-direction
            f.append(-c_acc*vel[i][0]    + spring_force[i][0] + lift_force[i][0])   # acc in x-direction
            f.append((1-c_vel)*vel[i][1])                                         # vel in y-direction
            f.append( -c_acc * vel[i][1] + spring_force[i][1] + lift_force[i][1])   # acc in y-direction
            f.append((1-c_vel)*vel[i][2])                                         # vel in z-direction
            f.append(-c_acc*vel[i][2]    + spring_force[i][2] + lift_force[i][2] -0.62*g)   # acc in z-direction

        f[:6]       = np.zeros(6) #Keeping the origin at [0,0,0]

        return f

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
    
    ### Running the whole damn thing
    for i in np.arange(0,.1,0.1):
    # for i in np.arange(0,1.1,0.1):
        u_p = np.round(1-i,1)

        ## case settings 
        u_s = 0     # Steering
        Vw = 20     # [m/s] Wind speed
        alpha_0 = 0
        #delta_ld_used  = 8                     # 8 or 13 [%] (depower-tape extension, in percentage)

        equal_boolean = False   # solving unequal lift-distribution
        direct_boolean = False  # solving direct or dynamic
        #billowing_boolean = True   # True means applying billowing

        g = 9.81 # [m/s^2] gravity

        ## solver settings
        solver = 'Radau'
        abserr = 1e-4 #1.0e-4             # Input to solver as atol -> Determines error control of solver
        relerr = 1e-2 #1.0e-2             # Input to solver as etol -> Determines error control of solver
        stoptime = 30                     # Limit of time
        max_step_size = 5e-3 ##TODO: interesting parameter to play with
        #initial_step = 1e-5  # first step size TODO: interesting parameter to play with    

        ## GEOMETRY AND CONNECTIVITY
        ## geometry and connectivity
        sys.path.insert(0, '../AerostructuralModel_v'+date+'/geometry_and_connectivity/') 
        import functions_connectivity as conn

        ## Getting the initial position for LENGTHS
        #CAD shape
        points_CAD = np.loadtxt('../AerostructuralModel_v'+date+'/geometry_and_connectivity/Geometry_modified_kcu.csv',delimiter = ',')

        ### Initial Guess of POSITION
        if u_p > 0.99: #if u_p equal to 1
            filename = '../AerostructuralModel_v'+date+'/run_results/pos_up_'+str(int(100*(u_p)))+'_presimulated.csv'
        else:
            filename = '../AerostructuralModel_v'+date+'/run_results/pos_up_'+str(int(100*(u_p+0.1)))+'_'+str(delta_ld_used)+'_'+str(billowing_boolean)+'.csv'
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
        if billowing_boolean == True:  #billowing applied to the central points 
            points_CAD[14,1] = points_CAD[14,1]*(ballooning_up(u_p)[0]) # Move points too to avoid asymmetry
            points_CAD[15,1] = points_CAD[15,1]*(ballooning_up(u_p)[0]) # Move points too to avoid asymmetry   

        ## DYNAMIC SOLVER
        ## solving the system

        ## SETTING SPRING LENGTHS
        springL = PSM.get_springL_bridle(ci,cj,points_CAD) 
        # KITE = spring lengths + billowing applied
        #ballooning_boolean = False
        springL_kite = PSM.get_springL_kite(ci_kite,cj_kite,points_CAD,billowing_boolean,u_p,ballooning_up)

        ## ACTUATION INPUT
        delta_ld = up_to_ld(u_p,delta_ld_used) # [mm] depower-tape extension
        springL[1] += delta_ld #springL[1] corresponds to the depower tape

        ##TODO: start with slightly shorter springs, for numerical reasons
        springL = np.array(springL)*0.999
        springL_kite = np.array(springL_kite)*0.999

        #the initital guess
        w0 = [] 
        for i in np.arange(0,len(pos_initial_guess)): #Rewrite to match vector w structure = [vel1_x,pos1_x,vel1_y,pos1_y,...,veln_z,posn_z]
            w0.append(pos_initial_guess[i][0])
            w0.append(1)  # appending the zero initial velocity behind
            w0.append(pos_initial_guess[i][1])
            w0.append(0)  # appending the initial y-velocity behind
            w0.append(pos_initial_guess[i][2])
            w0.append(10) #TODO:  # appending the zero initial velocity behind

        print('Solving the system... (u_p = '+str(u_p)+')')

        #DYNAMIC
        c_acc = 0.1 #0.1 #0.5
        c_vel = 0.1
        K_dynamic   = 2e5 #2e2
        K           = K_dynamic
        K_tube      = K_dynamic #/2 
        K_TE        = K_dynamic #/2
        K_diag      = K_dynamic #/2

        # ## lengthening the TE springs by x% (billowing mimick)
        # perc_TE_lengthen = 2 #2 * (1-u_p)
        # TE_idx = [2,8,14,20,26,32,38,44,50]
        # for i in range(0,len(springL_kite)):
        #     if i in TE_idx:
        #         springL_kite[i] = springL_kite[i]*(1+(perc_TE_lengthen/100))
        
        ## change the middle panels TE, as it stretches to much
        #springL_kite[26] = 0.95*springL_kite[26]

        if billowing_boolean == True:  #billowing applied to the central points
            ## allowing the diagonals to move by 5%
            TE_idx = [2,8,14,20,26,32,38,44,50]
            for i in range(0,len(springL_kite)):
                # if i is a diagonal spring
                if i not in tube_idx and i not in TE_idx: 
                    springL_kite[i] = 1.05*springL_kite[i]

        '''
        #### Define the termination event callback function
        tolerance_solver = 1e-20
        y_previous = w0
        convergence_event.terminal = True #terminate the event if it happens
        convergence_event.direction = 0   #trigger event when switching sign
        # Define the termination event callback function
        def convergence_event(t, y):
            if 'previous_solution' not in convergence_event.__dict__:
                convergence_event.y_previous = y
                return 0  # Termination condition not met
            
            #y_diff > 0 before convergence and < 0 after convergence
            y_diff = abs(np.max((y - convergence_event.y_previous))) - tolerance_solver
            print(y_diff)
            convergence_event.y_previous = y

            return y_diff
        '''
                
        # Solve the ODEs with termination event
        time_start = time.time()
        #wsol = solve_ivp(fun_dynamic, [0,stoptime], w0, method='Radau',events=termination_event,vectorized=False,rtol=relerr,atol=abserr,max_step=max_step_size)#,first_step = initial_step) 
        wsol = solve_ivp(fun_dynamic, [0,stoptime], w0, method=solver, vectorized=False,rtol=relerr,atol=abserr,max_step=max_step_size)#,first_step = initial_step) 
        time_end = time.time()
        print('Time to solve: ',np.round(time_end-time_start,2),'s')

        # Reading out the solution
        w_step = wsol.y[:,-1]
        q_lst,u_lst,i = [],[], 0 #Initializing
        while i < 222: #converting supervector back to coordinate vector
            q_lst.append([w_step[i  ],w_step[i+2],w_step[i+4]]) #positions
            u_lst.append([w_step[i+1],w_step[i+3],w_step[i+5]]) #velocities
            i += 6 #Incrementing it per 3

        pos,vel = np.array(q_lst),np.array(u_lst) #Renaming for consistency

        ### STORING RUNS
        filename = '../AerostructuralModel_v'+date+'/run_results/pos_up_'+str(int(100*u_p))+'_'+str(delta_ld_used)+'_'+str(billowing_boolean)+'.csv'
        np.savetxt(filename,pos,delimiter = ',')


    print("Done with: delta_ld_used = "+str(delta_ld_used)+" and billowing_boolean = "+str(billowing_boolean))
    return 


#%% Running the analysis
## running the analysis

run_full_analysis(8,True)
# run_full_analysis(8,False)
# run_full_analysis(13,True)
# run_full_analysis(13,False)


#%% Plotting width of the results

import numpy as np 

date = '28jun'

w_lst,u_p_lst = [],[]
for i in np.arange(0.0,1.1,0.1):
    u_p = np.round(1-i,1)
    u_p_lst.append(u_p)

    filename = '../AerostructuralModel_v'+date+'/run_results/pos_up_'+str(int(100*u_p))+'_8_billowing.csv'
    pos = np.loadtxt(filename,delimiter = ',')
    w_lst.append(pos[1][1]-pos[10][1])

import matplotlib.pyplot as plt
plt.plot(np.array(u_p_lst),np.array(w_lst).T)

#%% Printing solution
## printing out solution
## Loading Result


import sys
sys.path.insert(0, '../AerostructuralModel_v'+date+'/functions/') 
import functions_VSM_LLT as VSM
import functions_PSM as PSM
import functions_plot as functions_plot
import functions_plate_aero as plate_aero

points_CAD = np.loadtxt('../AerostructuralModel_v'+date+'/geometry_and_connectivity/Geometry_modified_kcu.csv',delimiter = ',')

import functions_connectivity as conn
ci, cj = conn.get_bridle_line_system_connectivity_KCU() #imports the defined connectivity from the geometry file
ci_kite,cj_kite, plates = conn.get_kite_plate_connectivity()
tube_idx = conn.inflatable_tubes_idx()
springL = PSM.get_springL_bridle(ci,cj,points_CAD) 
# KITE = spring lengths + billowing applied
#ballooning_boolean = False
sys.path.insert(0, '../AerostructuralModel_v'+date+'/billowing/') 
from Photogrammetry import ballooning_up 
delta_ld = 8
springL_kite = PSM.get_springL_kite(ci_kite,cj_kite,points_CAD,True,u_p,ballooning_up)


u_p = 0
date = '28jun'
filename = '../AerostructuralModel_v'+date+'/run_results/pos_up_'+str(int(100*u_p))+'_8_True.csv'
pos = np.loadtxt(filename,delimiter = ',')

#print(wsol)
#t_steps = []
#for i in np.arange(0,len(wsol.t)-1):
#    t_steps.append(wsol.t[i+1]-wsol.t[i])
#print("max t_step: ", np.max(t_steps))
#print("ave t_step: ",np.average(t_steps))
#lift_force = np.zeros(pos.shape) #Defining a matrix, filled with zeros of similar shape as the pos matrix
#print('lift_force',plate_aero.get_lift_force_orientation(plates,pos,lift_force,alpha_0,Vw,equal_boolean))
#print('Number of steps              :', wsol.t.size )
#print('Depower tape extension       :', delta_ld, 'mm (u_p: ', u_p, ')')
#print('Diff depower tape(end - beg) :', np.round(np.linalg.norm(pos[22]-pos[21]) - springL[1],2),'mm')

dL_bridle_stretch_lst,dL_bridle_slack_lst,dL_tubular_frame_stretch,dL_tubular_frame_slack,dL_TE_stretch,dL_TE_slack,dL_diagonal_stretch,dL_diagonal_slack = PSM.get_solution_print_no_damp(pos,points_CAD,ci,cj,ci_kite,cj_kite,springL,springL_kite,delta_ld,u_p,tube_idx)


print('Width LE  (OLD)        :', np.round(((-points_CAD[20][1]+points_CAD[19][1])),2),'mm')
print('Width LE (NEW)         :', np.round(((-pos[20][1]+pos[19][1])),2),'mm')
print('Width TE (OLD)         :', np.round(((points_CAD[1][1]-points_CAD[10][1])),2),'mm')
print('Width TE (NEW)         :', np.round(((pos[1][1]-pos[10][1])),2),'mm')
print('Diff Width (New - Old) :', np.round(((pos[1][1]-pos[10][1]) - ((points_CAD[1][1]-points_CAD[10][1]))),2),'mm')

#print('LE vs TE height        :', np.round(pos[4][2],1),np.round(pos[16][2],1))
#print(max(dL_bridle_stretch_lst))


#%% Plotting solution, geometry

plt.rcParams.update({'font.size': 10})
elev = 0
azim = 180

fig = plt.figure(figsize= (10,10))
# set up the axes for the first plot
ax = fig.add_subplot(1, 1, 1, projection='3d')
# PSM.plot_kite(pos/1000,ci,cj,ci_kite,cj_kite,plates,ax,True,'black',elev)
PSM.plot_kite_pretty(pos/1000,ci,cj,ax,True,'black',elev,tube_idx,ci_kite,cj_kite,plates)
#PSM.plot_kite_pretty(points_CAD/1000,ci,cj,ax,True,'red',elev,tube_idx,ci_kite,cj_kite,plates)
ax.grid(False)
ax.axis('off')
ax.view_init(elev = elev, azim = azim)
bbox = fig.bbox_inches.from_bounds(2, 2, 8, 6)

plt.rcParams.update({'font.size': 10})
elev = 90
azim = 180

fig = plt.figure(figsize= (10,10))
# set up the axes for the first plot
ax = fig.add_subplot(1, 1, 1, projection='3d')
# PSM.plot_kite(pos/1000,ci,cj,ci_kite,cj_kite,plates,ax,True,'black',elev)
PSM.plot_kite_pretty(pos/1000,ci,cj,ax,True,'black',elev,tube_idx,ci_kite,cj_kite,plates)
#PSM.plot_kite_pretty(points_CAD/1000,ci,cj,ax,True,'red',elev,tube_idx,ci_kite,cj_kite,plates)
ax.grid(False)
ax.axis('off')
ax.view_init(elev = elev, azim = azim)
bbox = fig.bbox_inches.from_bounds(2, 2, 8, 6)


plt.rcParams.update({'font.size': 10})
elev = 0
azim = 90

fig = plt.figure(figsize= (10,10))
# set up the axes for the first plot
ax = fig.add_subplot(1, 1, 1, projection='3d')
# PSM.plot_kite(pos/1000,ci,cj,ci_kite,cj_kite,plates,ax,True,'black',elev)
PSM.plot_kite_pretty(pos/1000,ci,cj,ax,True,'black',elev,tube_idx,ci_kite,cj_kite,plates)
#PSM.plot_kite_pretty(points_CAD/1000,ci,cj,ax,True,'red',elev,tube_idx,ci_kite,cj_kite,plates)
ax.grid(False)
ax.axis('off')
ax.view_init(elev = elev, azim = azim)
bbox = fig.bbox_inches.from_bounds(2, 2, 8, 6)

# plt.savefig('./plots/structural_model.pdf',bbox_inches = bbox)

## BOTH KITES 
plt.rcParams.update({'font.size': 10})
elev = 90
azim = 90

fig = plt.figure(figsize= (10,10))
# set up the axes for the first plot
ax = fig.add_subplot(1, 1, 1, projection='3d')
# PSM.plot_kite(pos/1000,ci,cj,ci_kite,cj_kite,plates,ax,True,'black',elev)
PSM.plot_kite_pretty(pos/1000,ci,cj,ax,True,'black',elev,tube_idx,ci_kite,cj_kite,plates)
#PSM.plot_kite_pretty(points_CAD/1000,ci,cj,ax,True,'red',elev,tube_idx,ci_kite,cj_kite,plates)
PSM.plot_kite_pretty(points_CAD/1000,ci,cj,ax,True,'green',elev,tube_idx,ci_kite,cj_kite,plates)
ax.grid(False)
ax.axis('off')
ax.view_init(elev = elev, azim = azim)
bbox = fig.bbox_inches.from_bounds(2, 2, 8, 6)

plt.rcParams.update({'font.size': 10})
elev = 0
azim = 90

fig = plt.figure(figsize= (10,10))
# set up the axes for the first plot
ax = fig.add_subplot(1, 1, 1, projection='3d')
# PSM.plot_kite(pos/1000,ci,cj,ci_kite,cj_kite,plates,ax,True,'black',elev)
PSM.plot_kite_pretty(pos/1000,ci,cj,ax,True,'black',elev,tube_idx,ci_kite,cj_kite,plates)
#PSM.plot_kite_pretty(points_CAD/1000,ci,cj,ax,True,'red',elev,tube_idx,ci_kite,cj_kite,plates)
PSM.plot_kite_pretty(points_CAD/1000,ci,cj,ax,True,'green',elev,tube_idx,ci_kite,cj_kite,plates)
ax.grid(False)
ax.axis('off')
ax.view_init(elev = elev, azim = azim)
bbox = fig.bbox_inches.from_bounds(2, 2, 8, 6)

plt.rcParams.update({'font.size': 10})
elev = 0
azim = 180

fig = plt.figure(figsize= (10,10))
# set up the axes for the first plot
ax = fig.add_subplot(1, 1, 1, projection='3d')
# PSM.plot_kite(pos/1000,ci,cj,ci_kite,cj_kite,plates,ax,True,'black',elev)
PSM.plot_kite_pretty(pos/1000,ci,cj,ax,True,'black',elev,tube_idx,ci_kite,cj_kite,plates)
#PSM.plot_kite_pretty(points_CAD/1000,ci,cj,ax,True,'red',elev,tube_idx,ci_kite,cj_kite,plates)
PSM.plot_kite_pretty(points_CAD/1000,ci,cj,ax,True,'green',elev,tube_idx,ci_kite,cj_kite,plates)
ax.grid(False)
ax.axis('off')
ax.view_init(elev = elev, azim = azim)
bbox = fig.bbox_inches.from_bounds(2, 2, 8, 6)

# plt.savefig('./plots/structural_model.pdf',bbox_inches = bbox)

# %%
