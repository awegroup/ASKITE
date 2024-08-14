#%% -*- coding: utf-8 -*-
#Created on Dec 6 2022
#Python Script to calculate the trajectory of a uniform circular flight of a plate
#@author: Jelle Poland

import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.integrate import solve_ivp

#%% DEFINING FUNCTIONS, initial conditions, connectivity matrix, gravity vector 

def get_w_initial_and_m_vector(length_w,length_pos,N_p,l_t,rho_rope,l_b,m_kcu,x_ini_plate,y_ini_plate,z_ini_plate,m_plate):
    ''' Defines w_initial and m_vector 
        output: w_initial =  [posx1,velx1,posy1,vely1,posz1,velz1,...]
                m_vector = [m1,m1,m1, m2,m2,m2,...]'''

    w_initial = np.zeros(length_w) # making an Nx3 matrix filled with zeros (N rows and 3 columns)
    m_vector = np.zeros(length_pos) #initialising mass_vector

    for i in range(0,N_p):     
        if i == 0:      # origin 
            # initial coordinates
            w_initial[i*6+0] = 0 # x-coordinate
            w_initial[i*6+2] = 0 # y-coordinate
            w_initial[i*6+4] = 0 # z-coordinate
            # mass vector
            m_vector[i*3+0] = 0.5*l_t*rho_rope
            m_vector[i*3+1] = 0.5*l_t*rho_rope
            m_vector[i*3+2] = 0.5*l_t*rho_rope

        elif i == 1:    # p1 
            #(kcu) initial coordinates
            w_initial[i*6+0] = 0 # x-coordinate
            w_initial[i*6+2] = 0 # y-coordinate
            w_initial[i*6+4] = l_t # z-coordinate
            # mass vector
            m_vector[i*3+0] = 0.5*l_t*rho_rope + 4*0.5*l_b+m_kcu
            m_vector[i*3+1] = 0.5*l_t*rho_rope + 4*0.5*l_b+m_kcu
            m_vector[i*3+2] = 0.5*l_t*rho_rope + 4*0.5*l_b+m_kcu

        elif i == 2:    # p2 
            # initial coordinates
            w_initial[i*6+0] = x_ini_plate # x-coordinate
            w_initial[i*6+2] = -y_ini_plate # y-coordinate
            w_initial[i*6+4] = z_ini_plate # z-coordinate
            # mass vector
            m_vector[i*3+0] = 0.5*l_b*rho_rope + 0.25*m_plate
            m_vector[i*3+1] = 0.5*l_b*rho_rope + 0.25*m_plate
            m_vector[i*3+2] = 0.5*l_b*rho_rope + 0.25*m_plate

        elif i == 3:    # p3 
            # initial coordinates
            w_initial[i*6+0] = -x_ini_plate # x-coordinate
            w_initial[i*6+2] = -y_ini_plate # y-coordinate
            w_initial[i*6+4] = z_ini_plate # z-coordinate
            # mass vector           
            m_vector[i*3+0] = 0.5*l_b*rho_rope + 0.25*m_plate
            m_vector[i*3+1] = 0.5*l_b*rho_rope + 0.25*m_plate
            m_vector[i*3+2] = 0.5*l_b*rho_rope + 0.25*m_plate

        elif i == 4:    # p4 
            # initial coordinates
            w_initial[i*6+0] = -x_ini_plate # x-coordinate
            w_initial[i*6+2] = y_ini_plate # y-coordinate
            w_initial[i*6+4] = z_ini_plate # z-coordinate
            # mass vector            
            m_vector[i*3+0] = 0.5*l_b*rho_rope + 0.25*m_plate
            m_vector[i*3+1] = 0.5*l_b*rho_rope + 0.25*m_plate
            m_vector[i*3+2] = 0.5*l_b*rho_rope + 0.25*m_plate

        elif i == 5:    # p5 
            # initial coordinates
            w_initial[i*6+0] = x_ini_plate # x-coordinate
            w_initial[i*6+2] = y_ini_plate # y-coordinate
            w_initial[i*6+4] = z_ini_plate # z-coordinate
            # mass vector            
            m_vector[i*3+0] = 0.5*l_b*rho_rope + 0.25*m_plate
            m_vector[i*3+1] = 0.5*l_b*rho_rope + 0.25*m_plate
            m_vector[i*3+2] = 0.5*l_b*rho_rope + 0.25*m_plate

    return w_initial,m_vector

def get_sd_variable_lst(N_sd,l_t,l_b,l_s):
    ''' Defines spring-damper variables
        input:  N_sd = number of spring-damper elements
                l_t = length of tether
                l_b = length of bridle
                l_s = length of spring
        output: sd_variable_lst = list with all spring-damper variables
                sd_variable_lst = [[sd-1],[sd-2],...]
                [sd-1] = [index_particle_a, index_particle_b, rest_length]'''

    sd_variable_lst = np.zeros((N_sd,3)) # initialising lst with all spring-damper variables
    for i in range(0, N_sd): #looping through all particles

        if i == 0 :                     # sd => tether between p0 and p1
            sd_variable_lst[i][0] = 0               #index of particle a within the pos array
            sd_variable_lst[i][1] = 1               #index of particle b within the pos array
            sd_variable_lst[i][2] = l_t             #rest length of the spring

        elif i >= 1 and i <= 4 :        # sd => bridle between p1 and (p2, p3, p4, p5)
            sd_variable_lst[i][0] = 1        # index of particle a within the pos array
            if   i==1:                  # sd => bridle between p1 and p2
                sd_variable_lst[i][1] = 2    # index of particle b within the pos array
                sd_variable_lst[i][2] = l_b  # rest length of the spring
            elif i==2:                  # sd => bridle between p1 and p3
                sd_variable_lst[i][1] = 3    # index of particle b within the pos array
                sd_variable_lst[i][2] = l_b  # rest length of the spring
            elif i==3:                  # sd => bridle between p1 and p4
                sd_variable_lst[i][1] = 4    # index of particle b within the pos array
                sd_variable_lst[i][2] = l_b  # rest length of the spring
            elif i==4:                  # sd => bridle between p1 and p5
                sd_variable_lst[i][1] = 5    # index of particle b within the pos array
                sd_variable_lst[i][2] = l_b  # rest length of the spring

        elif i >= 5 and i <= 8 :        # sd => plate edges, p2-p3,p3-p4,p4-p5,p5-p2 
            if   i==5: 
                sd_variable_lst[i][0] = 2           #index of particle a within the pos array
                sd_variable_lst[i][1] = 3           #index of particle b within the pos array
            elif i==6: 
                sd_variable_lst[i][0] = 3           #index of particle a within the pos array
                sd_variable_lst[i][1] = 4           #index of particle b within the pos array
            elif i==7: 
                sd_variable_lst[i][0] = 4           #index of particle a within the pos array
                sd_variable_lst[i][1] = 5           #index of particle b within the pos array
                
            elif i==8: 
                sd_variable_lst[i][0] = 5           #index of particle a within the pos array
                sd_variable_lst[i][1] = 2           #index of particle b within the pos array
            
            sd_variable_lst[i][2] = l_s             #rest length of the spring

    return sd_variable_lst

def get_Fg_acc(length_pos,N_p,g):
    ''' Defines gravitational force vector
        input : length_pos
        output: Fg_vector = [0,-g,0,0,-g,0,...]'''
    Fg_acc = np.zeros(length_pos) #initialising
    for i in range(0, N_p): #looping through all particles
        Fg_acc[(i*3)+2] = -g
    return Fg_acc

#%% Functions that change each loop

def get_Fsd_vector(length_pos,N_sd,sd_variable_lst,pos,c,k_AE,vel):
    ''' Defines spring-damper force vector
        input : length_pos,N,pos,delta_L,c,k_AE,vel
        output: Fsd_vector
                Fsd_vector = [Fsd1x,Fsd1y,Fsd1z,Fsd2x,Fsd2y,Fsd2z, ...]'''

    Fsd_vector = np.zeros(length_pos) #initialising
    for i in range(0, N_sd): #looping through all spring-dampers
        
        # defining spring-damper variables in terms of particle indices
        idx_a   = int(sd_variable_lst[i][0])     # index of particle a within the pos array
        idx_b   = int(sd_variable_lst[i][1])     # index of particle b within the pos array
        rest_l  = sd_variable_lst[i][2]          # rest length of spring between particle a and b

        # determining the vector between particles
        sep_vec = pos[idx_b] - pos[idx_a]   # substracting p1 from p0
        sep = np.linalg.norm(sep_vec)       # absolute magnitude of the vector (works both ways), indicating strength
        dL = (sep - rest_l)                 # springL is defined on a range(0,len(ci)) loop (works both ways)
        unit_vector = sep_vec / sep         # define the unit_vector

        # viscous damping coefficient (currently not scaled by the length of the rope)
        Fd_coeff = c #*(rest_l)         # scale damping by length (neglect static friction) 

        # spring-force (currently not making a smooth decrease in spring_force, for stability)
        if dL >= 0 :                        # if the spring is stretched
            Fs = dL * (k_AE /rest_l)      # k = AE/rest_length, here k_AE is the input therefore # Smooth decreas    -> (K*dL+1)
        elif dL< 0 :                        # slack means no spring force
            Fs = 0                          # smooth decrease   -> ((-1 / (dL - 1))) #might need to change from m to mm

        # Apply spring force to a
        Fsd_vector[idx_a + 0] +=  Fs * unit_vector[0] - Fd_coeff*vel[idx_a][0]
        Fsd_vector[idx_a + 1] +=  Fs * unit_vector[1] - Fd_coeff*vel[idx_a][1]
        Fsd_vector[idx_a + 2] +=  Fs * unit_vector[2] - Fd_coeff*vel[idx_a][2]

        # Apply spring force in opposite direction to b
        Fsd_vector[idx_b + 0] += -Fs * unit_vector[0]  + Fd_coeff*vel[idx_b][0]   
        Fsd_vector[idx_b + 1] += -Fs * unit_vector[1]  + Fd_coeff*vel[idx_b][1] 
        Fsd_vector[idx_b + 2] += -Fs * unit_vector[2]  + Fd_coeff*vel[idx_b][2] 

    return Fsd_vector

def get_Fa_vector(length_pos,pos,vel,v_wind,l_s,rho_air,N_p):
    ''' Defines aerodynamic force vector
        input : length_pos,pos,vel,v_wind,l_s,rho_air,N_p
        output: Fa_vector
                Fa_vector = [Fa1x,Fa1y,Fa1z,Fa2x,Fa2y,Fa2z, ...]'''


    # mid_panel_vec, vector describing middle line of panel
    sep_vec23   = pos[3] - pos[2]         # separation vector 'from p2 to p3'
    LE_mid      = 0.5*sep_vec23 + pos[2]  # middle of leading-edge
    sep_vec54   = pos[4] - pos[5]         # separation vector 'from p2 to p3'
    TE_mid      = 0.5*sep_vec54 + pos[5]  # middle of leading-edge
    mid_vec     = TE_mid - LE_mid         # middle vector from TE_mid to LE_mid  

    # vel_a, apparent wind speed vector
    vel_k_LE = 0.5*vel[2]+0.5*vel[3]      # average velocity of p2 and p3 on the leading edge
    vel_w = [0,v_wind,0]                  # wind velocity as a vector
    vel_a = vel_w-vel_k_LE                # apparent wind speed as felt by the kite's LE

    # CL, lift coefficient
    AoA = np.arccos( (np.dot(mid_vec,vel_a)) / (np.linalg.norm(mid_vec) * np.linalg.norm(vel_a) ) ) # angle of attack
    CL = 2*np.pi*AoA             # lift coefficient, assuming a flat plate
    
    # L, lift force
    S = l_s*l_s                                       # surface area of the panel
    L_mag = 0.5*rho_air*(np.linalg.norm(vel_a)**2)*S*CL # lift force

    # L_unit_vec, unit vector of lift force
    diag24_vec = pos[4] - pos[2]                # diagonal vector 'from p2 to p4'
    diag25_vec = pos[5] - pos[2]                # diagonal vector 'from p2 to p5'
    L_vec = np.cross(diag25_vec,diag24_vec)     # vector describing the direction of the lift force
    L_unit_vec = L_vec / np.linalg.norm(L_vec)  # unit vector of the lift force

    Fa_vector = np.zeros(length_pos) #initialising
    for i in range(0,N_p):     
        if i > 1: # if the particle is one of the panel corners
            Fa_vector[(i*3)+0] = L_mag * L_unit_vec[0] # apply the lift force in x-direction
            Fa_vector[(i*3)+1] = L_mag * L_unit_vec[1] # apply the lift force in y-direction
            Fa_vector[(i*3)+2] = L_mag * L_unit_vec[2] # apply the lift force in z-direction   

    return Fa_vector

#%% Defining the spring-damper system function
def vectorfield(t,w):
    """
    Defines the differential equations for the coupled spring-mass system.
    Calculates acceleration per particle, which should equal
            a = Fg/m + Fs/m + Fd/m = Fg/m + Fsd/m
        m  : mass of particle
            m[0] = m_block + 0.5*m_rope_section
            m[i] = m_rope_section ( = 0.5*m_rope_section_above+0.5*m_rope_section_below)
        Fg/m : Gravitational acceleration
            Fg/m[i_y] = g
        Fsd/m : acceleration due to spring-damping force
            Fsd/m[i] = [K * (springL_new - springL_old) - C *vel ]/ m

    input:  t = time
            w = [posx1,velx1,posy1,vely1,posz1,velz1,...]

    output: f = [velx1,accx1,vely1,accy1,velz1,accz1,...]
    """

    ### Changing the input vector back to a usable structure
    pos_lst,vel_lst,i = [],[], 0                    #Initializing
    while i < (length_w-5):                         #converting supervector-w back to coordinate vector-pos
        pos_lst.append([w[i],w[i+2],w[i+4]])        #positions
        vel_lst.append([w[i+1],w[i+3],w[i+5]])      #velocities
        i += 6                                      #Incrementing it per 6, because we have 3pos and 3vel per particle
    pos,vel = np.array(pos_lst),np.array(vel_lst)   #Making lists into arrays

    ### Applying the actuation
    sd_variable_lst[3][2] += -u_p + u_s
    sd_variable_lst[4][2] += -u_p - u_s

    ### Calculating the acceleration per particle
    # m_vector = already defined by get_w_initial_and_m_vector
    # Fg_acc   = NOT NEEDED, we do without gravity first                                         
    Fsd_acc = get_Fsd_vector(length_pos,N_sd,sd_variable_lst,pos,c,k_AE,vel)/m_vector  # getting spring-damper forces vector and transform to acceleration 
    Fa_acc  = get_Fa_vector(length_pos,pos_ini,vel_ini,v_wind,l_s,rho_air,N_p)/m_vector                                              # getting aerodynamic forces vector and transform to acceleration
    acc     = np.reshape( (Fsd_acc+Fa_acc), (N_p,3))                                        # summing up the accelerations

    print('vel =',vel)
    print('acc = ',acc)
    ### Filling the output vector f, with the velocity and acceleration
    f = [0,0,0,0,0,0]            # initializing, keeping the origin fixed
    for i in range(1,N_p):       # i loops over each point, except the origin
        f.append(vel[i][0])      # vel in x-direction
        f.append(acc[i][0])      # acc in x-direction
        f.append(vel[i][1])      # vel in y-direction
        f.append(acc[i][1])      # acc in y-direction
        f.append(vel[i][2])      # vel in z-direction
        f.append(acc[i][2])      # acc in z-direction
    return f


#%% INPUT PARAMETERS

# physical constants
g = 9.81            #[m/s^2]    gravitational acceleration
k_AE = 4E3          #[N]        cross-sectional area (A) * modulus of elasticity (E), k = AE/L
c = 4E2             #[N/m*s^-1] viscous damping coefficient (static-friction is neglected)
rho_rope = 0.1      #[kg/m]     LINEAR density of the rope
m_plate  = 10       #[kg]       mass of the plate
m_kcu   = 2         #[kg]       mass of the kcu
l_t = 10            #[m]        length of the tether
l_s = 1             #[m]        length of the plate sides
l_b = 1.5           #[m]        length of the bridles
N_p = 6             #[-]        number of particles
N_sd = 9            #[-]        number of spring-damper elements
v_wind = 10         #[m/s]      wind velocity
rho_air = 1.225     #[kg/m^3]   air density

# kite settings
u_p = 0.2 #[m] depower tape ,u_p = 0 = fully-powered kite
u_s = 0.0 #[m] steering tape,u_s = no steering, u_s =0.1 steering to the right

# solver settings
stoptime = 6           # limit of time
abserr = 1e-4          # absolute tolerance
relerr = 1e-3          # relative tolerance 
#t_eval_lst = [0,0.25*stoptime,0.5*stoptime,stoptime] # deteriming the time steps that you want to see the output of
solver_type = 'BDF'


#%% Initial conditions

# calculated coordinates
x_ini_plate = l_s/2 # initial x-coordinates of the plate particles
y_ini_plate = l_s/2 # initial y-coordinates of the plate particles
plate_diag = np.sqrt(l_s**2 + l_s**2) # diagonal of the plate
z_ini_plate = l_t + np.sqrt(l_b**2-(plate_diag/2)**2) # initial z-coordinates of the plate particles

# calculated input parameters
length_pos  = (N_p)*3     #[-]  length of the vectors: pos,vel,acc
length_w    = 2*length_pos  #[-]  length of w  

w_initial,m_vector = get_w_initial_and_m_vector(length_w,length_pos,N_p,l_t,rho_rope,l_b,m_kcu,x_ini_plate,y_ini_plate,z_ini_plate,m_plate)
sd_variable_lst    = get_sd_variable_lst(N_sd,l_t,l_b,l_s)

pos_lst,vel_lst,i = [],[], 0 #Initializing
while i < (length_w): #converting supervector w_step back to coordinate pos. & vel. vectors
    pos_lst.append([w_initial[i],w_initial[i+2],w_initial[i+4]]) #positions
    vel_lst.append([w_initial[i+1],w_initial[i+3],w_initial[i+5]]) #velocities
    i += 6 #Incrementing it per 3
pos_ini,vel_ini = np.array(pos_lst),np.array(vel_lst) #Renaming for consistency

print('f',vectorfield(0,w_initial)) #checking if the vectorfield works


#%% Solving the ODE
# SOLVER WORKINGS
# "vectorfield" has:            (type : 1st order differential equation)
#             input: pos & vel  (shape: [pos1x,vel1x,pos1y,vel1y,pos1z,vel1z, pos2x,vel2x,pos2y,vel2y,pos2z,vel2z, ...] )
#            output: vel & acc  (shape: [vel1x,acc1x,vel1y,acc1y,vel1z,acc1z, vel2x,acc2x,vel2y,acc2y,vel2z,acc2z, ...] )
# Solver has:                   (type : 4th order Runga-Kutta method, with adaptive step size)
#            input: vel & acc   (shape: [vel1x,acc1x,vel1y,acc1y,vel1z,acc1z, vel2x,acc2x,vel2y,acc2y,vel2z,acc2z, ...] )
#           output: pos & vel   (shape: [pos1x,vel1x,pos1y,vel1y,pos1z,vel1z, pos2x,vel2x,pos2y,vel2y,pos2z,vel2z, ...] )

#Solving
start_time = time.time() #define start time
wsol = solve_ivp(vectorfield, [0,stoptime], w_initial, method=solver_type,atol=abserr,rtol=relerr)#,t_eval = t_eval_lst)
end_time = time.time()  

#%% Solution

# Reading out the solution
w_step = wsol.y[:,-1]

#converting supervector w_step back to coordinate pos. & vel. vectors
pos_lst,vel_lst,i = [],[], 0 #Initializing
while i < (length_w): #converting supervector w_step back to coordinate pos. & vel. vectors
    pos_lst.append([w_step[i],w_step[i+2],w_step[i+4]]) #positions
    vel_lst.append([w_step[i+1],w_step[i+3],w_step[i+5]]) #velocities
    i += 6 #Incrementing it per 3

pos_new,vel_new = np.array(pos_lst),np.array(vel_lst) #Renaming for consistency
diff = np.amax(np.linalg.norm(pos_new-pos_ini,axis = 1)) #Defining the max difference in position of the points

# printing out the solution
print("-------------------------------------------------------------")
print("-------------------------- NEW RUN --------------------------")
print("-------------------------------------------------------------")
print('Total time                   : ', round( (end_time-start_time),6),'[s]')
print('Evaluations required         : ', wsol.nfev, '    [-]')
print('Max difference in point pos. : ',np.round(diff,6),'[m]')
print('pos_ini (y-coordinates)      : ',np.round(pos_ini[:,1],3))
print('pos_new (y-coordinates)      : ',np.round(pos_new[:,1],3))
#print('Evaluation time points       : ')
#print(np.round(wsol.t,3))
#print('Evaluation  points           : ')
#print(np.round(wsol.y,3))

'''
#%% PLOTTING the solution
import matplotlib.pyplot as plt
block_offset_for_print = np.round(block_offset,3)
plt.plot(wsol.t, wsol.y[2]-L*0.1,label='block (offset = '+str(block_offset_for_print)+'m)')
for i in range(1,N-1):
    plt.plot(wsol.t, wsol.y[2+i*6]- (L*0.1+i*delta_L),label='particle '+str(i))
plt.xlabel('time /s')
plt.ylabel('Particle Deviations /m')
plt.grid()
plt.legend()
plt.show()
'''