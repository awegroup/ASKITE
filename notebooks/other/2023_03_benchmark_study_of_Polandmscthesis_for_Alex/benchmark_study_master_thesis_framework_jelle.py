# %% -*- coding: utf-8 -*-
"""
Created on Nov 28 10:35

@author: Jelle Poland
"""

import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.integrate import solve_ivp

#  DEFINING FUNCTIONS


def get_mass(length_pos, m_block, delta_L, rho_rope):
    ''' Defines mass of each particle
            input : length_pos,m_block,delta_L,rho_rope
            output: m_vector
                    m_vector = [m1,m1,m1 ,m2,m2,m2, ...]'''

    m_vector = np.zeros(length_pos)
    for i in range(0, (length_pos)):
        if i < 3:  # Defining the mass of the first particle when i = 0,1,2
            m_vector[i] = m_block + 0.5*delta_L*rho_rope
        else:  # Defining the mass of the other particles
            m_vector[i] = delta_L*rho_rope
    return m_vector


def get_Fg_acc(length_pos):
    ''' Defines gravitational force vector
        input : length_pos
        output: Fg_vector = [0,-g,0,0,-g,0,...]'''
    Fg_acc = np.zeros(length_pos)  # initialising
    for i in range(0, (N)):  # looping through all particles
        Fg_acc[(i*3)+1] = -g
    return Fg_acc


def get_Fsd_vector(length_pos, N, pos, delta_L, c_rek, k_AE, vel):
    ''' Defines spring-damper force vector
            input : length_pos,N,pos,delta_L,c,K,vel
            output: Fsd_vector
                    Fsd_vector = [Fsd1x,Fsd1y,Fsd1z,Fsd2x,Fsd2y,Fsd2z, ...]'''

    Fsd_vector = np.zeros(length_pos)  # initialising
    for i in range(0, (N)):  # looping through all particles
        if i < (N-1):  # looping through all the particles, except the ceiling particle, as its not needed

            # determining the vector between particles
            # substracting the block mass from the particle above
            sep_vec = pos[i+1] - pos[i]
            # absolute magnitude of the vector (works both ways), indicating strength
            l_element = np.linalg.norm(sep_vec)
            # springL is defined on a range(0,len(ci)) loop (works both ways)
            dL = (l_element - delta_L)
            # define the unit_vector n the longitudinal direc.
            unit_long = sep_vec / l_element

            # spring_force ( currently without smooth decrease )
            if dL >= 0:                 # if the spring is stretched
                # k = AE/rest-length # Smooth decrease   -> (K*dL+1)
                Fs = dL * (k_AE/delta_L)
            elif dL < 0:                 # slack means no spring force
                # Smooth decrease   -> ((-1 / (dL - 1))) #might need to change from m to mm
                Fs = 0

            # material damping coefficient (not scaled by the length of the rope)
            # trying to get the damping only in the direction of the relative velocity
            dot_vel_1 = np.dot(vel[i], unit_long)
            dot_vel_2 = np.dot(vel[i+1], unit_long)

            # the relative velocity between the two particles than becomes
            rel_vel = dot_vel_1-dot_vel_2
            Fd = c_rek*rel_vel  # *(delta_L)

            ''' NOT NEEDED IN BENCHMARK
            ### viscous aerodynamic damping, air-friction essentially (MSc Thesis N.Geschiere is excellent reference)
            V_b = 0.5*(vel[i] + vel[i+1])       # velocity of the bridle = average velocity of the two particles
            V_b_app = v_wind - vel_ave          # apparent velocity of bridle
            V_b_norm = np.linalg.norm(V_b_app)  # norm of apparent velocity
            # derivation of equation below, see "Bridle Particle pdf"
            S_eff_bridle = diam_bridle * np.sqrt( l_element**2 - (np.dot(V_b_app,sep_vec)/V_b_norm)**2  ) 
            Fa_drag = 0.5 * rho_air * (V_b_app*V_b_norm) * S_eff_bridle * C_d_bridle  # Drag force, includes the direction of the velocity
            '''

            # apply spring force to i
            # (NOT NEEDED BENCHMARK) #-0.5*Fa_drag*unit_trans[0]
            Fsd_vector[(i * 3) + 0] += Fs * unit_long[0] - Fd*unit_long[0]
            # (NOT NEEDED BENCHMARK) #-0.5*Fa_drag*unit_trans[1]
            Fsd_vector[(i * 3) + 1] += Fs * unit_long[1] - Fd*unit_long[1]
            # (NOT NEEDED BENCHMARK) #-0.5*Fa_drag*unit_trans[2]
            Fsd_vector[(i * 3) + 2] += Fs * unit_long[2] - Fd*unit_long[2]

            # apply spring force in opposite direction to i+1
            # (NOT NEEDED BENCHMARK) #-0.5*Fa_drag*unit_trans[0]
            Fsd_vector[((i+1) * 3) + 0] += -Fs * unit_long[0] + Fd*unit_long[0]
            # (NOT NEEDED BENCHMARK) #-0.5*Fa_drag*unit_trans[1]
            Fsd_vector[((i+1) * 3) + 1] += -Fs * unit_long[1] + Fd*unit_long[1]
            # (NOT NEEDED BENCHMARK) #-0.5*Fa_drag*unit_trans[2]
            Fsd_vector[((i+1) * 3) + 2] += -Fs * unit_long[2] + Fd*unit_long[2]

    return Fsd_vector


def get_w_initial(length_w, h_ceiling, L, block_offset, N):
    ''' Defines w_initial 
        input : length_w,h_ceiling,L,block_offset
        output: w_initial
                w_initial =  [posx1,velx1,posy1,vely1,posz1,velz1,...]'''

    # making an Nx3 matrix filled with zeros (N rows and 3 columns)
    w_initial = np.zeros(length_w)
    h_block_rest = h_ceiling-L  # calculating height of the block in rest
    for i in range(0, N):
        if i == 0:      # determining the initial z-coordinate of the block
            w_initial[i*6+2] = h_block_rest + block_offset/N
        elif i < (N-1):  # dealing with all the rope particles, except the last one
            w_initial[i*6+2] = h_block_rest + delta_L*i + block_offset/N
        else:           # setting the z-coordinate of the ceiling particle
            w_initial[i*6+2] = h_ceiling
    return w_initial

#  Defining the spring-damper system function


def vectorfield(t, w):
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

    # Changing the input vector back to a usable structure
    pos_lst, vel_lst, i = [], [], 0  # Initializing
    while i < (length_w-5):  # converting supervector-w back to coordinate vector-pos
        pos_lst.append([w[i], w[i+2], w[i+4]])  # positions
        vel_lst.append([w[i+1], w[i+3], w[i+5]])  # velocities
        i += 6  # Incrementing it per 6, because we have 3pos and 3vel per particle
    pos, vel = np.array(pos_lst), np.array(vel_lst)  # Making lists into arrays

    # Calculating the acceleration per particle
    m_vector = get_mass(length_pos, m_block, delta_L,
                        rho_rope)             # getting mass vector
    # getting vector filled with g [g,g,g,..]
    Fg_acc = get_Fg_acc(length_pos)
    # getting spring-damper forces vector and transform to acceleration
    # Fsd_acc = get_Fsd_vector(length_pos, N, pos, delta_L,
    #                         c, k_AE, vel, rho_air, C_d_bridle) / m_vector
    Fsd_acc = get_Fsd_vector(length_pos, N, pos, delta_L,
                             c_rek, k_AE, vel) / m_vector
    # summing up the accelerations
    acc = np.reshape((Fg_acc+Fsd_acc), (N, 3))

    # Filling the output vector f, with the velocity and acceleration
    f = []                       # Initializing
    for i in range(0, N):         # i loops over each point, including the bridle line point
        f.append(vel[i][0]*(1-c_snel))      # vel in x-direction
        f.append(acc[i][0])      # acc in x-direction
        f.append(vel[i][1]*(1-c_snel))      # vel in y-direction
        f.append(acc[i][1])      # acc in y-direction
        f.append(vel[i][2]*(1-c_snel))      # vel in z-direction
        f.append(acc[i][2])      # acc in z-direction

    f[-6:] = np.zeros(6)   # keeping the top particle, attached to the ceiling.

    return f


# INPUT PARAMETERS

# physical constants
g = 9.81  # [m/s^2]    gravitational acceleration
k_AE = 6E4  # [N/m]      stiffness of the rope
# 98846.8  # [N/m*s^-1] viscous damping coefficient (excludes static friction)
c_rek = 0
c_snel =0
L = 10  # [m]        length of the rope
rho_rope = 0.1  # [kg/m]     LINEAR density of the rope
m_block = 1E2  # [kg]       Mass of the block
# [m]        Height of the ceiling where the rope is attached to, 10% margin with rope length is assumed
h_ceiling = L*1.1
rho_air = 1.225  # [kg/m^3]   ISA-standard air density
diam_bridle = 0.02  # [m]        diameter of bridle lines
# C_d_bridle = 1.05  # [-]        drag-coefficient of bridles

# simulation parameters
N = 4  # [-]          #particles, [0=block,1:N-1 = rope, N = ceiling particle]
# [m]         #INITIAL block offset, that starts the whole vibration motion
block_offset = 0.1


# solver settings
stoptime = 5         # limit of time
abserr = 1e-4          # absolute tolerance
relerr = 1e-3          # relative tolerance
# t_eval_lst = [0,0.25*stoptime,0.5*stoptime,stoptime] # deteriming the time steps that you want to see the output of
solver_type = 'Radau'

# CALCULATED input parameters
delta_L = L/(N-1)  # [m]  distance between the particles in the rope
length_pos = (N)*3  # [-]  length of the vectors: pos,vel,acc
length_w = 2*length_pos  # [-]  length of w

# % Initial Conditions

# Initial conditions
w_initial = get_w_initial(length_w, h_ceiling, L, block_offset, N)
pos_lst, vel_lst, i = [], [], 0  # Initializing
while i < (length_w):  # converting supervector w_step back to coordinate pos. & vel. vectors
    pos_lst.append([w_initial[i], w_initial[i+2], w_initial[i+4]])  # positions
    vel_lst.append([w_initial[i+1], w_initial[i+3],
                   w_initial[i+5]])  # velocities
    i += 6  # Incrementing it per 3
pos_ini, vel_ini = np.array(pos_lst), np.array(
    vel_lst)  # Renaming for consistency


#  Solving
# SOLVER WORKINGS
# "vectorfield" has:            (type : 1st order differential equation)
#             input: pos & vel  (shape: [pos1x,vel1x,pos1y,vel1y,pos1z,vel1z, pos2x,vel2x,pos2y,vel2y,pos2z,vel2z, ...] )
#            output: vel & acc  (shape: [vel1x,acc1x,vel1y,acc1y,vel1z,acc1z, vel2x,acc2x,vel2y,acc2y,vel2z,acc2z, ...] )
# Solver has:                   (type : 4th order Runga-Kutta method, with adaptive step size)
#            input: vel & acc   (shape: [vel1x,acc1x,vel1y,acc1y,vel1z,acc1z, vel2x,acc2x,vel2y,acc2y,vel2z,acc2z, ...] )
#           output: pos & vel   (shape: [pos1x,vel1x,pos1y,vel1y,pos1z,vel1z, pos2x,vel2x,pos2y,vel2y,pos2z,vel2z, ...] )

start_time = time.time()  # define start time
wsol = solve_ivp(vectorfield, [0, stoptime], w_initial, method=solver_type,
                 atol=abserr, rtol=relerr)  # ,t_eval = t_eval_lst)
end_time = time.time()

## other approach
delta_t = 0.01
w_new = w_initial
i = 0
while i < stoptime:
    wsol = solve_ivp(vectorfield, [0, delta_t], w_new, method=solver_type,
                 atol=abserr, rtol=relerr)  # ,t_eval = t_eval_lst)
    i = i + delta_t
    w_new = wsol.y[:,-1]

# Solution

# Reading out the solution
w_step = wsol.y[:, -1]
w_step = w_new

# converting supervector w_step back to coordinate pos. & vel. vectors
pos_lst, vel_lst, i = [], [], 0  # Initializing
while i < (length_w):  # converting supervector w_step back to coordinate pos. & vel. vectors
    pos_lst.append([w_step[i], w_step[i+2], w_step[i+4]])  # positions
    vel_lst.append([w_step[i+1], w_step[i+3], w_step[i+5]])  # velocities
    i += 6  # Incrementing it per 3

pos_new, vel_new = np.array(pos_lst), np.array(
    vel_lst)  # Renaming for consistency
pos_ini[0][1] = L*0.1  # Removing the initial block offset
# Defining the max difference in position of the points
diff = np.amax(np.linalg.norm(pos_new-pos_ini, axis=1))

# printing out the solution
print("-------------------------------------------------------------")
print("-------------------------- NEW RUN --------------------------")
print("-------------------------------------------------------------")
print('Total time                   : ',
      round((end_time-start_time), 6), '[s]')
print('Evaluations required         : ', wsol.nfev, '    [-]')
print('Max difference in point pos. : ', np.round(diff, 6), '[m]')
print('pos_ini (y-coordinates)      : ', np.round(pos_ini[:, 1], 3))
print('pos_new (y-coordinates)      : ', np.round(pos_new[:, 1], 3))
#print('Evaluation time points       : ')
# print(np.round(wsol.t,3))
#print('Evaluation  points           : ')
# print(np.round(wsol.y,3))

# PLOTTING the solution
block_offset_for_print = np.round(block_offset, 3)
plt.plot(wsol.t, wsol.y[2]-L*0.1,
         label='block (offset = '+str(block_offset_for_print)+'m)')
for i in range(1, N-1):
    plt.plot(wsol.t, wsol.y[2+i*6] - (L*0.1+i*delta_L),
             label='particle '+str(i))
plt.xlabel('time /s')
plt.ylabel('Particle Deviations /m')
plt.grid()
plt.legend()
plt.show()

# %%
