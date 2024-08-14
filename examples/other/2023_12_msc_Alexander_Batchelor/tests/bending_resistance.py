# %%

"""
Script for PS framework validation, benchmark case where tether is fixed at top end and exhibits longitudal oscillations
due to a dropped mass fixed at its other end at t = 0.

It plots out the runtime comparison between the different solvers
"""
%load_ext autoreload
%autoreload 2
import sys
import os
# TODO: remove this hardcoding
folder_path = '/home/jellepoland/surfdrive/phd/code/phd_Jelle_Poland/Simulations'
os.chdir(folder_path)  # This should not be needed
sys.path.append(os.getcwd())
                
import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt
import pandas as pd
import sys
import time
from src.particleSystem.ParticleSystem import ParticleSystem
from scipy.integrate import solve_ivp

def connectivity_matrix(n: int):
    matrix = np.eye(n, k=1) + np.eye(n, k=-1)
    return matrix

def initial_conditions(L: float, l0: float, n: int, m_segment: float, m_block: float, m_tether: float):
    
    m_segment = m_tether / (n - 1)
    
    ## Vertical Tether
    # conditions = [
    #     [[0, 0, (n - 1) * l0[i] - i * l0[i]], [0, 0, 0], m_segment, False]
    #     for i in range(n-1)
    # ]
    # conditions.append([[0,0,0],[0,0,0],m_block,False])

    ## Horizontal Tether
    conditions = [
        [[L + i * l0[i] - (n - 1) * l0[i], 0, 0], [0, 0, 0], m_segment, False]
        for i in range(n-1)
    ]
    conditions.append([[L, 0, 0], [0, 0, 0], m_block, False])


    conditions[0][-1] = True  # Set top end of tether to fixed
    conditions[-1][-1] = params['is_last_node_fixed'] # Set bottom end of tether to free
    conditions[0][-2] -= 0.5 * m_segment
    conditions[-1][-2] += m_block - 0.5 * m_segment
    return conditions

def initialize_tether_cross_over(k_bend,params):
 
    # calculated parameters
    params["l0"] = (params["L"] / (params["n"] - 1))*np.ones(params["n"]-1)
    params["m_segment"] = params["L"] * params["rho_tether"] / (params["n"] - 1)
    params["k"] = params["k"] * (params["n"] - 1)  # segment stiffness

    # instantiate connectivity matrix and initial conditions array
    b = np.nonzero(np.triu(connectivity_matrix(params["n"])))
    b_conn = np.column_stack((b[0], b[1]))

    init_cond = initial_conditions(
    params["L"],params["l0"], params["n"], params["m_segment"], params["m_block"], params["m_tether"]
    )

    params["is_fixed"] = np.array([init_cond[i][-1] for i in range(params["n"])])
    m = [init_cond[i][-2] for i in range(params["n"])]
    # params["c"] = 2 * np.sqrt(params["k"][0] / np.sum(m))

    #TODO: magic is happening here
    ## Violin's "crossover spring approach"
    if k_bend != 0:
        rotational_particles = {}
        for i in range(1,params["n"]-1):
            rotational_particles[f'{i}'] = [i-1,i+1]

        for key in rotational_particles.keys():
            # finding the particles nearby
            p1_idx,p2_idx = rotational_particles[key]
            p1 = np.array(init_cond[p1_idx][0])
            p2 = np.array(init_cond[p2_idx][0])

            # increasing the length of the variables
            b_conn = np.append(b_conn,[[p1_idx,p2_idx]],axis=0)
            params["k"] = np.append(params["k"],k_bend)
            params["l0"] = np.append(params["l0"],np.linalg.norm(p1 - p2))
            params["is_compression"] = np.append(params["is_compression"],True)
            params["is_tension"] = np.append(params["is_tension"],True)
            params["is_pulley"] = np.append(params["is_pulley"],False)
            params["is_rotational"] = np.append(params["is_rotational"],False)

    # print(f'params[l0]: {params["l0"]}')
    # print(f'b_conn: {b_conn}')
    # print(f'k: {params["k"]}')
    # print(f'l0: {params["l0"]}')

    return b_conn, init_cond, params

def initialize_tether_bending_spring(k_bend,params):

    # calculated parameters
    params["l0"] = (params["L"] / (params["n"] - 1))*np.ones(params["n"]-1)
    params["m_segment"] = params["L"] * params["rho_tether"] / (params["n"] - 1)
    params["k"] = params["k"] * (params["n"] - 1)  # segment stiffness

    # instantiate connectivity matrix and initial conditions array
    b = np.nonzero(np.triu(connectivity_matrix(params["n"])))
    b_conn = np.column_stack((b[0], b[1]))

    init_cond = initial_conditions(
    params["L"],params["l0"], params["n"], params["m_segment"], params["m_block"], params["m_tether"]
    )

    params["is_fixed"] = np.array([init_cond[i][-1] for i in range(params["n"])])
    m = [init_cond[i][-2] for i in range(params["n"])]
    # params["c"] = 2 * np.sqrt(params["k"][0] / np.sum(m))


    #TODO: magic is happening here
    ## "bending spring approach"
    if k_bend != 0:
        
        ## Identify indices of the last spring
        # to be able to add these as keys to the new bending_params  
        idx_last_normal_spring = len(params["k"]) - 1
        idx_of_this_spring = idx_last_normal_spring
        bending_params = {} 

        ## Define a rotational particle list
        # each node with resistance is an index/key
        # each value represents the links 
        rotational_particles = {}
        for i in range(1,params["n"]-1):
            rotational_particles[f'{i}'] = [i-1,i+1]

        ## Looping over each key, and defining the surroundings
        # print(f'rotational_particles: {len(rotational_particles)} {rotational_particles}')
        for key in rotational_particles.keys():
            # finding the particles nearby
            p1_idx,p2_idx = rotational_particles[key]
            pe_idx = int(key)
            p1 = np.array(init_cond[p1_idx][0])#extracting only the initial position 
            p2 = np.array(init_cond[p2_idx][0])#extracting only the initial position

            # defining the particles
            Pa = p1
            Pb = p2
            Pe = np.array(init_cond[int(key)][0])
            # print(f'Pa: {p1_idx} - {Pa}')
            # print(f'Pb: {p2_idx} - {Pb}')
            # print(f'Pe: {pe_idx} - {Pe}')

            # calculating the projection vector
            P_proj = Pa + np.dot(Pe - Pa, Pb - Pa) / np.dot(Pb - Pa, Pb - Pa) * (
                                Pb - Pa
                            )
            # distance from A and B to the projected point
            h_a = np.linalg.norm(P_proj - Pa)
            h_b = np.linalg.norm(P_proj - Pb)

            # defining the alpha's 
            alpha_a = h_b / (h_a + h_b)
            alpha_b = h_a / (h_a + h_b)
            alpha_e = -1

            # making a dict of the bending_parameters
            # print(f'idx_of_this_spring: {idx_of_this_spring}')
            bending_params[f'{idx_of_this_spring}'] = [alpha_a,alpha_b,alpha_e,pe_idx] 
            idx_of_this_spring += 1

            
            #TODO: this essentially reduces k_bend to mu by a factor of 10
            thickness_of_rotational_resistance_beam = 0.1 #represents the thickness of the strut tubes 
            ## calculating lambda_bend, fed in as k[idx] for the rotational True springs
            lambda_bend = (2/3)* ((h_a+h_b)/((h_a*h_b)**2)) * k_bend* thickness_of_rotational_resistance_beam

            # increasing the length of the variables, by 1 - done for each additional rotational spring element 
            b_conn = np.append(b_conn,[[p1_idx,p2_idx]],axis=0)
            params["k"] = np.append(params["k"],lambda_bend)
            params["l0"] = np.append(params["l0"],np.linalg.norm(p1 - p2))
            params["is_compression"] = np.append(params["is_compression"],False)
            params["is_tension"] = np.append(params["is_tension"],False)
            params["is_pulley"] = np.append(params["is_pulley"],False)
            params["is_rotational"] = np.append(params["is_rotational"],True)

        # adding the alpha_dict to params
        # print(f'params[l0]: {params["l0"]}') 
        params['bending_params'] = bending_params       
        print(f'params:{params}')
        print(f'bending_params: {bending_params}')
        # print(f'b_conn: {b_conn}')

    return b_conn, init_cond, params


def plot(psystem, init_cond,plt_number,f_shear,k_bend,params,plot_title):
    n = params["n"]
    t_vector = np.linspace(0, params["t_steps"] * params["dt"], params["t_steps"] + 1)

    x = {}
    v = {}
    for i in range(n):
        x[f"x{i + 1}"] = np.zeros(len(t_vector))
        x[f"y{i + 1}"] = np.zeros(len(t_vector))
        x[f"z{i + 1}"] = np.zeros(len(t_vector))
        v[f"vx{i + 1}"] = np.zeros(len(t_vector))
        v[f"vy{i + 1}"] = np.zeros(len(t_vector))
        v[f"vz{i + 1}"] = np.zeros(len(t_vector))

    position = pd.DataFrame(index=t_vector, columns=x)
    velocity = pd.DataFrame(index=t_vector, columns=v)

    g = params["g"]
    n = params["n"]

    m = [init_cond[i][-2] for i in range(n)]
    f_ext = np.array([np.array([0, 0, -g * m[i]])+f_shear[i] for i in range(n)])
    f_ext = f_ext.flatten()
    masses = np.array([[init_cond[i][-2] for j in range(3)] for i in range(n)])
    masses = masses.flatten()
    m = np.diag(masses)

    start_time = time.time()
    for step in t_vector:  # propagating the simulation for each timestep and saving results
        if step == 0:
            x, v = psystem.x_v_current
            position.loc[step], velocity.loc[step] = x, v

        x_next, v_next = psystem.kin_damp_sim(f_ext)

        position.loc[step], velocity.loc[step] = x_next, v_next

        residual_f_with_fixed_nodes = f_ext + psystem.f_int
        residual_f = []
        for is_fixed in params["is_fixed"]:
            if not is_fixed:
                residual_f.append(residual_f_with_fixed_nodes)

        if np.linalg.norm(residual_f) <= 1e-3:
            print(f"Convergences at time: {step:.2f}")
            print(step)
            break
    
    stop_time = time.time()
    # print(f'residual with fixed nodes: {np.linalg.norm(residual_f_with_fixed_nodes)}')
    print(f'{plot_title} runtime: {stop_time-start_time:.2f}s, f_residual: {np.linalg.norm(residual_f):.2f}N')
    points_x = [point[0] for point in x_next.reshape((params["n"], 3))]
    points_y = [point[1] for point in x_next.reshape((params["n"], 3))]
    points_z = [point[2] for point in x_next.reshape((params["n"], 3))]

    plt.subplot(3,1,plt_number)
    plt.plot(points_x,points_z,'o-',label=f'{plot_title} k_bend: {k_bend}')
    plt.legend(loc=[0.65,0.05])
    # plt.title(f'{plot_title} springs')
    plt.xlim(-0.05*params['xlim'],params["xlim"])
    plt.ylim(-0.5*params["ylim"],0.1*params["ylim"])


def plot_all(params,n_nodes,k_bend_list,f_shear,is_last_node_fixed):


    params["is_last_node_fixed"] = is_last_node_fixed

    # plotting
    fig = plt.figure(figsize=(12, 7))

    k_bend = k_bend_list[0]
    params1 = params.copy()
    b_conn1, init_cond1, params1 = initialize_tether_cross_over(k_bend,params1)
    ps1 = ParticleSystem(b_conn1, init_cond1,params1)
    plot(ps1, init_cond1,1,f_shear,k_bend,params1,'Crossover')

    k_bend = k_bend_list[1]
    params2 = params.copy()
    b_conn2, init_cond2, params2 = initialize_tether_cross_over(k_bend,params2)
    ps2 = ParticleSystem(b_conn2, init_cond2,params2)
    plot(ps2, init_cond2,1,f_shear,k_bend,params2,'Crossover')

    k_bend = k_bend_list[2]
    params3 = params.copy()
    b_conn3, init_cond3, params3 = initialize_tether_cross_over(k_bend,params3)
    ps3 = ParticleSystem(b_conn3, init_cond3,params3)
    plot(ps3, init_cond3,1,f_shear,k_bend,params3,'Crossover')

    k_bend = k_bend_list[3]
    params4 = params.copy()
    b_conn4, init_cond4, params4 = initialize_tether_cross_over(k_bend,params4)
    ps4 = ParticleSystem(b_conn4, init_cond4,params4)
    plot(ps4, init_cond4,1,f_shear,k_bend,params4,'Crossover')


    ############################

    k_bend = k_bend_list[0]
    params5 = params.copy()
    b_conn5, init_cond5, params5 = initialize_tether_bending_spring(k_bend,params5)
    ps5 = ParticleSystem(b_conn5, init_cond5,params5)
    plot(ps5, init_cond5,2,f_shear,k_bend,params5,'Bending')

    k_bend = k_bend_list[1]
    params6 = params.copy()
    b_conn6, init_cond6, params6 = initialize_tether_bending_spring(k_bend,params6)
    ps6 = ParticleSystem(b_conn6, init_cond6,params6)
    plot(ps6, init_cond6,2,f_shear,k_bend,params6,'Bending')

    k_bend = k_bend_list[2]
    params7 = params.copy()
    b_conn7, init_cond7, params7 = initialize_tether_bending_spring(k_bend,params7)
    ps7 = ParticleSystem(b_conn7, init_cond7,params7)
    plot(ps7, init_cond7,2,f_shear,k_bend,params7,'Bending')

    k_bend = k_bend_list[3]
    params8 = params.copy()
    b_conn8, init_cond8, params8 = initialize_tether_bending_spring(k_bend,params8)
    ps8 = ParticleSystem(b_conn8, init_cond8,params8)
    plot(ps8, init_cond8,2,f_shear,k_bend,params8,'Bending')

    ############################

    k_bend = k_bend_list[0]
    params9 = params.copy()
    b_conn9, init_cond9, params9 = initialize_tether_cross_over(k_bend,params9)
    ps9 = ParticleSystem(b_conn9, init_cond9,params9)
    plot(ps9, init_cond9,3,f_shear,k_bend,params9,'Torsion')

    # k_bend = k_bend_list[1]
    # params10 = params.copy()
    # b_conn10, init_cond10, params10 = initialize_tether_bending_spring(k_bend,params10)
    # ps10 = ParticleSystem(b_conn10, init_cond10,params10)
    # plot(ps10, init_cond10,3,f_shear,k_bend,params10,'Torsion')

    # k_bend = k_bend_list[2]
    # params11 = params.copy()
    # b_conn11, init_cond11, params11 = initialize_tether_bending_spring(k_bend,params11)
    # ps11 = ParticleSystem(b_conn11, init_cond11,params11)
    # plot(ps11, init_cond11,3,f_shear,k_bend,params11,'Torsion')

    # k_bend = k_bend_list[3]
    # params12 = params.copy()
    # b_conn12, init_cond12, params12 = initialize_tether_bending_spring(k_bend,params12)
    # ps12 = ParticleSystem(b_conn12, init_cond12,params12)
    # plot(ps12, init_cond12,3,f_shear,k_bend,params12,'Torsion')
    

n_nodes = 15
k = 1e2
mag_f_shear = 7e1 # 2e2
xlim = 21.
ylim = 35.
is_last_node_fixed = False

f_shear = np.zeros((n_nodes,3))

## applying it equally over each node
# f_shear = np.array([[mag_f_shear/n_nodes,0,0] for i in range(n_nodes)])
# f_shear[-2:-1] *= 0
# f_shear[0:2] *= 0

## applying it only to the last-node upwards
f_shear[(n_nodes-1)][2] = mag_f_shear
# f_shear[(n_nodes-2)][2] = mag_f_shear/2
# f_shear[(n_nodes-3)][2] = mag_f_shear/3


# dictionary of required parameters
params = {
    # model parameters
    "n": n_nodes,  # [-] number of particles
    "k": k*np.ones(n_nodes),  # [N/m] spring stiffness
    "c": 1e3,  # [N s/m] damping coefficient
    "L": 20,  # [m] tether length
    "m_block": 1.,  # [kg] mass attached to end of tether
    "m_tether": 1.*n_nodes,
    "rho_tether": .1,# 0.1,  # [kg/m] mass density tether
    # simulation settings
    "dt": 0.05,  # [s] simulation timestep
    "t_steps": int(1e2),  # [-] number of simulated time steps
    "abs_tol": 1e-5,  # [m/s] absolute error tolerance iterative solver
    "rel_tol": 1e-4,  # [-] relative error tolerance iterative solver
    "max_iter": 1e5,  # [-] maximum number of iterations
    # physical parameters
    "g": 9.807,  # [m/s^2] gravitational acceleration
    "is_compression": np.array([True for i in range(n_nodes)]),
    "is_tension": np.array([True for i in range(n_nodes)]),
    "is_pulley": np.array([False for i in range(n_nodes)]),
    "is_rotational": np.array([False for i in range(n_nodes)]),
    "xlim": xlim,
    "ylim": ylim}

# plot_all(params,n_nodes,[0,1e2,1e3,5e4],f_shear,is_last_node_fixed=True)
plot_all(params,n_nodes,[1e5,1e6,1e8,1e10],f_shear,is_last_node_fixed=False)
plt.show()

























#%% 
## trying to implement simple linear bending spring

# defining the input
Pa = np.array([0.5,0,2.0])
Pb = np.array([0.0,0,0.0])
Pe = np.array([0.2,0,1.5])
thickness_of_rotational_resistance_beam = 1. # diameter of tube, sort of width of the cell
mu_bending = 1e2 # bending stiffness of the tube

### a bit longer, below is in 1-line
# # Find the direction vector of the line from Pa to Pb
# D = Pb - Pa

# # Find the vector from Pa to Pe
# V = Pe - Pa

# # Calculate the projection
# projection = np.dot(V, D) / np.dot(D, D) * D

# # Add the projection to point Pa to get the projected point
# P_proj = Pa + projection


P_proj = Pa + np.dot(Pe - Pa, Pb - Pa) / np.dot(Pb - Pa, Pb - Pa) * (
                    Pb - Pa
                )
# distance from A and B to the projected point
h_a = np.linalg.norm(P_proj - Pa)
h_b = np.linalg.norm(P_proj - Pb)

# defining the alpha's
## they will vary depending on the 
alpha_a = h_b / (h_a + h_b)
alpha_b = h_a / (h_a + h_b)
alpha_e = -1

# # Calculate the normals for edges (A-B) and (A-E)
# N_AB = np.cross(Pb - Pa, Pe - Pa)
# N_AE = np.cross(Pe - Pa, Pb - Pa)

# # Calculate the alpha coefficients
# alpha_a = np.linalg.norm(N_AE) / (np.linalg.norm(N_AE) + np.linalg.norm(N_AB))
# alpha_b = np.linalg.norm(N_AB) / (np.linalg.norm(N_AE) + np.linalg.norm(N_AB))
# alpha_e = -0.5

# defining lambda
K_bend = mu_bending * thickness_of_rotational_resistance_beam
lambda_bend = (2/3)* ((h_a+h_b)/((h_a*h_b)**2)) * K_bend

plt.figure()
# plt.scatter(P_proj[0],P_proj[2],label='Projected point',c='r')
plt.plot([Pa[0],Pb[0]],[Pa[2],Pb[2]],'-',label='Pa to Pb',c='k')
plt.plot([Pa[0],Pe[0],Pb[0]],[Pa[2],Pe[2],Pb[2]],'o--','blue',label='Original line')

## at runtime
Pa_new = np.array([0.5,0,2.0])
Pb_new = np.array([0.0,0,0.0])
Pe_new = np.array([1.0,0,1.0])

# calculate bending vector R

# old definitino doesn't work
# R = alpha_a * Pa_new + alpha_b * Pb_new + alpha_e * Pe_new

# new wrt original orientation
R = Pe - Pe_new

# new wrt being a straight line
P_proj = Pa_new + np.dot(Pe_new - Pa_new, Pb_new - Pa_new) / np.dot(Pb_new - Pa_new, Pb_new - Pa_new) * (Pb_new - Pa_new)
R = P_proj - Pe_new

# calculate the forces
Fa = -lambda_bend * alpha_a * R
Fb = -lambda_bend * alpha_b * R
Fe = -lambda_bend * alpha_e * R

# printing
print(f'alpha_a: {alpha_a}')
print(f'alpha_b: {alpha_b}')
print(f'alpha_e: {alpha_e}')
print(f'lambda_bend: {lambda_bend}')
print(f'R: {R}')
print(f'F_a: {Fa}')
print(f'F_b: {Fb}')
print(f'F_e: {Fe}')

# plotting
plt.title('Bending spring')

# Plot the vectors as arrows
Fa_plot = 0.004*Fa
Fb_plot = 0.004*Fb
Fe_plot = 0.004*Fe 
plt.quiver(Pa_new[0], Pa_new[2], Fa_plot[0], Fa_plot[2], angles='xy', scale_units='xy', scale=1, color='r', label='F_a')
plt.quiver(Pb_new[0], Pb_new[2], Fb_plot[0], Fb_plot[2], angles='xy', scale_units='xy', scale=1, color='g', label='F_b')
plt.quiver(Pe_new[0], Pe_new[2], Fe_plot[0], Fe_plot[2], angles='xy', scale_units='xy', scale=1, color='b', label='F_e')


plt.plot([Pe_new[0],(Pe_new[0]+R[0])],[Pe_new[2],(Pe_new[2]+R[2])],'--',c='g',label='R',lw = 2, ms = 8)
plt.plot([Pa_new[0],Pe_new[0],Pb_new[0]],[Pa_new[2],Pe_new[2],Pb_new[2]],'o-',c='orange',label='Bent line')
plt.text(Pa_new[0],Pa_new[2],'Pa',fontsize=15)
plt.text(Pb_new[0],Pb_new[2],'Pb',fontsize=15)
plt.text(Pe_new[0],Pe_new[2],'Pe',fontsize=15)
plt.legend()
# plt.axis('equal')
plt.xlim(-0.5,3)
plt.ylim(-1,2.5)
plt.show()


# # %% Trying with 4 points

# import numpy as np
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# from mpl_toolkits.mplot3d.art3d import Poly3DCollection
# %matplotlib widget

# # Define the positions of the 4 points
# Pa = np.array([0.0, 0, 1.0])
# Pb = np.array([2.0, 0, 1.0])
# Pc = np.array([1.0, -0.5, .8])
# Pd = np.array([1.0, 0.5, .8])

# thickness_of_rotational_resistance_beam = 1.0  # Diameter of the tube, sort of width of the cell
# mu_bending = 1e2  # Bending stiffness of the tube

# # Find the direction vector of the line from Pa to Pb
# D = Pb - Pa

# # Find the vector from Pa to Pc
# V = Pc - Pa

# # Calculate the projection
# projection = np.dot(V, D) / np.dot(D, D) * D

# # Add the projection to point Pa to get the projected point
# P_proj = Pa + projection

# # Distance from A and B to the projected point
# h_a = np.linalg.norm(P_proj - Pa)
# h_b = np.linalg.norm(P_proj - Pb)

# # Calculate the normals for edges (A-B) and (A-C)
# N_AB = np.cross(Pb - Pa, Pc - Pa)
# N_AC = np.cross(Pc - Pa, Pb - Pa)

# # Calculate the alpha coefficients
# alpha_a = h_b / (h_a + h_b)
# alpha_b = h_a / (h_a + h_b)
# alpha_c = np.linalg.norm(N_AC) / (np.linalg.norm(N_AC) + np.linalg.norm(N_AB))
# alpha_d = np.linalg.norm(N_AB) / (np.linalg.norm(N_AC) + np.linalg.norm(N_AB))

# # Define the lambda coefficient
# K_bend = mu_bending * thickness_of_rotational_resistance_beam
# lambda_bend = (2/3) * ((h_a + h_b) / (h_a * h_b)**2) * K_bend

# # Calculate the bending vector R
# R = alpha_a * Pa + alpha_b * Pb + alpha_c * Pc + alpha_d * Pd

# # Calculate the forces
# Fa = -lambda_bend * alpha_a * R
# Fb = -lambda_bend * alpha_b * R
# Fc = -lambda_bend * alpha_c * R
# Fd = -lambda_bend * alpha_d * R

# # Print the results
# print(f'alpha_a: {alpha_a}')
# print(f'alpha_b: {alpha_b}')
# print(f'alpha_c: {alpha_c}')
# print(f'alpha_d: {alpha_d}')
# print(f'lambda_bend: {lambda_bend}')
# print(f'R: {R}')
# print(f'F_a: {Fa}')
# print(f'F_b: {Fb}')
# print(f'F_c: {Fc}')
# print(f'F_d: {Fd}')

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# # Define the vertices for the two triangular surfaces
# vertices1 = [Pa, Pc, Pd]
# vertices2 = [Pb, Pc, Pd]

# # Define the triangles for the two surfaces
# triangles1 = [[0, 1, 2]]
# triangles2 = [[0, 1, 2]]

# # Plot the first surface with color
# surface1 = Poly3DCollection([vertices1], facecolors='orange', edgecolors='k', alpha=0.7)
# ax.add_collection3d(surface1)

# # Plot the second surface with color
# surface2 = Poly3DCollection([vertices2], facecolors='orange', edgecolors='k', alpha=0.7)
# ax.add_collection3d(surface2)

# # # Plot the lines
# # ax.plot([Pa[0], Pb[0], Pc[0], Pd[0], Pa[0]], [Pa[1], Pb[1], Pc[1], Pd[1], Pa[1]],
# #         [Pa[2], Pb[2], Pc[2], Pd[2], Pa[2]], 'k-')

# # Plot the bending vector R
# ax.quiver(R[0], R[1], R[2], R[0], R[1], R[2], color='g', label='R', arrow_length_ratio=0.1)


# # Labels for points
# ax.text(Pa[0], Pa[1], Pa[2], 'Pa', fontsize=12, color='k', label='Point Labels')
# ax.text(Pb[0], Pb[1], Pb[2], 'Pb', fontsize=12, color='k')
# ax.text(Pc[0], Pc[1], Pc[2], 'Pc', fontsize=12, color='k')
# ax.text(Pd[0], Pd[1], Pd[2], 'Pd', fontsize=12, color='k')

# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.set_zlabel('Z')
# ax.legend()

# plt.show()

# %%
