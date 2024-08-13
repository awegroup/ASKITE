# -*- coding: utf-8 -*-
#%% Setup
import time
import logging

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.constants import c

from LightSailSim.particleSystem.ParticleSystem import ParticleSystem
from LightSailSim.Sim.simulations import Simulate_Lightsail
import LightSailSim.Mesh.mesh_functions as MF
import LightSailSim.ExternalForces.optical_interpolators.interpolators as interp
from LightSailSim.ExternalForces.LaserBeam import LaserBeam
from LightSailSim.ExternalForces.OpticalForceCalculator import OpticalForceCalculator
from LightSailSim.ExternalForces.OpticalForceCalculator import ParticleOpticalPropertyType

global_start_time = time.time()

# Setup parameters
params = {
    # model parameters
    "c": 1,  # [N s/m] damping coefficient
    "m_segment": 1, # [kg] mass of each node
    "thickness":100e-9, # [m] thickness of PhC

    # simulation settings
    "dt": 0.1,  # [s]       simulation timestep
    'adaptive_timestepping':1e-2, # [m] max distance traversed per timestep
    "t_steps": 1e6,  # [-]      max number of simulated time steps
    "abs_tol": 1e-10,  # [m/s]     absolute error tolerance iterative solver
    "rel_tol": 1e-5,  # [-]       relative error tolerance iterative solver
    "max_iter": int(1e2),  # [-]       maximum number of iterations for the bicgstab solver

    # Simulation Steps
    "convergence_threshold": 1e-6,
    "min_iterations":30,

    # Mesh_dependent_settings
    "midstrip_width": 1,
    "boundary_margin": 0.175
    }

params['E'] = 470e9
params['G'] = 0
params['E_x'] = params['E']*7/100
params['E_y'] = params['E']*18/100


# Setup mesh
n_segments = 15 # make sure this is uneven so there are no particles on the centerline
length = 1
mesh = MF.mesh_phc_square_cross(length,
                                mesh_edge_length=length/n_segments,
                                params = params,
                                noncompressive=True)
# We have to add some particles to act as a support structure.
stiffness_support = 1e9+1 # [n/m*m] line stiffness
k_support = stiffness_support / (length / n_segments)
l_support = length/n_segments/50

simulate_3D = False

for i in range((n_segments+1)**2):
    # calculate coordinates
    xyz = mesh[1][i][0].copy()
    if xyz[1] ==0 and simulate_3D:
        xyz[1]-=l_support
    elif xyz[1] == length and simulate_3D:
        xyz[1]+=l_support
    elif xyz[0] == 0:
        xyz[0]-=l_support
    elif xyz[0] == length:
        xyz[0]+=l_support


    if np.any(xyz != mesh[1][i][0]):
        particle = [xyz, np.zeros(3), params['m_segment'], True]
        link = [i, len(mesh[1]), k_support, 1]
        mesh[1].append(particle)
        mesh[0].append(link)

    xyz = mesh[1][i][0].copy()
    if (np.all(xyz == [0,0,0])
        or np.all(xyz == [0,length,0])
        or np.all(xyz == [length,0,0])
        or np.all(xyz == [length,length,0])) and simulate_3D:

        if xyz[0] ==0:
            xyz[0]-=l_support
        elif xyz[0] == length:
            xyz[0]+=l_support

        if np.any(xyz != mesh[1][i][0]):
            particle = [xyz, np.zeros(3), params['m_segment'], True]
            link = [i, len(mesh[1]), k_support, 1]
            mesh[1].append(particle)
            mesh[0].append(link)

# init particle system
PS = ParticleSystem(*mesh, params, clean_particles=False)
starting_postions = PS.x_v_current_3D[0]
# Setup the optical sytem
I_0 = 100e9 /(10*10)
mu_x = 0.5
mu_y = 0.5
sigma = 1/2
w=2*length
if simulate_3D:
    LB = LaserBeam(lambda x, y: I_0 * np.exp(-1/2 *((x-mu_x)/sigma)**2 # gaussian laser
                                             -1/2 *((y-mu_y)/sigma)**2),
                   lambda x,y: np.outer(np.ones(x.shape),[0,1]))
else:
    LB = LaserBeam(lambda x, y: I_0 * np.exp(-2*((x-mu_x)/w)**2), # gaussian laser
                   lambda x,y: np.outer(np.ones(x.shape),[0,1]))
# Import the crystal
fname = interp.PhC_library['Gao']
#fname = interp.PhC_library['dummy']
interp_right_side = interp.create_interpolator(fname,np.pi)
interp_left_side = interp.create_interpolator(fname, 0)


# set the correct boundary conditions and crystals on the particle system
for p in PS.particles:
    if simulate_3D:
        if p.x[1] == 0 or p.x[1] == length:
            p.set_fixed(True, [0,0,1], 'plane')

    if p.x[0] == 0 or p.x[0] == length:
        p.set_fixed(True, [0,0,1], 'plane')

    p.optical_type = ParticleOpticalPropertyType.ARBITRARY_PHC
    if p.x[0]>length/2:
        p.optical_interpolator = interp_right_side
    else:
        p.optical_interpolator = interp_left_side

OFC = OpticalForceCalculator(PS, LB)
SIM = Simulate_Lightsail(PS,OFC,params)

#%% Plot displaced PS with distributed and net forces

plot_check = True
deform = True

if plot_check:
    PS.displace([0,0,0,0,3,0])
    if deform:
        SIM.run_simulation(plotframes=0, printframes=50, simulation_function='kinetic_damping',file_id='_check_')

        fig_convergence = plt.figure()
        ax_kin = fig_convergence.add_subplot(211)
        ax_kin.semilogy(PS.history["E_kin"])
        ax_f = fig_convergence.add_subplot(212)
        ax_f.semilogy(PS.history["net_force"])
        fig_convergence.show()

    forces = OFC.force_value()

    net_force, net_moments = OFC.calculate_restoring_forces()
    fig = plt.figure(figsize = (10,6))

    ax = fig.add_subplot(projection='3d')
    ax = PS.plot(ax)

    COM = PS.calculate_center_of_mass()
    x,_ = PS.x_v_current_3D

    a_u = forces[:,0]
    a_v = forces[:,1]
    a_w = forces[:,2]

    x,y,z = x[:,0], x[:,1], x[:,2]

    ax.quiver(x,y,z,a_u,a_v,a_w, length = 5)
    ax.quiver(COM[0],COM[1],COM[2],
              net_force[0],net_force[1],net_force[2],
              length = 1, label ='Net Force', color='r')
    ax.quiver(COM[0],COM[1],COM[2],
              net_moments[0],net_moments[1],net_moments[2],
              length = 2, label ='Net Moment', color='magenta')
    fig.tight_layout()
    plt.show()

    # Very important to put it back where it came from!
    PS.un_displace()

rot_and_trans = False
if rot_and_trans:
    #%% Getting translation and rotation data
    print('Starting rotations and translations')
    translations = np.linspace(-length,length,17)
    rotations = np.linspace(-10,10,17)

    translation_plot = []
    trans_plot = False
    if trans_plot:
        fig0 = plt.figure(figsize = [20,16])

    # Deform the PS for each step?
    resimulate_on_displacement = True

    print("\n\nCalculating rotations")
    for i, t in enumerate(translations):
        print(f'\nTranslation {t=}')
        PS.displace([t,0,0,0,0,0])

        if resimulate_on_displacement:
            SIM.run_simulation(plotframes=0, printframes=50, simulation_function='kinetic_damping',file_id=f'_{t}_')

            # The force data is a little sensitive to random fluctuation, so instead I'm pulling the last
            # 50 entries from a ring buffer I added to the history.
            net_force = np.array([np.sum(forces,axis=0) for forces in PS.history['forces_ringbuffer']])
            net_force = np.sum(net_force, axis=0)/len(PS.history['forces_ringbuffer'])
            _, net_moments = OFC.calculate_restoring_forces()
        else:
            net_force, net_moments = OFC.calculate_restoring_forces()
        if trans_plot:
            fig0.clear()
            ax0 = fig0.add_subplot(projection='3d')
            ax0.set_xlim([-0.5,1.5])
            ax0.set_ylim([0,1])
            ax0.set_zlim([0,0.5])
            ax0.set_aspect('equal')
            ax0 = PS.plot_forces(OFC.force_value(),ax0)
            ax0.set_title(f'{t:.1f}')
            COM = PS.calculate_center_of_mass()
            ax0.quiver(COM[0],COM[1],COM[2],
                      net_force[0],net_force[1],net_force[2],
                      length = 1/2, label ='Net Force', color='r')
            ax0.quiver(COM[0],COM[1],COM[2],
                      net_moments[0],net_moments[1],net_moments[2],
                      length = 2, label ='Net Moment', color='magenta')
            ax0.quiver(COM[0],COM[1],COM[2],
                      net_moments[0],net_moments[1],net_moments[2],
                      length = 1.4, color='magenta')
            fig0.tight_layout()
            fig0.savefig(f'temp/translation-{i}-{t:.2f}.jpg', dpi = 200, format = 'jpg')
        PS.un_displace()

        translation_plot.append([t, *net_force, *net_moments])

    rotation_plot=[]
    rot_plot = False
    if rot_plot:
        fig0 = plt.figure(figsize = [20,16])

    print("\n\nCalculating rotations")
    for i, r in enumerate(rotations):
        print(f'\nRotation {r=}')
        PS.displace([0,0,0,0,r,0])

        if resimulate_on_displacement:
            SIM.run_simulation(plotframes=0, printframes=50, simulation_function='kinetic_damping', file_id=f'_{r}_')

            # The force data is a little sensitive to random fluctuation, so instead I'm pulling the last
            # 30 entries from a ring buffer I added to the history.
            net_force = np.array([np.sum(forces,axis=0) for forces in PS.history['forces_ringbuffer']])
            net_force = np.sum(net_force, axis=0)/len(PS.history['forces_ringbuffer'])
            _, net_moments = OFC.calculate_restoring_forces()
        else:
            net_force, net_moments = OFC.calculate_restoring_forces()

        if rot_plot:
            fig0.clear()
            ax0 = fig0.add_subplot(projection='3d')
            ax0 = PS.plot_forces(OFC.force_value(),ax0)
            ax0.set_title(f'{r}')
            ax0.quiver(COM[0],COM[1],COM[2],
                      net_force[0],net_force[1],net_force[2],
                      length = 0.5, label ='Net Force', color='r')
            ax0.quiver(COM[0],COM[1],COM[2],
                      net_moments[0],net_moments[1],net_moments[2],
                      length = 5, label ='Net Moment', color='magenta')
            ax0.quiver(COM[0],COM[1],COM[2],
                      net_moments[0],net_moments[1],net_moments[2],
                      length = 3.5, color='magenta')
            ax0.legend()
            ax0.set_xlim([0,1])
            ax0.set_ylim([0,1])
            ax0.set_zlim([-0.1,0.1])
            ax0.set_aspect('equal')
            fig0.tight_layout()
            fig0.savefig(f'temp/rotation-{i}-{r:.1f}.jpg', dpi = 200, format = 'jpg')

        PS.un_displace()

        rotation_plot.append([r, *net_force, *net_moments])


    translation_plot= np.array(translation_plot)
    rotation_plot = np.array(rotation_plot)
    header = ['displacement', 'Fx', 'Fy', 'Fz', 'Rx', 'Ry', 'Rz']
    pd.DataFrame(translation_plot,columns=header).to_csv(f"temp/translation_{stiffness_support=}.csv", header=True, index = False)
    pd.DataFrame(rotation_plot,columns=header).to_csv(f"temp/rotation_{stiffness_support=}.csv", header=True, index = False)



#%% Reproducing Fig. 4 from Gao et al 2022

gao_et_al_figure_four = False
if gao_et_al_figure_four:
    fig = plt.figure(figsize=(10,6))

    ax1 = fig.add_subplot(221)
    ax1.plot(rotation_plot[:,0], rotation_plot[:,3]/(I_0/c))
    ax1.set_title('Tilt angle versus vertical force')
    ax1.set_xlabel('Tilt angle [deg]')
    ax1.set_ylabel("$F_z [I_0D/c]$")
    ax1.set_ylim([0,rotation_plot[:,3].max()/(I_0/c)*1.2])
    ax1.set_xlim([-10,10])
    ax1.grid()

    ax2 = fig.add_subplot(222)
    ax2.plot(rotation_plot[:,0], rotation_plot[:,5]/(I_0/c))
    ax2.set_title('Tilt angle versus torque')
    ax2.set_xlabel('Tilt angle [deg]')
    ax2.set_ylabel("$\tau_y [I_0D^2/c]$")
    ax2.set_xlim([-10,10])
    ax2.grid()

    ax3 = fig.add_subplot(223)
    ax3.plot(rotation_plot[:,0], -rotation_plot[:,1]/(I_0/c))
    ax3.set_title('Tilt angle versus lateral force')
    ax3.set_xlabel('Tilt angle [deg]')
    ax3.set_ylabel("$F_x [I_0D/c]$")
    ax3.set_xlim([-10,10])
    ax3.grid()

    ax4 = fig.add_subplot(224)
    ax4.plot(translation_plot[:,0], translation_plot[:,1]/(I_0/c))
    ax4.set_title('Translation versus lateral force')
    ax4.set_xlabel('Translation [D]')
    ax4.set_ylabel("$F_x [I_0D/c]$")
    ax4.set_xlim([-1,1])
    ax4.grid()

    fig.tight_layout()


# %% Let's run some simulations!
simulate = False
if simulate:
    SIM.run_simulation(plotframes=0, printframes=10, simulation_function='kinetic_damping', plot_forces=True)
    PS.plot()


# %% print time it took
delta_time = time.time()-global_start_time
print(f'All in all that took {delta_time//60:.0f}m {delta_time%60:.2f}s')
