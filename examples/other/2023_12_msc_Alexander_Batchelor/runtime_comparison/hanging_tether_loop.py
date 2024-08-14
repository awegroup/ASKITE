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

n_test_list = [3, 4, 10]

t_visc_list, t_kin_damp_list, t_kin_damp_with_q_list, t_solve_ivp_list = [], [], [], []
for num_nodes in n_test_list:
    print(f"----------------------------------------------")
    print(f" Starting New Loop with {num_nodes} nodes ")
    print(f"----------------------------------------------")
    print(f" ")

    # dictionary of required parameters
    params = {
        # model parameters
        "n": num_nodes,  # [-] number of particles
        "k": 1e3*np.ones(num_nodes),  # [N/m] spring stiffness
        "c": 20,  # [N s/m] damping coefficient
        "L": 10,  # [m] tether length
        "m_block": 100,  # [kg] mass attached to end of tether
        "rho_tether": 0.1,  # [kg/m] mass density tether
        # simulation settings
        "dt": 0.05,  # [s] simulation timestep
        "t_steps": int(2.5e3),  # [-] number of simulated time steps
        "abs_tol": 1e-5,  # [m/s] absolute error tolerance iterative solver
        "rel_tol": 1e-4,  # [-] relative error tolerance iterative solver
        "max_iter": 1e5,  # [-] maximum number of iterations
        # physical parameters
        "g": 9.807,  # [m/s^2] gravitational acceleration
        "is_compression": np.array([True for i in range(num_nodes)]),
        "is_tension": np.array([True for i in range(num_nodes)]),
        "is_pulley": np.array([False for i in range(num_nodes)]),
        "is_rotational": np.array([False for i in range(num_nodes)]),
        "is_fixed": np.array([False for i in range(num_nodes)]),
        "pulley_other_line_pair": {},
    }

    ## reducing number of t_steps, when n-increases
    if params["n"] > 8:
        params["t_steps"] = int(params["t_steps"] / (2))

    def connectivity_matrix(n: int):
        matrix = np.eye(n, k=1) + np.eye(n, k=-1)
        return matrix

    def initial_conditions(l0: float, n: int, m_segment: float, m_block: float):
        conditions = [
            [[0, 0, (n - 1) * l0[i] - i * l0[i]], [0, 0, 0], m_segment, False]
            for i in range(n)
        ]
        conditions[0][-1] = True  # Set top end of tether to fixed
        conditions[0][-2] -= 0.5 * m_segment
        conditions[-1][-2] += m_block - 0.5 * m_segment
        return conditions

    # calculated parameters
    params["l0"] = (params["L"] / (params["n"] - 1))*np.ones(params["n"])
    params["m_segment"] = params["L"] * params["rho_tether"] / (params["n"] - 1)
    params["k"] = params["k"] * (params["n"] - 1)  # segment stiffness

    # instantiate connectivity matrix and initial conditions array
    b = np.nonzero(np.triu(connectivity_matrix(params["n"])))
    b_conn = np.column_stack((b[0], b[1]))

    init_cond = initial_conditions(
        params["l0"], params["n"], params["m_segment"], params["m_block"]
    )

    m = [init_cond[i][-2] for i in range(params["n"])]
    params["c"] = 2 * np.sqrt(params["k"][0] / np.sum(m))

    def instantiate_ps():
        return ParticleSystem(b_conn, init_cond, params)

    ps1, ps2, ps3 = instantiate_ps(), instantiate_ps(), instantiate_ps()

    def system_for_solve_ivp(t, state, args_system):
        n, masses, grav_constant, k, c, l0, b_conn = args_system

        positions = state[:n]
        velocities = state[n:]

        def calculate_f_spring(p1, p2, k, l0):
            relative_pos = np.array([p1 - p2])
            norm_pos = np.linalg.norm(relative_pos)

            if norm_pos != 0:
                unit_vector = relative_pos / norm_pos
            else:
                unit_vector = np.array([0, 0, 0])

            delta_length = norm_pos - l0
            f_spring = -k * delta_length * unit_vector

            return np.squeeze(f_spring)

        f_int = np.zeros_like(positions)
        for idx,conn in enumerate(b_conn):
            i, j = conn

            p1 = np.array([0, 0, positions[i][0]])
            p2 = np.array([0, 0, positions[j][0]])

            fs = calculate_f_spring(p1, p2, k[idx], l0[idx])

            f_int[i] += fs[2]
            f_int[j] -= fs[2]

        accelerations = np.zeros_like(positions)
        for i, (pos, vel, mass) in enumerate(zip(positions, velocities, masses)):
            ##TODO: flatten points, forces and stuff
            f_s = f_int[i]  # [i*3: i*3+3]
            f_g = -grav_constant * mass
            f_d = -c * vel
            acc = (f_s + f_g + f_d) / mass

            if i == 0:
                acc = [0]  # np.zeros(3)

            accelerations[i] = acc

        velocities = velocities.flatten()
        accelerations = accelerations.flatten()
        return np.concatenate([velocities, accelerations])

    def plot(
        psystem: ParticleSystem, psystem2: ParticleSystem, psystem3: ParticleSystem
    ):
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
        position2 = pd.DataFrame(index=t_vector, columns=x)
        velocity2 = pd.DataFrame(index=t_vector, columns=v)
        position3 = pd.DataFrame(index=t_vector, columns=x)
        velocity3 = pd.DataFrame(index=t_vector, columns=v)

        g = params["g"]
        n = params["n"]

        m = [init_cond[i][-2] for i in range(n)]
        f_ext = np.array([[0, 0, -g * m[i]] for i in range(n)])
        f_ext = f_ext.flatten()
        f_check = f_ext.copy()
        f_check[:3] = 0
        masses = np.array([[init_cond[i][-2] for j in range(3)] for i in range(n)])
        masses = masses.flatten()
        m = np.diag(masses)

        start_time = time.time()
        for step in t_vector:  # propagating the simulation for each timestep and saving results
            if step == 0:
                x, v = psystem.x_v_current
                position.loc[step], velocity.loc[step] = x, v

            x_next, v_next = psystem.simulate(f_ext)

            position.loc[step], velocity.loc[step] = x_next, v_next

            residual_f = f_check + psystem.f_int
            if np.linalg.norm(residual_f) <= 1e-3:
                print(f"viscous convergences at time: {step:.2f}")
                break
        stop_time = time.time()
        t_visc = stop_time - start_time

        start_time2 = time.time()
        for step in t_vector:  # propagating the simulation for each timestep and saving results

            if step == 0:
                x, v = psystem2.x_v_current
                position2.loc[step], velocity2.loc[step] = x, v
                continue

            x_next, v_next = psystem2.kin_damp_sim(f_ext)
            position2.loc[step], velocity2.loc[step] = x_next, v_next

            residual_f = f_check + psystem2.f_int
            if np.linalg.norm(residual_f) <= 1e-3:
                print(f"kin_damp convergences at time: {step:.2f}")
                break
        stop_time2 = time.time()
        t_kin_damp = stop_time2 - start_time2

        start_time3 = time.time()
        for step in t_vector:  # propagating the simulation for each timestep and saving results

            if step == 0:
                x, v = psystem3.x_v_current
                position3.loc[step], velocity3.loc[step] = x, v
                continue

            x_next, v_next = psystem3.kin_damp_sim(f_ext, q_correction=True)
            position3.loc[step], velocity3.loc[step] = x_next, v_next
            residual_f = f_check + psystem3.f_int
            if np.linalg.norm(residual_f) <= 1e-3:
                print(f"kin_damp with q convergences at time: {step:.2f}")
                break
        stop_time3 = time.time()
        t_kin_damp_with_q = stop_time3 - start_time3

        ### SOLVE IVP
        # initialising the system
        mass_input = np.array([mass[0] for mass in masses.reshape(n, 3)])
        args_system = [n, mass_input, g, params["k"], params["c"], params["l0"], b_conn]
        t_span = (0, params["t_steps"] * params["dt"])
        t_eval = t_vector

        # finding the initial state
        pos_1d = np.array([cond[0][2] for cond in init_cond])
        vel_1d = np.zeros(n)
        initial_state = np.concatenate([pos_1d, vel_1d])

        # initial run to determine when it convergences
        solution = solve_ivp(
            system_for_solve_ivp,
            t_span,
            initial_state,
            t_eval=t_eval,
            vectorized=True,
            args=(args_system,))
        
        position4 = solution.y
        convergence_step = len(t_vector)-1
        for i, y in enumerate(position4[n - 1]):
            if np.abs(y - x_next[-1]) < 1e-3:
                print(f"solve_ivp converges at time: {i:.2f}")
                convergence_step = i
                break

        ## adjusting the t_span to the point at which the solve_ivp converges (could be done better ofcourse)
        t_span = (0, convergence_step * params["dt"])
        t_eval = t_vector[: convergence_step + 1]
        start_time4 = time.time()
        solution = solve_ivp(
            system_for_solve_ivp,
            t_span,
            initial_state,
            t_eval=t_eval,
            vectorized=True,
            args=(args_system,))
        
        stop_time4 = time.time()
        t_solve_ivp = stop_time4 - start_time4

        print(f"PS classic: {(stop_time - start_time):.4f} s")
        print(f"PS kinetic w/o q: {(stop_time2 - start_time2):.4f} s")
        print(f"PS kinetic with q: {(stop_time3 - start_time3):.4f} s")
        print(f"Solve_ivp: {(stop_time4 - start_time4):.4f} s")

        ## correction the position2 and position3 delay
        position2 = position2.shift(-1)
        position3 = position3.shift(-1)

        # generating analytical solution for the same time vector
        # exact, decay = exact_solution(t_vector)

        # plotting & graph configuration
        # for i in range(1, n):
        #     position[f"z{i + 1}"] -=   init_cond[i][0][-1]
        #     position[f"z{i + 1}"].plot()
        #     position2[f"z{i + 1}"] -=   init_cond[i][0][-1]
        #     position2[f"z{i + 1}"].plot()
        #     position3[f"z{i + 1}"] -=   init_cond[i][0][-1]
        #     position3[f"z{i + 1}"].plot()

        ## below are the plot functions
        i = n - 1
        position[f"z{i + 1}"] -= init_cond[i][0][-1]
        position[f"z{i + 1}"].plot()
        position2[f"z{i + 1}"] -= init_cond[i][0][-1]
        position2[f"z{i + 1}"].plot(lw=3)
        position3[f"z{i + 1}"] -= init_cond[i][0][-1]
        position3[f"z{i + 1}"].plot()

        plt.plot(t_vector, position4[n - 1], lw=1)

        # x_length = position["z2"].count()
        # plt.plot(t_vector[:x_length], exact.iloc[:x_length], ls='--')

        # for i in range(n - 1):      # setting particle colors equivalent to their analytical solution
        #     color = plt.gca().lines[i].get_color()
        #     plt.gca().lines[i + n - 1].set_color(color)

        plt.xlabel("time [s]")
        plt.ylabel("position of last node [m]")
        plt.title(" Simulation of longitudal tether oscillations")
        plt.legend(["viscous", "kinetic w/o q", "kinetic with q", "solve_ivp"])
        plt.grid()
        plt.show()

        t_visc_list.append(t_visc)
        t_kin_damp_list.append(t_kin_damp)
        t_kin_damp_with_q_list.append(t_kin_damp_with_q)
        t_solve_ivp_list.append(t_solve_ivp)

        return t_visc_list, t_kin_damp_list, t_kin_damp_with_q_list, t_solve_ivp_list

    t_visc, t_kin_damp, t_kin_damp_with_q, t_solve_ivp = plot(ps1, ps2, ps3)


def plotting_results(
    t_visc_list, t_kin_damp_list, t_kin_damp_with_q_list, t_solve_ivp_list
):
    plt.figure()
    plt.plot(n_test_list, t_solve_ivp_list, "o-")
    plt.plot(n_test_list, t_visc_list, "o-")
    plt.plot(n_test_list, t_kin_damp_list, "o-")
    plt.plot(n_test_list, t_kin_damp_with_q_list, "o-")
    plt.legend([
            "Black-box: scipy.integrate.solve_ivp",
            "Own Solver: Viscous",
            "Own Solver: Kinetic w/o q",
            "Own Solver: Kinetic with q",
                ])
    plt.xlabel("Number of nodes")
    plt.ylabel("Runtime [s]")
    plt.title("Runtime Comparison using a hanging tether")
    plt.grid()
    plt.show()


if __name__ == "__main__":
    plotting_results(t_visc_list, t_kin_damp_list, t_kin_damp_with_q_list, t_solve_ivp_list)

# %%
