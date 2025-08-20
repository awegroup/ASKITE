"""
ParticleSystem framework
...
"""

import copy
import numpy as np
import numpy.typing as npt
from scipy.sparse.linalg import bicgstab
from kitesim.particleSystem.Particle import Particle
from kitesim.particleSystem.SpringDamper import SpringDamper


class ParticleSystem:
    def __init__(
        self,
        connectivity_matrix: npt.ArrayLike,
        initial_conditions: npt.ArrayLike,
        sim_param: dict,
    ):
        """
        Constructor for ParticleSystem object, model made up of n particles
        :param connectivity_matrix: sparse n-by-n matrix, where an 1 at index (i,j) means
                                    that particle i and j are connected
        :param initial_conditions: Array of n arrays to instantiate particles. Each array must contain the information
                                   required for the particle constructor: [initial_pos, initial_vel, mass, fixed: bool]
        :param sim_param: Dictionary of other parameters required (k, l0, dt, ...)
        """

        ##TODO: don't need to define these in here, give this an update
        ## can also just call them within the functions
        # e.g. self__k = sim_param["k"] is only used in the __one_d_force_vector function
        # so we can just call it there and leave it out of the beginning
        # THIS does create import problems, just calling self is easier for now

        ##TODO: added copies of arrays, to prevent changes in essence
        # creating an 'immutability illusion', be careful though output must also be copied
        # one can still change variables using an update function!

        # not all variables are made "immutable" because they should actually change throughout
        # the simulation, e.g. the position and velocity of the particles

        # particles
        # self.__position
        # self.__velocity
        # self.__mass
        # self.__fixed
        ##TODO: add a particle specific damping coefficient

        # springdampers, np.copy always returns an np.arrays
        self.__connectivity_matrix = np.copy(connectivity_matrix)
        self.__k = np.copy(sim_param["k"])
        self.__l0 = np.copy(sim_param["l0"])
        self.__is_compression = np.copy(sim_param["is_compression"])
        self.__is_tension = np.copy(sim_param["is_tension"])

        # special for pulleys
        if "is_pulley" in sim_param:
            self.__is_pulley = np.copy(sim_param["is_pulley"])
        # if this entry is not provided, create an empty array where each index is set to False
        else:
            self.__is_pulley = np.zeros((len(self.__connectivity_matrix),), dtype=bool)
        if "pulley_other_line_pair" in sim_param:
            # copy.deepcopy is used, because np.copy only works for np.arrays -this must remain an dict
            self.__pulley_other_line_pair = copy.deepcopy(
                sim_param["pulley_other_line_pair"]
            )

        # special for bending-resistance
        self.__is_rotational = np.copy(sim_param["is_rotational"])
        if "bending_params" in sim_param:
            self.__bending_params = np.copy(sim_param["bending_params"])

        # other
        self.__dt = copy.deepcopy(sim_param["dt"])
        self.__n = copy.deepcopy(sim_param["n"])

        ##TODO: split damping particle & spring-damping specific
        self.__c = copy.deepcopy(sim_param["c"])

        self.__rtol = copy.deepcopy(sim_param["rel_tol"])
        self.__atol = copy.deepcopy(sim_param["abs_tol"])
        self.__maxiter = copy.deepcopy(sim_param["max_iter"])

        # allocate memory
        self.__particles = []
        self.__springdampers = []
        self.__f = np.zeros((self.__n * 3,))
        self.__jx = np.zeros((self.__n * 3, self.__n * 3))
        self.__jv = np.zeros((self.__n * 3, self.__n * 3))

        self.__instantiate_particles(initial_conditions)
        self.__m_matrix = self.__construct_m_matrix()
        self.__instantiate_springdampers()

        # Variables required for kinetic damping & residual force damping
        self.__w_kin = self.__calc_kin_energy()
        self.__w_kin_min1 = self.__calc_kin_energy()
        self.__w_kin_min2 = self.__calc_kin_energy()
        self.__vis_damp = copy.deepcopy(sim_param["is_with_visc_damping"])
        self.__x_min1 = np.zeros(
            self.__n,
        )
        self.__x_min2 = np.zeros(
            self.__n,
        )
        return

    def __str__(self):
        print(
            "ParticleSystem object instantiated with attributes\nConnectivity matrix:"
        )
        print(self.__connectivity_matrix)
        print("Instantiated particles:")
        n = 1
        for particle in self.__particles:
            print(f" p{n}: ", particle)
            n += 1
        return ""

    def __instantiate_particles(self, initial_conditions):
        for set_of_initial_cond in initial_conditions:
            x = set_of_initial_cond[0]
            v = set_of_initial_cond[1]
            m = set_of_initial_cond[2]
            f = set_of_initial_cond[3]
            self.__particles.append(Particle(x, v, m, f))
        return

    def __instantiate_springdampers(self):
        for i, index in enumerate(self.__connectivity_matrix):
            # TODO: check if pulley
            self.__springdampers.append(
                SpringDamper(
                    self.__particles[index[0]],
                    self.__particles[index[1]],
                    self.__k[i],
                    self.__l0[i],
                    self.__c,
                    self.__dt,
                    self.__is_compression[i],
                    self.__is_tension[i],
                    self.__is_pulley[i],
                    self.__is_rotational[i],
                )
            )

        return

    def __construct_m_matrix(self):
        matrix = np.zeros((self.__n * 3, self.__n * 3))

        for i in range(self.__n):
            matrix[i * 3 : i * 3 + 3, i * 3 : i * 3 + 3] += (
                np.identity(3) * self.__particles[i].m
            )

        return matrix

    def __calc_kin_energy(self):
        v = self.__pack_v_current()
        w_kin = np.matmul(
            np.matmul(v.T, self.__m_matrix), v
        )  # Kinetic energy, 0.5 constant can be neglected
        return w_kin

    def simulate(self, f_external: npt.ArrayLike = ()):
        if not len(f_external):
            f_external = np.zeros(
                self.__n * 3,
            )
        f = self.__one_d_force_vector() + f_external

        v_current = self.__pack_v_current()
        x_current = self.__pack_x_current()

        jx, jv = self.__system_jacobians()

        ## old-code
        # # constructing A matrix and b vector for solver
        # A = self.__m_matrix - self.__dt * jv - self.__dt**2 * jx
        # b = self.__dt * f + self.__dt**2 * np.matmul(jx, v_current)

        # Constructing A matrix and b vector for solver
        # TODO: check papers to determine the right-values here
        is_with_rayleigh_damping = False
        is_with_bulk_viscosity = False
        alpha = 0  # (mass, low-freq) Rayleigh damping coefficient, adjust as needed
        beta = 0.001  # (material, high-freq) between 0.06 and 0.1 from Flores2013 - Parachutes) # Rayleigh damping coefficient, adjust as needed
        gamma = 0.001  # Bulk viscosity damping coefficient, adjust as needed

        A = self.__m_matrix - self.__dt * jv - self.__dt**2 * jx

        if is_with_rayleigh_damping:
            A += alpha * self.__m_matrix + beta * self.__dt * jx

        if is_with_bulk_viscosity:
            A += gamma * self.__m_matrix
            b = (
                self.__dt * f
                + self.__dt**2 * np.matmul(jx, v_current)
                - gamma * self.__dt * v_current
            )
        else:
            b = self.__dt * f + self.__dt**2 * np.matmul(jx, v_current)

        # checking conditioning of A and b
        # print("conditioning A:", np.linalg.cond(A))

        for i in range(self.__n):
            if self.__particles[i].fixed:
                A[i * 3 : (i + 1) * 3] = 0  # zeroes out row i to i + 3
                A[:, i * 3 : (i + 1) * 3] = 0  # zeroes out column i to i + 3
                b[i * 3 : (i + 1) * 3] = 0  # zeroes out row i

        # BiCGSTAB from scipy library
        dv, _ = bicgstab(
            A, b, rtol=self.__rtol, atol=self.__atol, maxiter=self.__maxiter
        )

        v_next = v_current + dv

        ##TODO: this is where one would apply damping?
        v_next = self.__c * v_next

        ##TODO: this shouls be coming through the config object class
        is_with_vel_cap = False
        psm_vel_cap = 5
        if is_with_vel_cap and np.max(v_next) > psm_vel_cap:
            v_next = v_next / np.max(v_next)

        x_next = x_current + self.__dt * v_next

        # function returns the pos. and vel. for the next timestep, but for fixed particles this value doesn't update!
        self.__update_x_v(x_next, v_next)
        return x_next, v_next

    def kin_damp_sim(
        self, f_ext: npt.ArrayLike, q_correction: bool = False
    ):  # kinetic damping alghorithm
        # TODO: this was not commented in Alex's versions
        # if self.__vis_damp:  # Condition resetting viscous damping to 0
        #     self.__c = 0
        #     self.__springdampers = []
        #     self.__instantiate_springdampers()
        #     self.__vis_damp = False

        if len(f_ext):  # condition checking if an f_ext is passed as argument
            self.__save_state()
            x_next, v_next = self.simulate(f_ext)
        else:
            self.__save_state()
            x_next, v_next = self.simulate()

        w_kin_new = self.__calc_kin_energy()

        if (
            w_kin_new > self.__w_kin
        ):  # kin damping algorithm, takes effect when decrease in kin energy is detected
            self.__update_w_kin(w_kin_new)
        else:
            v_next = np.zeros(
                self.__n * 3,
            )

            if (
                q_correction
            ):  # statement to check if q_correction is desired, standard is turned off
                q = (self.__w_kin - w_kin_new) / (
                    2 * self.__w_kin - self.__w_kin_min1 - w_kin_new
                )
                # print(q)
                # print(self.__w_kin, w_kin_new)
                # !!! Not sure if linear interpolation between states is the way to determine new x_next !!!
                if q < 0.5:
                    x_next = self.__x_min2 + (q / 0.5) * (self.__x_min1 - self.__x_min2)
                elif q == 0.5:
                    x_next = self.__x_min1
                elif q < 1:
                    x_next = self.__x_min1 + ((q - 0.5) / 0.5) * (
                        x_next - self.__x_min1
                    )

                # Can also use this q factor to recalculate the state for certain timestep h

            self.__update_x_v(x_next, v_next)
            self.__update_w_kin(0)

        return x_next, v_next

    def __pack_v_current(self):
        return np.array([particle.v for particle in self.__particles]).flatten()

    def __pack_x_current(self):
        return np.array([particle.x for particle in self.__particles]).flatten()

    def __one_d_force_vector(self):
        self.__f[self.__f != 0] = 0

        # calling this once, instead of for each springdamper
        x_current = self.__pack_x_current()
        # print(f"len(x_current): {len(x_current)} ")

        for idx in range(
            len(self.__springdampers)
        ):  # TODO: for idx,springdamper in enumarata(self.__springdampers):
            if self.__is_pulley[idx]:  # if pulley #TODO: springdamper.is_pulley
                idx_p3, idx_p4, rest_length_p3p4 = self.__pulley_other_line_pair[
                    str(idx)
                ]
                # points are flat, so we need to configure a bit
                p3 = np.array([x_current[int(idx_p3) * 3 : int(idx_p3) * 3 + 3]])
                p4 = np.array([x_current[int(idx_p4) * 3 : int(idx_p4) * 3 + 3]])
                norm_p3p4 = np.linalg.norm(p3 - p4)
                delta_length_pulley_other_line = norm_p3p4 - rest_length_p3p4
                fs, fd = self.__springdampers[idx].force_value(
                    delta_length_pulley_other_line
                )

                i, j = self.__connectivity_matrix[idx]
                self.__f[i * 3 : i * 3 + 3] += fs + fd
                self.__f[j * 3 : j * 3 + 3] -= fs + fd

            # TODO: put rotational back
            # elif self.__is_rotational[
            #     idx
            # ]:  # elif springdamper.is_rotational # if in need for a linear bending spring
            #     # print(f" ")
            #     # print(f"-----------------------------------------------")
            #     # print(f"spring idx: {idx}")
            #     # getting the additional info
            #     alpha_a, alpha_b, alpha_e, idx_pe = self.__bending_params[f"{idx}"]
            #     # print(f"alpha_a: {alpha_a}")
            #     # print(f"alpha_b: {alpha_b}")
            #     # print(f"alpha_e: {alpha_e}")
            #     # print(f"idx_pe: {idx_pe}")

            #     # getting the connectivity, idx_pa and idx_pb are node 1 and 2 and pe is the point in between
            #     idx_pa, idx_pb = self.__connectivity_matrix[idx]
            #     # print(f"self.__connectivity_matrix]: {self.__connectivity_matrix}")
            #     # print(
            #     #     f"self.__connectivity_matrix[{idx}]: {self.__connectivity_matrix[idx]}"
            #     # )
            #     # print(f"idx_pa: {idx_pa}")
            #     # print(f"idx_pb: {idx_pb}")
            #     pa = np.array(x_current[int(idx_pa) * 3 : int(idx_pa) * 3 + 3])
            #     pb = np.array(x_current[int(idx_pb) * 3 : int(idx_pb) * 3 + 3])
            #     pe = np.array(x_current[int(idx_pe) * 3 : int(idx_pe) * 3 + 3])
            #     # print(f"Pa: {pa}")
            #     # print(f"pb: {pb}")
            #     # print(f"pe: {pe}")

            #     # Calculating the projection vector of e onto the line a-b
            #     pe_proj = pa + np.dot(pe - pa, pb - pa) / np.dot(pb - pa, pb - pa) * (
            #         pb - pa
            #     )

            #     # print(f"pe_proj: {pe_proj}")

            #     # calculating the bending vector bending_vector, which should make the line straight again
            #     bending_vector = pe_proj - pe
            #     # calculating the forces on each point
            #     # where alpha is dict, that contains the relative contributions of each point
            #     # alpha_a and alpha_b are positive and alpha_e is negative
            #     fa = -self.__k[idx] * alpha_a * bending_vector
            #     fb = -self.__k[idx] * alpha_b * bending_vector
            #     fc = -self.__k[idx] * alpha_e * bending_vector

            #     # print(f"bending_vector: {bending_vector}")
            #     # print(f"fa: {fa}")
            #     # print(f"fb: {fb}")
            #     # print(f"fc: {fc}")

            #     # # addding the forces to the 1d vector
            #     self.__f[idx_pa * 3 : idx_pa * 3 + 3] += fa
            #     self.__f[idx_pb * 3 : idx_pb * 3 + 3] += fb
            #     self.__f[idx_pe * 3 : idx_pe * 3 + 3] += fc

            else:  # if just spring-damper
                fs, fd = self.__springdampers[idx].force_value()
                i, j = self.__connectivity_matrix[idx]
                self.__f[i * 3 : i * 3 + 3] += fs + fd
                self.__f[j * 3 : j * 3 + 3] -= fs + fd

        return self.__f

    def __system_jacobians(self):
        self.__jx[self.__jx != 0] = 0
        self.__jv[self.__jv != 0] = 0

        for n in range(len(self.__springdampers)):
            jx, jv = self.__springdampers[n].calculate_jacobian()
            i, j = self.__connectivity_matrix[n]

            self.__jx[i * 3 : i * 3 + 3, i * 3 : i * 3 + 3] += jx
            self.__jx[j * 3 : j * 3 + 3, j * 3 : j * 3 + 3] += jx
            self.__jx[i * 3 : i * 3 + 3, j * 3 : j * 3 + 3] -= jx
            self.__jx[j * 3 : j * 3 + 3, i * 3 : i * 3 + 3] -= jx

            self.__jv[i * 3 : i * 3 + 3, i * 3 : i * 3 + 3] += jv
            self.__jv[j * 3 : j * 3 + 3, j * 3 : j * 3 + 3] += jv
            self.__jv[i * 3 : i * 3 + 3, j * 3 : j * 3 + 3] -= jv
            self.__jv[j * 3 : j * 3 + 3, i * 3 : i * 3 + 3] -= jv

        # i = 0
        # j = 1
        # for springdamper in self.__springdampers:
        #     jx, jv = springdamper.calculate_jacobian()
        #     if j > 40:
        #         print(i, j)
        #         print(j * 3, j * 3 + 3, j * 3, j * 3 + 3)
        #         print()
        #     self.__jx[i * 3:i * 3 + 3, i * 3:i * 3 + 3] += jx
        #     self.__jx[j * 3:j * 3 + 3, j * 3:j * 3 + 3] += jx
        #     self.__jx[i * 3:i * 3 + 3, j * 3:j * 3 + 3] -= jx
        #     self.__jx[j * 3:j * 3 + 3, i * 3:i * 3 + 3] -= jx
        #
        #     self.__jv[i * 3:i * 3 + 3, i * 3:i * 3 + 3] += jv
        #     self.__jv[j * 3:j * 3 + 3, j * 3:j * 3 + 3] += jv
        #     self.__jv[i * 3:i * 3 + 3, j * 3:j * 3 + 3] -= jv
        #     self.__jv[j * 3:j * 3 + 3, i * 3:i * 3 + 3] -= jv
        #
        #     i += 1
        #     j += 1

        return self.__jx, self.__jv

    def __update_x_v(self, x_next: npt.ArrayLike, v_next: npt.ArrayLike):
        for i in range(self.__n):
            self.__particles[i].update_pos(x_next[i * 3 : i * 3 + 3])
            self.__particles[i].update_vel(v_next[i * 3 : i * 3 + 3])
        return

    def __update_w_kin(self, w_kin_new: float):
        self.__w_kin_min2 = self.__w_kin_min1
        self.__w_kin_min1 = self.__w_kin
        self.__w_kin = w_kin_new
        return

    # def __update_f_res(self, f_res_new: float):
    #     self.__f_res_min2 = self.__f_res_min1
    #     self.__f_res_min1 = self.__f_res
    #     self.__f_res = f_res_new
    #     return

    def __save_state(self):
        self.__x_min2 = self.__x_min1
        self.__x_min1 = self.__pack_x_current()
        return

    ##TODO: note to self, below is a regular method, rather than a calleble property
    # and thus we leave it without the @property
    # function to update the rest length of a springdamper
    # IMPORTANT, only the rest length of the springdamper is updated, not the springdamper itself!!
    def update_rest_length(self, idx: int, delta_rest_length: float):
        self.__springdampers[idx].update_rest_length(delta_rest_length)
        return

    def update_damping(self, new_damping_ratio: float):
        self.__c = new_damping_ratio
        return

    @property
    def extract_rest_length(self):
        return np.array(
            [springdamper.return_rest_length for springdamper in self.__springdampers]
        )

    @property
    def f_int(self):
        return self.__f

    ##TODO: not sure if this is the most-efficient
    @property
    def x_current_2D(self):
        return np.array([particle.x for particle in self.__particles])

    @property
    def particles(self):
        return self.__particles

    ##TODO: these are commented out for now, as they are not used.
    # @property
    # def particles(self):  # Temporary solution to calculate external aerodynamic forces
    #     return self.__particles

    # @property
    # def springdampers(self):
    #     return self.__springdampers

    # @property
    # def stiffness_m(self):
    #     self.__system_jacobians()
    #     return self.__jx

    @property
    def x_v_current(self):
        return self.__pack_x_current(), self.__pack_v_current()


if __name__ == "__main__":
    c_matrix = [[0, 1], [1, 0]]
    init_cond = [[[0, 0, 0], [0, 0, 0], 1, True], [[0, 0, 0], [0, 0, 0], 1, False]]

    params = {
        # model parameters
        "n": 2,  # [-] number of particles
        "k": 2e4,  # [N/m] spring stiffness
        "c": 0,  # [N s/m] damping coefficient
        "l0": 0,  # [m] rest length
        # simulation settings
        "dt": 0.001,  # [s] simulation timestep
        "t_steps": 1000,  # [-] number of simulated time steps
        "abs_tol": 1e-50,  # [m/s] absolute error tolerance iterative solver
        "rel_tol": 1e-5,  # [-] relative error tolerance iterative solver
        "max_iter": 1e5,  # [-] maximum number of iterations
        # physical parameters
        "g": 9.81,  # [m/s^2] gravitational acceleration
    }

    ps = ParticleSystem(c_matrix, init_cond, params)
    print(ps)
    print(ps.system_energy)
    pass
