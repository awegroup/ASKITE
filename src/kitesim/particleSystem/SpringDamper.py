"""
Child Class 'SpringDamper', for spring-damper objects to be instantiated in ParticleSystem
"""

import numpy as np

from kitesim.particleSystem.ImplicitForce import ImplicitForce
from kitesim.particleSystem.Particle import Particle


class SpringDamper(ImplicitForce):
    def __init__(
        self,
        p1: Particle,
        p2: Particle,
        k: float,
        l0: float,
        c: float,
        dt: float,
        is_compression: bool = True,
        is_tension: bool = True,
        is_pulley: bool = False,
        is_rotational: bool = False,
    ):
        self.__k = k
        self.__c = c
        self.__l0 = l0
        self.__is_compression = is_compression
        self.__is_tension = is_tension
        self.__is_pulley = is_pulley
        self.__is_rotational = is_rotational
        self.__dt = dt
        super().__init__(p1, p2)
        return

    def __str__(self):
        return (
            f"SpringDamper object, spring stiffness [n/m]: {self.__k}, rest length [m]: {self.__l0}\n"
            f"Damping coefficient [N s/m]: {self.__c}"
            f"Assigned particles\n  p1: {self.p1}\n  p2: {self.p2}"
        )

    def __relative_pos(self):
        return np.array([self.p1.x - self.p2.x])

    def __relative_vel(self):
        return np.array([self.p1.v - self.p2.v])

    def force_value(self, extra_rest_length: float = 0):
        # TODO: HERE the check, what kind of "edge" this is.
        # IF pulley
        # ELIF rotational
        # ElIF compression
        # ELIF tension
        # ELSE error!
        return (
            self.__calculate_f_spring(extra_rest_length),
            self.__calculate_f_damping(),
        )

    # To allow actuation mid-loop
    def update_rest_length(self, delta_rest_length: float = 0):
        self.__l0 += delta_rest_length
        return

    # def calculate_pulley:

    # TODO: write multiple version for each "edge" type
    def __calculate_f_spring(self, delta_length_pulley_other_line: float = 0):
        relative_pos = self.__relative_pos()
        norm_pos = np.linalg.norm(relative_pos)

        rest_length = self.__l0
        delta_length = norm_pos - rest_length

        if self.__is_pulley:  # checking for pulleys
            delta_length += delta_length_pulley_other_line

        # if extended AND no tension-resistance
        # or if compressed AND no compression-resistance
        # or if delta_length is 0 OR distance is zero
        if delta_length > 0 and not self.__is_tension:
            unit_vector = np.array([0, 0, 0])
        elif delta_length < 0 and not self.__is_compression:
            unit_vector = np.array([0, 0, 0])
        elif delta_length == 0 or norm_pos == 0:
            unit_vector = np.array([0, 0, 0])
        else:  # if compressed / stretch with resistance
            unit_vector = relative_pos / norm_pos

        f_spring = -self.__k * delta_length * unit_vector

        # if self.__is_rotational:
        #     print("rotational_element detected")

        return np.squeeze(f_spring)

    def __calculate_f_damping(self):
        relative_pos = self.__relative_pos()
        relative_vel = np.squeeze(self.__relative_vel())
        norm_pos = np.linalg.norm(relative_pos)

        if norm_pos != 0:
            unit_vector = np.squeeze(relative_pos / norm_pos)
        else:
            unit_vector = np.squeeze(np.array([0, 0, 0]))

        f_damping = -self.__c * np.dot(relative_vel, unit_vector) * unit_vector
        return np.squeeze(f_damping)

    def calculate_jacobian(self):
        # TODO: current implementation assumes linear spring-damper
        relative_pos = self.__relative_pos()
        norm_pos = np.linalg.norm(relative_pos)

        if norm_pos != 0:
            unit_vector = relative_pos / norm_pos
        else:
            norm_pos = 1
            unit_vector = np.array([0, 0, 0])

        i = np.identity(3)
        T = np.matmul(np.transpose(unit_vector), unit_vector)
        jx = -self.__k * ((self.__l0 / norm_pos - 1) * (T - i) + T)

        jv = -self.__c * i

        return jx, jv

    @property  # TODO: now you can call springdamper.is_pulley! :)
    def is_pulley():
        return self.__is_pulley

    @property
    def is_rotational():
        return self.__is_rotational

    @property
    def return_rest_length(self):
        return self.__l0


if __name__ == "__main__":
    particle1 = Particle([0, 0, 0], [0, 0, 0], 1, False)
    particle2 = Particle([0, 0, 1], [0, 0, 1], 1, False)
    stiffness = 1e5
    damping = 10
    rest_length = 0
    timestep = 0.001

    springdamper = SpringDamper(
        particle1, particle2, stiffness, rest_length, damping, timestep
    )

    print(springdamper)
    print()
    print(springdamper.force_value())
    # print((np.sqrt(3)-1)*1e5/np.sqrt(3))  # value check
    print()
    print(springdamper.calculate_jacobian())
    pass
