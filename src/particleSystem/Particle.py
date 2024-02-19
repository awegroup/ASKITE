"""
Child Class 'Particle', for particle objects to be instantiated in ParticleSystem
"""
from src.particleSystem.SystemObject import SystemObject
import numpy as np
import numpy.typing as npt


class Particle(SystemObject):

    def __init__(self, x: npt.ArrayLike, v: npt.ArrayLike, m: float, fixed: bool):
        self.__x = np.array(x)
        self.__v = np.array(v)
        self.__m = m
        self.__fixed = fixed
        super().__init__()
        return

    def __str__(self):
        return f"Particle Object, position [m]: [{self.__x[0]}, {self.__x[1]}, {self.__x[2]}], " \
               f"velocity [m/s]: [{self.__v[0]}, {self.__v[1]}, {self.__v[2]}], mass [kg]: {self.__m}" \
               f", fixed: {self.__fixed}"

    def update_pos(self, new_pos: npt.ArrayLike):
        if not self.__fixed:
            self.__x = np.array(new_pos)
        return

    def update_vel(self, new_vel: npt.ArrayLike):
        if not self.__fixed:
            self.__v = np.array(new_vel)
        return

    @property
    def x(self):
        return self.__x

    @property
    def v(self):
        return self.__v

    @property
    def m(self):
        return self.__m

    @property
    def fixed(self):
        return self.__fixed


if __name__ == "__main__":
    position = [0, 0, 0]
    velocity = [0, 0, 0]
    mass = 1
    fixed1 = False
    fixed2 = True
    p1 = Particle(position, velocity, mass, fixed1)
    p2 = Particle(position, velocity, mass, fixed2)
    print(p1)
    print(p2)
    updated_pos = [0, 0, 1]
    updated_vel = [0, 0, 1]
    p1.update_pos(updated_pos)
    p1.update_vel(updated_vel)
    p2.update_pos(updated_pos)
    p2.update_vel(updated_vel)
    print(p1)
    print(p2)
    pass
