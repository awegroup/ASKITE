"""
Child Abstract Base Class 'ImplicitForce', for implicit force objects to be instantiated in ParticleSystem
"""
from src.particleSystem.Force import Force
from src.particleSystem.Particle import Particle
from abc import abstractmethod
from abc import abstractproperty


class ImplicitForce(Force):
    def __init__(self, p1: Particle, p2: Particle):
        self.__p1 = p1
        self.__p2 = p2
        super().__init__()
        # TODO: when instantiating a base class, you also need to call
        # the base class's __init__ method. This is done by calling super().__init__().
        return

    def __str__(self):
        return

    @abstractmethod
    def calculate_jacobian(self):
        return

    @property
    def p1(self):
        return self.__p1

    @property
    def p2(self):
        return self.__p2


if __name__ == "__main__":
    pass
