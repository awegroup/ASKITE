"""
Child Abstract Base Class 'Force', for force objects to be instantiated in ParticleSystem
"""

from abc import abstractmethod

from .SystemObject import SystemObject


class Force(SystemObject):

    def __init__(self):
        super().__init__()
        return

    def __str__(self):
        return

    @abstractmethod
    def force_value(self):
        return


if __name__ == "__main__":
    pass
