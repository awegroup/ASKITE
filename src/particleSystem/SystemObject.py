"""
Parent Abstract Base Class 'SystemObject', for objects to be instantiated in ParticleSystem
"""
from abc import ABC

# TODO: does not contain particlesytem, to avoid circular import.


class SystemObject(ABC):
    def __init__(self):
        return

    def __str__(self):
        return


if __name__ == "__main__":
    pass
