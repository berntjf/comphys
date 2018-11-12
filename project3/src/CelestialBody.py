import numpy as np
from src.UnitSet import UnitSet

unit = UnitSet({'AU' : 1, 'day' : 1, 'kg' : 1})

class CelestialBody:
    def __init__(self, name, pos, vel, mass):
        self.name = name
        self.pos = pos
        self.vel = vel
        self.mass = mass

    def kinetic_energy(self):
        v = self.vel
        return 0.5 * self.mass * (v[0]**2 + v[1]**2 + v[2]**2)

    def potential_energy(self, other, G):
        r = self.pos - other.pos
        return - G * self.mass * other.mass / (r[0]**2 + r[1]**2 + r[2]**2)**0.5

if __name__=='__main__':
    CelestialBody(np.array([1., 1., 1.]), np.array([1., 1., 1.]), 6e24)
