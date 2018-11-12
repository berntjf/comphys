import numpy as np
import matplotlib.pyplot as plt

from scipy.linalg import norm

from src.integrate import *
from src.read import read_planets, read_ic
from src.SolarSystem import SolarSystem
from src.UnitSet import UnitSet

unit = UnitSet({'AU' : 1, 'day' : 1, 'kg' : 1})
objects = read_planets('data/planets.dat', unit)
x0, v0 = read_ic('data/ic.dat', unit)

v_unit_vec = np.array([-x0[1], x0[0], np.zeros(objects.size)]) / (x0[0]**2 + x0[1]**2)**0.5
v0 =  (self.G * objects[0].mass / norm(r, axis=0))**0.5 * v_unit_vec
v0[:,0] = np.zeros(3)

solar_system = SolarSystem(objects=objects, x0=x0, v0=v0, unit=unit)

for i in range(1, 8):
    N = 3**i
    t = np.linspace(0, unit.yr, N)
    plt.title("%g time steps" % N)
    solar_system.calculate_orbits(Verlet(), t)
    solar_system.plot(solar_system.xi[:,0,3], label="x")
    solar_system.plot(solar_system.xi[:,1,3], label="y")
    solar_system.plot(solar_system.xi[:,2,3], label="z")
    plt.show()
