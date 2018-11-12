import numpy as np

from scipy.linalg import norm

from src.integrate import *
from src.read import read_planets, read_ic
from src.SolarSystem import SolarSystem
from src.UnitSet import UnitSet

unit = UnitSet({'AU' : 1, 'day' : 1, 'kg' : 1})
objects = read_planets('data/planets.dat', unit)
x0, v0 = read_ic('data/ic.dat', unit)

v0 = SolarSystem.velocity_for_circular_orbit(unit.G, objects[0].mass, x0)
v0[:,0] = np.zeros(3)

solar_system = SolarSystem(objects=objects, x0=x0, v0=v0, unit=unit)

t = np.linspace(0, unit.yr, 37)
solar_system.calculate_orbits(VelocityVerlet(), t, steps=100)
import matplotlib.pyplot as plt
pot_energy = solar_system.pot_energy()
kin_energy = np.sum(0.5 * solar_system.masses * solar_system.vi**2, axis=(1,2))


plt.title('Potential energy')
solar_system.plot(pot_energy / unit.Joule)
plt.savefig('plot/3c_potential_energy.png')
plt.clf()

plt.title('Kinetic energy')
solar_system.plot(kin_energy / unit.Joule)
plt.savefig('plot/3c_kinetic_energy.png')
plt.clf()

plt.title('Angular momentum')
solar_system.plot(np.sum((norm(solar_system.xi, axis=1) * solar_system.masses * norm(solar_system.vi, axis=1) / (unit.kg * unit.meter**2 / unit.sec))[:,1:], axis=1))
plt.savefig('plot/3c_angular_momentum.png')
