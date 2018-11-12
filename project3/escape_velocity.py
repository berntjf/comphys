import numpy as np
import matplotlib.pyplot as plt

from scipy.linalg import norm

from src.integrate import *
from src.read import read_planets, read_ic
from src.SolarSystem import a_g
from src.UnitSet import UnitSet

unit = UnitSet({'AU' : 1, 'day' : 1, 'kg' : 1})
objects = read_planets('data/planets.dat', unit)
solver = VelocityVerlet()

sun = objects[0]
earth = objects[3]

neg_G_mass = - unit.G * sun.mass

def a_g(neg_G_mass, r, beta):
    r_square = r[0]**2
    for x in r[1:]:
        r_square += x**beta
    return neg_G_mass * r_square**(-1.5) * r

v0 = np.array([0, np.sqrt(2 * unit.G * sun.mass / unit.AU)])


for beta in [2, 2.05, 2.2, 2.5, 3]:
    def f(r):
        return a_g(neg_G_mass, r, beta)

    solver(f, x0=np.array([unit.AU, 0]), v0=v0, t=np.linspace(0, 10*unit.yr, 30), steps=1000)

    axis = 25
    plt.ion()
    for i in range(solver.ti.size):
        pos = solver.xi[i] / unit.AU
        plt.scatter(0, 0, color='yellow')
        plt.scatter(pos[0], pos[1])
        plt.title('Escape velocity, F=GMm/r^%g, %3.1f years' % (beta, (solver.ti[i] / unit.yr)))
        plt.xlabel('x pos [AU]')
        plt.ylabel('y pos [AU]')
        plt.axis([-axis, axis, -axis, axis])
        plt.pause(0.05)
        plt.clf()
    plt.ioff()
