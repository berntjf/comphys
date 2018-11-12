import sys
import numpy as np
import matplotlib.pyplot as plt

from numba import jit
from numpy import newaxis
from scipy.linalg import norm

from src.CelestialBody import CelestialBody
from src.integrate import *

def read_planets(file_name, unit_set):
    objects = []
    with open(file_name, 'r') as infile:
        iter = 0
        prev_time_stamp = None
        last_line = infile.readline()
        while last_line != '':
            name = last_line.split()[1]
            mass = infile.readline().split()[1]
            current_time_stamp = infile.readline()
            line = infile.readline()
            pos = np.array([float(line[4:26]), float(line[30:52]), float(line[56:78])]) * unit_set.AU
            line = infile.readline()
            vel = np.array([float(line[4:26]), float(line[30:52]), float(line[56:78])]) * unit_set.AU / unit_set.day
            if current_time_stamp != prev_time_stamp and not prev_time_stamp is None:
                print("Looks like sample is taken at different times! That does not work.")
                sys.exit(1)
            objects.append(CelestialBody(name, pos, vel, mass))
            last_line = infile.readline()
    return np.array(objects)

class SolarSystem:
    def __init__(self, unit, planets_no=None):
        self.objects = read_planets('data/planets.dat', unit)[planets_no]
        self.unit = unit
        self.x0 = np.array([object.pos for object in self.objects]).T
        self.v0 = np.array([object.vel for object in self.objects]).T
        self.no_objects = len(self.objects)

    def get_neg_G_masses(self):
        masses = np.empty(self.no_objects)
        for i in range(self.no_objects):
            masses[i] = self.objects[i].mass
        return - self.unit.G * masses

    @staticmethod
    @jit
    def a_g(neg_G_mass, r):
        r_square = r[0]**2
        for i in range(1, len(r)):
            r_square += r[i]**2
        return neg_G_mass * r_square**(-1.5) * r

    def calculate_orbits(self, solver, t, steps=1):
        neg_G_masses = self.get_neg_G_masses()
        no_objects = self.no_objects
        a_g = self.a_g
        @jit
        def f(positions):
            result = np.zeros((3, no_objects))
            for j in range(no_objects):
                for i in range(no_objects):
                    if i != j:
                        result[:,j] += a_g(neg_G_masses[i], positions[:,j]-positions[:,i])
            return result
        self.xi, self.vi, self.ti = solver(f, self.x0, self.v0, t, steps=steps)

    def calculate_orbits_simple(self, solver, t, steps=1):
        """Assuming the Sun is standing still at the center of mass at pos r_vec = 0."""
        neg_G_masses = self.get_neg_G_masses()
        no_objects = self.no_objects
        a_g = self.a_g
        @jit
        def f(positions):
            result = np.zeros((3, no_objects))
            for j in range(1, no_objects):
                result[:,j] += a_g(neg_G_masses[0], positions[:,j])
                for i in range(1, no_objects):
                    if i != j:
                        result[:,j] += a_g(neg_G_masses[i], positions[:,j]-positions[:,i])
            return result
        self.xi, self.vi, self.ti = solver(f, self.x0, self.v0, t, steps=steps)

    def calculate_orbits_one_body(self, solver, t, steps=1):
        neg_G_masses = self.get_neg_G_masses()
        no_objects = self.no_objects
        a_g = self.a_g
        @jit
        def f(positions):
            return a_g(neg_G_masses[0], positions)
        xi_planet, vi_planet, self.ti = solver(f, self.x0[:,1], self.v0[:,1], t, steps=steps)
        self.xi = np.zeros(xi_planet.shape+(2,))
        self.xi[:,:,1] = xi_planet
        self.vi = np.zeros(vi_planet.shape+(2,))
        self.vi[:,:,1] = vi_planet

    def calculate_orbits_relativistic(self, t, steps=1):
        solver = ForwardEulerVelocity()
        neg_G_masses = self.get_neg_G_masses()
        no_objects = self.no_objects
        a_g = self.a_g
        c = self.unit.c
        #@jit
        def f(r, v):
            result = np.zeros((3, no_objects))
            for j in range(no_objects):
                for i in range(no_objects):
                    if i != j:
                        result[:,j] += a_g(neg_G_masses[i], r[:,j]-r[:,i])
            l_sq = (r[1]*v[2])**2 + (r[2]*v[0])**2 + (r[0]*v[1])**2
            result = result * (1 + 3 * l_sq / (r[0]**2+r[1]**2+r[2]**2) / c**2)
            return result
        self.xi, self.vi, self.ti = solver(f, self.x0, self.v0, t, steps=steps)

    def show_animation(self, axis=10, title=""):
        plt.ion()
        for i in range(self.ti.size):
            pos = self.xi[i] / self.unit.AU
            for index, pos_planet in enumerate(pos.T):
                plt.scatter(pos_planet[0], pos_planet[1])
                plt.annotate(self.objects[index], (pos_planet[0], pos_planet[1]))
            plt.title(title)
            plt.xlabel('x pos [AU]')
            plt.ylabel('y pos [AU]')
            plt.axis([-axis, axis, -axis, axis])
            plt.pause(0.001)
            plt.clf()
        plt.ioff()

    def plot_by_time(self, y_array, label=""):
        plt.plot(self.ti / self.unit.day, y_array, label=label)
        plt.xlabel('Time [day]')
        if label != "":
            plt.legend()

    def show_parameter_plot(self, axis=10, title=""):
        for i in range(self.no_objects):
            plt.plot(self.xi[:,0,i], self.xi[:,1,i])
        plt.scatter(0, 0)
        plt.title(title)
        plt.axis([-axis, axis, -axis, axis])
        plt.show()
