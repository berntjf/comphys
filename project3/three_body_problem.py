import numpy as np
from src.UnitSet import UnitSet
from src.SolarSystem import SolarSystem
from src.integrate import *

unit = UnitSet({'AU' : 1, 'yr' : 1, 'kg' : 1})

ss = SolarSystem(unit, planets_no=(0, 3, 5))

real_Jupiter_mass = ss.objects[2].mass
# Real mass of Jupiter
ss.objects[2].mass = real_Jupiter_mass
ss.calculate_orbits_simple(solver=VelocityVerlet(), t=np.linspace(0, 10*unit.yr, 1000), steps=10)
ss.show_parameter_plot(axis=6)

# 10 times mass of Jupiter
ss.objects[2].mass = 10 * real_Jupiter_mass
ss.calculate_orbits_simple(solver=VelocityVerlet(), t=np.linspace(0, 10*unit.yr, 1000), steps=10)
ss.show_parameter_plot(axis=6)

# 1000 times mass of Jupiter
ss.objects[2].mass = 1000 * real_Jupiter_mass
ss.calculate_orbits_simple(solver=VelocityVerlet(), t=np.linspace(0, 10*unit.yr, 1000), steps=10)
ss.show_parameter_plot(axis=6)
