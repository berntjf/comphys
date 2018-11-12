import numpy as np
from src.UnitSet import UnitSet
from src.SolarSystem import SolarSystem
from src.integrate import *

unit = UnitSet({'AU' : 1, 'yr' : 1, 'kg' : 1})

# 3 objects
ss = SolarSystem(unit, planets_no=[0, 1])
ss.calculate_orbits_relativistic(t=np.linspace(0, 10*unit.yr, 1000), steps=10)
ss.show_parameter_plot(axis=0.5)
