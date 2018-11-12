import numpy as np
from src.UnitSet import UnitSet
from src.SolarSystem import SolarSystem
from src.integrate import *

unit = UnitSet({'AU' : 1, 'yr' : 1, 'kg' : 1})

ss = SolarSystem(unit, planets_no=[0, 3])
t = np.linspace(0, unit.yr, 52)
#ss.calculate_orbits(Verlet(), ic, t, steps=1000)
#ss.show_animation(4, title='Verlet integration')
#ss.show_animation(35, title='Verlet integration')
ss.calculate_orbits_one_body(ForwardEuler(), t, steps=500)
ss.show_animation(title='Forward Euler')
#ss.show_animation(35, title='Forward Euler')
#print(ss.x)
