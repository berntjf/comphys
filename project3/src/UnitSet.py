import math
import numpy as np
from math import tau
from numpy import exp
from fractions import Fraction

c_SI    = 299792458
G_SI    = Fraction(667428,     10**16)
hbar_SI = Fraction(10545718,   10**41)
kB_SI   = Fraction(138064852,  10**31)
e_SI    = Fraction(1602176565, 10**28)
secs_per_day = Fraction(86400, 1)
days_per_yr = Fraction(36524, 100) # Julian year
AU_SI = 149597870700
Avogadro = 602214129*10**15

class UnitSet:
    """
    c    = speed of light
    gg   = Gauss law of gracity constant
    hbar = reduced Planck's constant
    kB   = Boltzmann constant
    G    = Newton's gracity constant = gg / 4pi
    """

    def __init__(self, given_units):
        self._define_basic_units(given_units)

        meter, sec, kg, Coulomb = self.meter, self.sec, self.kg, self.Coulomb
        c, G, hbar = self.c, self.G, self.hbar


    def __call__(self, unit_name):
        # return self.units[unit_name]
        return getattr(self, unit_name)

    def _define_basic_units(self, given_units):
        # Making a dictionary
        self.units = {}
        for elem in ['meter', 'AU', 'sec', 'day', 'yr', 'kg', 'Kelvin', 'c', 'G', 'hbar', 'kB', 'e', 'Coulomb']:
            self.units[elem] = float('nan')
        for elem in given_units:
            self.units[elem] = given_units[elem]
        for elem in self.units:
            exec('self.%s = float(self.units["%s"])' % (elem, elem))
            exec('global %s; %s = float(self("%s"))' % (elem, elem, elem))

        def set(name, value):
            try:
                before = self(name)
            except AttributeError:
                exec('global %s; %s = float("nan")' % (name, name))
                exec('self.%s = float("nan")' % (name))
                before = self(name)
            if math.isnan(before) and not math.isnan(value):
                self.units[name] = value
                exec('self.%s = float(self.units["%s"])' % (name, name))
                exec('global %s; %s = float(self("%s"))' % (name, name, name))

        # Defining units
        prev_units = None
        while self.units != prev_units: # dvs. at det er framgang
            prev_units = self.units.copy()
            set('meter',    1 / AU_SI           * AU)
            set('meter',    1 / c_SI            * sec * c)
            set('AU',       AU_SI               * meter)
            set('sec',      c_SI                * meter / c)
            set('sec',      1 / secs_per_day    * day)
            set('day',      secs_per_day        * sec)
            set('day',      1 / days_per_yr     * yr)
            set('yr',       days_per_yr         * day)
            set('kg',       1 / hbar_SI         * meter**(-2) * sec * hbar )
            set('Kelvin',                         meter**2 / sec**2 * kg / kB)
            set('c',        c_SI                * meter / sec)
            set('G',        G_SI * meter**3 / sec**2 / kg)
            set('hbar',     hbar_SI * meter**2 / sec * kg)
            set('kB',       kB_SI * meter**2 / sec**2 * kg / Kelvin)
            set('e',        e_SI * Coulomb)
            set('Coulomb',  e / e_SI)

            # Time
            set('minute',   60                  * sec)
            set('hour',      3600.   * sec)

            # Length
            set('ly',      yr * c)

            # Mass
            set('me',      9.1093826e-31              * kg)
            set('mn',      1.6749e-27                 * kg)
            set('Dalton',      1e-3/Avogadro              * kg)
            set('gram',      1e-3                       * kg)
            set('tonne',      1e+3                       * kg)
            set('solar_mass',      tau**2*AU_SI**3*yr**-2  * meter**3 / G)

            # Force, energy, power, etc.
            set('Newton',      1.       * meter * sec**-2 * kg)
            set('Joule',      1.       * kg * meter**2 / sec**2)
            set('Watt',      1.       * kg * meter**2 * sec**-3)
            set('h',      tau      * hbar)

            # Electricricity
            set('mu0', 2*10**7 * tau * meter * kg / Coulomb**2)
            set('Ampere',      1.        * (Coulomb / sec))
            set('Volt',      1.        * (meter**2 / sec**2 * kg / Coulomb))
            set('epsilon0',      1.        /mu0 / c**2)
            set('ke',      0.5   * (mu0 / tau * c**2))
            set('Coulomb', Ampere * sec)
            set('Coulomb', meter**2 / sec**2 * kg / Volt)

            # Dimensionless
            set('deg',      tau/360)
            set('Avogadro',      Avogadro)

        if math.isnan(self.units['Coulomb']):
            self.Coulomb = 1

class SIUnitSet(UnitSet):
    def __init__(self):
        super(SIUnitSet, self).__init__({'meter':1, 'sec':1, 'kg':1, 'Kelvin':1, 'Coulomb':1})

class RandomUnitSet(UnitSet):
    def __init__(self):
        from numpy.random import exponential
        super(RandomUnitSet, self).__init__({'meter':exponential(), 'sec':exponential(), 'kg':exponential(), 'Kelvin':exponential(), 'Coulomb':exponential()})

if __name__=='__main__':
    unit = SIUnitSet()
    exp = 149597870700.0, 1.9885998275439966e+30
    com = unit.AU, unit.solar_mass
    for i in range(len(exp)):
        assert exp[i] == com[i], "Expected %s, computed %s." % (exp[i], com[i])
