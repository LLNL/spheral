#ATS:test(SELF, label="EquationOfState unit tests")

from SpheralTestUtilities import fuzzyEqual
import numpy as np
from math import *
import unittest

from Spheral1d import *

# Create a global random number generator.
import random
rangen = random.Random()

# We'll just use a gamma-law gas to base our tests on.
gamma = 5.0/3.0
mu = 5.0
units = CGS()
eos = GammaLawGas(gamma = gamma,
                  mu = mu,
                  constants = units)
gamma1 = gamma - 1.0

#-------------------------------------------------------------------------------
# Test harness.
#-------------------------------------------------------------------------------
class EOSTest(unittest.TestCase):

    def testSpecificThermalEnergyLookup(self):
        n = 1000000                      # How many tests?
        rhoMin, rhoMax = 1e-5, 1000.0    # range of densities
        Pmin, Pmax = 0.0, 1e5            # range of target pressures
        epsTol, Ptol = 1e-10, 1e-10      # tolerances
        Perrcheck = 1e-5                 # error check
        maxIterations = 100
        rho0 = np.random.uniform(rhoMin, rhoMax, n)
        P0 = np.random.uniform(Pmin, Pmax, n)

        # Now some cuteness to call a method across elements of numpy arrays.
        epsLookup = np.vectorize(eos.specificThermalEnergyForPressure)
        Plookup = np.vectorize(eos.pressure)

        # Find the estimated eps(P, rho)
        eps = epsLookup(P0, rho0, 0.0, 1e10, epsTol, Ptol, maxIterations)
        
        # The corresponding pressures.
        P = Plookup(rho0, eps)

        # Compute the error
        # reltol, abstol = 1e4*Ptol, 1e4*Ptol
        # def passfail(x, y):
        #     if abs(y) < abstol:
        #         if abs(x - y) > abstol:
        #             print "abs Blago: ", x, y
        #         return abs(x - y) < abstol
        #     else:
        #         if abs(x - y)/y > reltol:
        #             print "rel Blago: ", x, y
        #         return abs(x - y)/y < reltol
        # passfailLookup = np.vectorize(passfail)
        # Perr = np.minimum(P0, np.abs((P - P0)/np.maximum(1.0e-10, P0)))
        # self.failUnless(passfailLookup(P, P0).all(),
        #                 "Pressure error out of tolerance: %s vs %s" % (Perr.max(), Ptol))
        # assert passfailLookup(P, P0).all()
        Perr = np.minimum(P0, np.abs((P - P0)/np.maximum(1.0e-10, P0)))
        self.failUnless((Perr < Perrcheck).all(),
                        "Pressure error out of tolerance: %s > %s" % (Perr.max(), Perrcheck))
        return

#-------------------------------------------------------------------------------
# Run those tests.
#-------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
