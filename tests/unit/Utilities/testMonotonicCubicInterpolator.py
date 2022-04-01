#ATS:test(SELF, label="MonotonicCubicInterpolator unit tests")

from Spheral import *
from SpheralTestUtilities import *

import unittest
import matplotlib.pyplot as plt
import numpy as np

# Create a global random number generator.
import random
rangen = random.Random()

#===========================================================================
# A generator for creating a range of x values to test
#===========================================================================
def xgen(n, xmin, xmax):
    assert n > 2
    count = 0
    while count < n:
        if count == 0:
            yield xmin
        elif count == n - 1:
            yield xmax
        else:
            yield rangen.uniform(xmin, xmax)
        count += 1

#===============================================================================
# Functor to generate the problem
#===============================================================================
class Fquad(ScalarScalarFunctor):
    def __init__(self, A, B, C):
        ScalarScalarFunctor.__init__(self)
        self.A = A
        self.B = B
        self.C = C
        return
    def __call__(self, x):
        return self.A + self.B*x + self.C*x*x
    def prime(self, x):
        return self.B + 2.0*self.C*x
    def prime2(self, x):
        return 2.0*self.C

#===============================================================================
# TestMonotonicCubicInterpolator
#===============================================================================
class TestMonotonicCubicInterpolator(unittest.TestCase):

    #===========================================================================
    # Set up
    #===========================================================================
    def setUp(self):
        self.ntests = 100000
        self.n = 100
        self.xmin = -10.0
        self.xmax =  40.0
        self.dx = (self.xmax - self.xmin)/(self.n - 1)
        self.A = rangen.uniform(-100.0, 100.0)
        self.B = rangen.uniform(-100.0, 100.0)
        self.C = rangen.uniform(-100.0, 100.0)
        self.func = Fquad(self.A, self.B, self.C)
        self.F = MonotonicCubicInterpolator(self.xmin, self.xmax, self.n, self.func)
        self.fuzz = 1.0e-10
        return

    #===========================================================================
    # Plotting fun
    #===========================================================================
    def plotem(self, x, F, Finterp):
        nx = 100
        x = np.arange(self.xmin, 1.001*self.xmax, (self.xmax - self.xmin)/nx)
        y0 = np.array([F(xi) for xi in x])
        y1 = np.array([Finterp(xi) for xi in x])
        ydiff = y1 - y0

        fig0 = plt.figure()
        plt.plot(x, y0, label="Analytic")

        fig1 = plt.figure()
        plt.plot(x, y1, label="Interpolation")

        fig2 = plt.figure()
        plt.plot(x, ydiff, label="Error")

        plt.show()

    #===========================================================================
    # Interpolate y
    #===========================================================================
    def test_yinterp(self):
        for x in xgen(self.ntests, self.xmin, self.xmax):
            passing = fuzzyEqual(self.F(x), self.func(x), self.fuzz)
            if not passing:
                self.plotem(x, self.F, self.func)
            self.failUnless(passing,
                            "Error interpolating F(x): %g != %g" % (self.F(x), self.func(x)))

if __name__ == "__main__":
    unittest.main()
