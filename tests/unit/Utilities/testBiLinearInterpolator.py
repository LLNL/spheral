#ATS:test(SELF, label="BiLinearInterpolator unit tests")

from Spheral import *
from SpheralTestUtilities import *

import unittest
import numpy as np
from numpy.polynomial.polynomial import polyval2d as npP2D
import numpy.random

import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

# Create a global random number generator.
import random
rangen = random.Random()

#===============================================================================
# A generator for creating a range of x values to test
#===============================================================================
def xygen(n, xmin, xmax, ymin, ymax):
    yield xmin, ymin
    yield xmax, ymin
    yield xmax, ymax,
    yield xmin, ymax
    count = 4
    while count < n:
        yield rangen.uniform(xmin, xmax), rangen.uniform(ymin, ymax)
        count += 1

#===============================================================================
# A polynomial functor
#===============================================================================
class PolynomialFunctor(ScalarScalarScalarFunctor):
    def __init__(self, order, cmin, cmax):
        ScalarScalarScalarFunctor.__init__(self)
        self.order = order
        self.c = np.random.rand(order + 1, order + 1)
        self.c *= cmax - cmin
        self.c += cmin
        return

    def __call__(self, x, y):
        return npP2D(x, y, self.c)

#===============================================================================
# TestLinearInterpolator
#===============================================================================
class TestBiLinearInterpolator(unittest.TestCase):

    #===========================================================================
    # Set up
    #===========================================================================
    def setUp(self):
        self.ntests = 1
        self.n = 10
        return

    #===========================================================================
    # Interpolate a linear function (should be exact)
    #===========================================================================
    def test_interp_linear(self):
        for (nx, ny) in ((2, 2), (3, 3), (10, 10), (3, 20)):
            for itest in xrange(self.ntests):
                F = PolynomialFunctor(1, -10.0, 10.0)
                xmin, ymin = -100.0, -100.0
                xmax, ymax =  100.0,  100.0
                Finterp = BiLinearInterpolator(xmin, xmax, ymin, ymax, nx, ny, F)
                # print "Finterp.coeffs: ", list(Finterp.coeffs)
                # print "poly.c        : ", F.c
                # assert fuzzyEqual(Finterp.coeffs[0], F.c[0][0])
                # assert fuzzyEqual(Finterp.coeffs[1], F.c[1][0]/(xmax - xmin))
                # assert fuzzyEqual(Finterp.coeffs[3], F.c[1][1]/((xmax - xmin)*(ymax - ymin)))
                # assert fuzzyEqual(Finterp.coeffs[2], F.c[0][1]/(ymax - ymin))
                for x, y in xygen(self.n, xmin, xmax, ymin, ymax):
                    self.failUnless(fuzzyEqual(Finterp(x, y), F(x, y)),
                                    "Interpolation off: %g != %g, err=%g" % (Finterp(x, y), F(x, y), abs(Finterp(x,y) - F(x,y))/(abs(F(x,y)) + 1e-10)))

    #===========================================================================
    # Interpolate a quadratic function (not exact)
    #===========================================================================
    def test_interp_quadratic(self):
        for (nx, ny) in ((3, 3), (10, 10), (3, 20)):
            for itest in xrange(self.ntests):
                F = PolynomialFunctor(2, -10.0, 10.0)
                xmin, ymin = -100.0, -100.0
                xmax, ymax =  100.0,  100.0
                Finterp = BiLinearInterpolator(xmin, xmax, ymin, ymax, nx, ny, F)

                # Compute a tolerance based on the range of the function
                z0 = np.array([F(xmin, ymin), F(xmax, ymin), F(xmax, ymax), F(xmin, ymax)])
                tol = 1.0e-3*(z0.max() - z0.min())
                for x, y in xygen(self.n, xmin, xmax, ymin, ymax):
                    self.failUnless(fuzzyEqual(Finterp(x, y), F(x, y), tol),
                                    "Interpolation off: %g != %g, Err %g" % (Finterp(x, y), F(x, y), abs(Finterp(x, y) - F(x, y))/max(1e-10, F(x, y))))

    #===========================================================================
    # Interpolate a cubic function (not exact)
    #===========================================================================
    def test_interp_cubic(self):
        for (nx, ny, acc) in ((20, 20, 10.0),
                              (100, 100, 5.0e-1)):
            for itest in xrange(self.ntests):
                F = PolynomialFunctor(3, -10.0, 10.0)
                xmin, ymin = -100.0, -100.0
                xmax, ymax =  100.0,  100.0
                Finterp = BiLinearInterpolator(xmin, xmax, ymin, ymax, nx, ny, F)
                # sys.stderr.write("(%i, %i): %s\n" % (nx, ny, Finterp))
                
                # Compute a tolerance based on the range of the function
                z0 = np.array([F(xmin, ymin), F(xmax, ymin), F(xmax, ymax), F(xmin, ymax)])
                tol = 1.0e-3*(z0.max() - z0.min())
                for x, y in xygen(self.n, xmin, xmax, ymin, ymax):
                    self.failUnless(fuzzyEqual(Finterp(x, y), F(x, y), tol),
                                    "Interpolation off @ (%g,%g): %g != %g, err=%g" %
                                    (x, y, Finterp(x, y), F(x, y), 2.0*abs((F(x, y) - Finterp(x, y))/(F(x, y) + Finterp(x, y)))) + "\nCoefficients: " + str(F.c))


if __name__ == "__main__":
    unittest.main()
