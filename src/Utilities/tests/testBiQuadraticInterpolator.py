#ATS:test(SELF, label="BiQuadraticInterpolator unit tests")

from Spheral import *
from SpheralTestUtilities import *

import unittest
import numpy as np
from numpy.polynomial.polynomial import polyval2d as npP2D

import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

# Create a global random number generator.
import random
rangen = random.Random()

#===============================================================================
# TestQuadraticInterpolator
#===============================================================================
class TestBiQuadraticInterpolator(unittest.TestCase):

    #===========================================================================
    # Set up
    #===========================================================================
    def setUp(self):
        self.ntests = 20
        self.n = 20
        return

    #===========================================================================
    # A quadratic answer object
    #===========================================================================
    class QuadraticFunctor(ScalarScalarScalarFunctor):
        def __init__(self,
                     c00, c01, c02,
                     c10, c11, 
                     c20):
            ScalarScalarScalarFunctor.__init__(self)
            self.c = np.array([[c00, c01, c02],
                               [c10, c11, 0.0],
                               [c20, 0.0, 0.0]])
            return

        def __call__(self, x, y):
            return npP2D(x, y, self.c)

    #===========================================================================
    # A cubic answer object
    #===========================================================================
    class CubicFunctor(ScalarScalarScalarFunctor):
        def __init__(self,
                     c00, c01, c02, c03,
                     c10, c11, c12, 
                     c20, c21,
                     c30):
            ScalarScalarScalarFunctor.__init__(self)
            self.c = np.array([[c00, c01, c02, c03],
                               [c10, c11, c12, 0.0],
                               [c20, c21, 0.0, 0.0],
                               [c30, 0.0, 0.0, 0.0]])
            return

        def __call__(self, x, y):
            return npP2D(x, y, self.c)

    #===========================================================================
    # Interpolate a quadratic function (should be exact)
    #===========================================================================
    def test_interp_quadratic(self):
        for (nx, ny) in ((3, 3), (10, 10), (3, 20)):
            for itest in xrange(self.ntests):
                F = self.QuadraticFunctor(rangen.uniform(-10.0, 10.0),
                                          rangen.uniform(-10.0, 10.0),
                                          rangen.uniform(-10.0, 10.0),
                                          rangen.uniform(-10.0, 10.0),
                                          rangen.uniform(-10.0, 10.0),
                                          rangen.uniform(-10.0, 10.0))
                xmin, ymin = -100.0, -100.0
                xmax, ymax =  100.0,  100.0
                Finterp = BiQuadraticInterpolator(xmin, xmax, ymin, ymax, nx, ny, F)
                assert fuzzyEqual(Finterp.coeffs[0], F.c[0][0])
                assert fuzzyEqual(Finterp.coeffs[1], F.c[1][0])
                assert fuzzyEqual(Finterp.coeffs[2], F.c[0][1])
                assert fuzzyEqual(Finterp.coeffs[3], F.c[1][1])
                assert fuzzyEqual(Finterp.coeffs[4], F.c[2][0])
                assert fuzzyEqual(Finterp.coeffs[5], F.c[0][2])
                for i in xrange(self.n):
                    x = rangen.uniform(xmin, xmax)
                    y = rangen.uniform(ymin, ymax)
                    self.failUnless(fuzzyEqual(Finterp(x, y), F(x, y)),
                                    "Interpolation off: %g != %g" % (Finterp(x, y), F(x, y)))

    #===========================================================================
    # Interpolate a cubic function (not exact)
    #===========================================================================
    def test_interp_cubic(self):
        for (nx, ny, acc) in ((20, 20, 10.0),
                              (100, 100, 5.0e-1)):
            for itest in xrange(self.ntests):
                F = self.CubicFunctor(rangen.uniform(-10.0, 10.0),
                                      rangen.uniform(-10.0, 10.0),
                                      rangen.uniform(-10.0, 10.0),
                                      rangen.uniform(-10.0, 10.0),
                                      rangen.uniform(-10.0, 10.0),
                                      rangen.uniform(-10.0, 10.0),
                                      rangen.uniform(-10.0, 10.0),
                                      rangen.uniform(-10.0, 10.0),
                                      rangen.uniform(-10.0, 10.0),
                                      rangen.uniform(-10.0, 10.0))
                xmin, ymin = -100.0, -100.0
                xmax, ymax =  100.0,  100.0
                Finterp = BiQuadraticInterpolator(xmin, xmax, ymin, ymax, nx, ny, F)
                # sys.stderr.write("(%i, %i): %s\n" % (nx, ny, Finterp))
                
                # # Plotting fun
                # x = np.arange(xmin.x, 1.001*xmax.x, (xmax.x - xmin.x)/100)
                # y = np.arange(xmin.y, 1.001*xmax.y, (xmax.y - xmin.y)/100)
                # x, y = np.meshgrid(x, y)
                # NX, NY = x.shape
                # z0 = np.array([[F(Vector(x[j][i], y[j][i])) for j in xrange(NY)] for i in xrange(NX)])
                # z1 = np.array([[Finterp(Vector(x[j][i], y[j][i])) for j in xrange(NY)] for i in xrange(NX)])
                # zdiff = z1 - z0

                # fig0 = plt.figure()
                # ax0 = fig0.add_subplot(111, projection='3d')
                # surf0 = ax0.plot_surface(x, y, z0, cmap=cm.coolwarm,
                #                          linewidth=0, antialiased=False)
                # fig1 = plt.figure()
                # ax1 = fig1.add_subplot(111, projection='3d')
                # surf1 = ax1.plot_surface(x, y, z1, cmap=cm.coolwarm,
                #                          linewidth=0, antialiased=False)
                # fig2 = plt.figure()
                # ax2 = fig2.add_subplot(111, projection='3d')
                # surf2 = ax2.plot_surface(x, y, zdiff, cmap=cm.coolwarm,
                #                          linewidth=0, antialiased=False)
                # plt.show()
                # # Plotting fun

                for i in xrange(self.n):
                    x = rangen.uniform(xmin, xmax)
                    y = rangen.uniform(ymin, ymax)
                    self.failUnless(fuzzyEqual(Finterp(x, y), F(x, y), acc),
                                    "Interpolation off: %g != %g, err=%g" % (Finterp(x, y), F(x, y), 2.0*abs((F(x, y) - Finterp(x, y))/(F(x, y) + Finterp(x, y)))))


if __name__ == "__main__":
    unittest.main()
