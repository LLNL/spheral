#ATS:test(SELF, label="BiQuadraticInterpolator unit tests")

from Spheral import *
from SpheralTestUtilities import *

import unittest
import numpy as np
from numpy.polynomial.polynomial import polyval2d as npP2D

import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

from testBiLinearInterpolator import xygen, PolynomialFunctor

# Create a global random number generator.
import random
rangen = random.Random()

#===========================================================================
# A generator for creating a range of x values to test
#===========================================================================
def xygen(n, xmin, xmax, ymin, ymax):
    assert n > 4
    count = 0
    while count < n:
        if count == 0:
            yield xmin, ymin
        elif count == 1:
            yield xmax, ymin
        elif count == 2:
            yield xmax, ymax,
        elif count == 3:
            yield xmin, ymax
        else:
            yield rangen.uniform(xmin, xmax), rangen.uniform(ymin, ymax)
        count += 1

#===============================================================================
# TestQuadraticInterpolatorBase
#===============================================================================
class TestBiQuadraticInterpolatorBase:

    #===========================================================================
    # Set up
    #===========================================================================
    def setUp(self):
        self.ntests = 20
        self.n = 1000
        self.xmin, self.ymin = -100.0, -100.0
        self.xmax, self.ymax =  100.0,  100.0
        return

    #===========================================================================
    # Analytic coordinate of a grid point
    #===========================================================================
    def coord_ans(self, dx, ix, nx, xlog):
        assert ix <= nx
        if xlog:
            return self.xmin + dx*exp(ix - nx)
        else:
            return self.xmin + ix*dx

    #===========================================================================
    # Analytic normalized coordinate
    #===========================================================================
    def eta_ans(self, x, dx, ix, nx, xlog):
        x0 = self.coord_ans(dx, ix, nx, xlog)
        x1 = self.coord_ans(dx, ix + 1, nx, xlog)
        return max(0.0, min(1.0, (x - x0)/(x1 - x0)))

    #===========================================================================
    # Analytic lowerBound
    #===========================================================================
    def lowerBound_ans(self, x, y, dx, dy, nx, ny, xlog, ylog):
        if xlog:
            ix = nx - 1 + int(log(min(1, max(1e-10, (x - self.xmin)/dx))))
        else:
            ix = int((x - self.xmin)/dx)
        if ylog:
            iy = ny - 1 + int(log(min(1, max(1e-10, (y - self.ymin)/dy))))
        else:
            iy = int((y - self.ymin)/dy)
        ix = min(nx - 1, max(0, ix))
        iy = min(ny - 1, max(0, iy))
        i0 = 9*(nx*iy + ix)
        return ix, iy, i0

    #===========================================================================
    # lowerBound
    #===========================================================================
    def test_lowerBound(self):
        for (nx, ny) in ((2, 2), (3, 3), (10, 10), (3, 20)):
            F = PolynomialFunctor(2, -10.0, 10.0)
            Finterp = self.generateInterpolator(nx, ny, F)
            for itest in xrange(self.ntests):
                x = rangen.uniform(self.xmin, self.xmax)
                y = rangen.uniform(self.ymin, self.ymax)
                ix0, iy0, i0 = self.lowerBound_ans(x, y,
                                                   Finterp.xstep,
                                                   Finterp.ystep,
                                                   nx - 1, ny - 1,
                                                   Finterp.xlog, Finterp.ylog)
                ix, iy, i = Finterp.lowerBound(x, y)
                self.failUnless((ix == ix0) and (iy == iy0) and (i == i0),
                                "(nx,ny) = (%i,%i), (x,y) = (%g,%g)\n  lowerBound lookup: %i %i %i\n  lowerBound answer: %i %i %i" % (nx, ny, x, y, ix, iy, i, ix0, iy0, i0))

    #===========================================================================
    # Plotting fun
    #===========================================================================
    def plotem(self, x, y, F, Finterp):
        x = np.arange(self.xmin, 1.001*self.xmax, (self.xmax - self.xmin)/100)
        y = np.arange(self.ymin, 1.001*self.ymax, (self.ymax - self.ymin)/100)
        x, y = np.meshgrid(x, y)
        NX, NY = x.shape
        z0 = np.array([[F(x[j][i], y[j][i]) for j in xrange(NY)] for i in xrange(NX)])
        z1 = np.array([[Finterp(x[j][i], y[j][i]) for j in xrange(NY)] for i in xrange(NX)])
        zdiff = z1 - z0

        fig0 = plt.figure()
        ax0 = fig0.add_subplot(111, projection='3d')
        surf0 = ax0.plot_surface(x, y, z0, cmap=cm.coolwarm,
                                 linewidth=0, antialiased=False)
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111, projection='3d')
        surf1 = ax1.plot_surface(x, y, z1, cmap=cm.coolwarm,
                                 linewidth=0, antialiased=False)
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111, projection='3d')
        surf2 = ax2.plot_surface(x, y, zdiff, cmap=cm.coolwarm,
                                 linewidth=0, antialiased=False)
        plt.show()

    #===========================================================================
    # Interpolate a quadratic function (should be exact)
    #===========================================================================
    def test_interp_quadratic(self):
        for (nx, ny) in ((2, 2), (3, 3), (10, 10), (3, 20)):
            for itest in xrange(self.ntests):
                F = PolynomialFunctor(2, -10.0, 10.0)
                Finterp = self.generateInterpolator(nx, ny, F)

                # print "Coeffs: ", list(Finterp.coeffs)
                # print "Ans   : ", F.c

                # Should match a quadratic function near exactly
                for x, y in xygen(self.n, self.xmin, self.xmax, self.ymin, self.ymax):
                    if Finterp.xlog:
                        x = max(x, self.coord_ans(Finterp.xstep, 0, nx - 1, True))
                    if Finterp.ylog:
                        y = max(y, self.coord_ans(Finterp.ystep, 0, ny - 1, True))
                    passing = fuzzyEqual(Finterp(x, y), F(x, y))
                    if not passing:
                        self.plotem(x, y, F, Finterp)
                    self.failUnless(passing,
                                    "Interpolation off @(%g, %g) with grid (%i, %i): %g != %g" % (x, y, nx, ny, Finterp(x, y), F(x, y)))

    #===========================================================================
    # Interpolate a cubic function (not exact)
    #===========================================================================
    def test_interp_cubic(self):
        for (nx, ny) in ((2, 2), (3, 3), (10, 10), (3, 20)):
            for itest in xrange(self.ntests):
                F = PolynomialFunctor(3, -10.0, 10.0)
                Finterp = self.generateInterpolator(nx, ny, F)
                # sys.stderr.write("(%i, %i): %s\n" % (nx, ny, Finterp))
                
                # Compute a tolerance based on the range of the function
                z0 = np.array([F(self.xmin, self.ymin), F(self.xmax, self.ymin), F(self.xmax, self.ymax), F(self.xmin, self.ymax)])
                tol = 1.0e-6*(z0.max() - z0.min())

                for x, y in xygen(self.n, self.xmin, self.xmax, self.ymin, self.ymax):
                    if Finterp.xlog:
                        x = max(x, self.coord_ans(Finterp.xstep, 0, nx - 1, True))
                    if Finterp.ylog:
                        y = max(y, self.coord_ans(Finterp.ystep, 0, ny - 1, True))
                    passing = fuzzyEqual(Finterp(x, y), F(x, y), tol)
                    if not passing:
                        self.plotem(x, y, F, Finterp)
                    self.failUnless(passing,
                                    "Interpolation off @ (%g,%g): %g != %g, err=%g" %
                                    (x, y, Finterp(x, y), F(x, y), 2.0*abs((F(x, y) - Finterp(x, y))/(F(x, y) + Finterp(x, y)))) + "\nCoefficients: " + str(F.c))


#===============================================================================
# TestQuadraticInterpolatorLinearSpacing
#===============================================================================
class TestBiQuadraticInterpolatorLinearSpacing(TestBiQuadraticInterpolatorBase,
                                               unittest.TestCase):

    #===========================================================================
    # Generate a BiQuadraticInterpolator
    #===========================================================================
    def generateInterpolator(self,
                             nx,
                             ny,
                             func):
        return BiQuadraticInterpolator(xmin = self.xmin,
                                       xmax = self.xmax,
                                       ymin = self.ymin,
                                       ymax = self.ymax,
                                       nx = nx,
                                       ny = ny,
                                       F = func,
                                       xlog = False,
                                       ylog = False)

#===============================================================================
# TestQuadraticInterpolatorXlog
#===============================================================================
class TestBiQuadraticInterpolatorXlog(TestBiQuadraticInterpolatorBase,
                                      unittest.TestCase):

    #===========================================================================
    # Generate a BiQuadraticInterpolator
    #===========================================================================
    def generateInterpolator(self,
                             nx,
                             ny,
                             func):
        return BiQuadraticInterpolator(xmin = self.xmin,
                                       xmax = self.xmax,
                                       ymin = self.ymin,
                                       ymax = self.ymax,
                                       nx = nx,
                                       ny = ny,
                                       F = func,
                                       xlog = True,
                                       ylog = False)

#===============================================================================
# TestQuadraticInterpolatorYlog
#===============================================================================
class TestBiQuadraticInterpolatorYlog(TestBiQuadraticInterpolatorBase,
                                      unittest.TestCase):

    #===========================================================================
    # Generate a BiQuadraticInterpolator
    #===========================================================================
    def generateInterpolator(self,
                             nx,
                             ny,
                             func):
        return BiQuadraticInterpolator(xmin = self.xmin,
                                       xmax = self.xmax,
                                       ymin = self.ymin,
                                       ymax = self.ymax,
                                       nx = nx,
                                       ny = ny,
                                       F = func,
                                       xlog = False,
                                       ylog = True)

#===============================================================================
# TestQuadraticInterpolatorXYlog
#===============================================================================
class TestBiQuadraticInterpolatorXYlog(TestBiQuadraticInterpolatorBase,
                                       unittest.TestCase):

    #===========================================================================
    # Generate a BiQuadraticInterpolator
    #===========================================================================
    def generateInterpolator(self,
                             nx,
                             ny,
                             func):
        return BiQuadraticInterpolator(xmin = self.xmin,
                                       xmax = self.xmax,
                                       ymin = self.ymin,
                                       ymax = self.ymax,
                                       nx = nx,
                                       ny = ny,
                                       F = func,
                                       xlog = True,
                                       ylog = True)

#===============================================================================
# Run the tests
#===============================================================================
if __name__ == "__main__":
    unittest.main()

