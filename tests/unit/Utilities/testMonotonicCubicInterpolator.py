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
# Functors to generate the problem
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

class Fcubic(ScalarScalarFunctor):
    def __init__(self, A, B, C, D):
        ScalarScalarFunctor.__init__(self)
        self.A = A
        self.B = B
        self.C = C
        self.D = D
        return
    def __call__(self, x):
        return self.A + self.B*x + self.C*x*x + self.D*x*x*x
    def prime(self, x):
        return self.B + 2.0*self.C*x + 3.0*self.C*x*x
    def prime2(self, x):
        return 2.0*self.C + 6.0*self.C*x

class Fgrad(ScalarScalarFunctor):
    def __init__(self, F):
        ScalarScalarFunctor.__init__(self)
        self.F = F
        return
    def __call__(self, x):
        return self.F.prime(x)

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
        return

    #===========================================================================
    # Plotting fun
    #===========================================================================
    def plotem(self, x, xmin, xmax, F, Finterp):
        nx = 100
        x = np.arange(xmin, 1.001*xmax, (xmax - xmin)/nx)
        y0 = np.array([F(xi) for xi in x])
        y1 = np.array([Finterp(xi) for xi in x])
        ydiff = y1 - y0

        fig0 = plt.figure()
        plt.plot(x, y0, label="Analytic")
        plt.legend()

        fig1 = plt.figure()
        plt.plot(x, y1, label="Interpolation")
        plt.legend()

        fig2 = plt.figure()
        plt.plot(x, ydiff, label="Error")
        plt.legend()

        plt.show()

    #===========================================================================
    # Interpolate a quadratic function (func only)
    #===========================================================================
    def test_quad_interp(self):
        xmin = -10.0
        xmax =  40.0
        A = rangen.uniform(-100.0, 100.0)
        B = rangen.uniform(-100.0, 100.0)
        C = rangen.uniform(-100.0, 100.0)
        func = Fquad(A, B, C)
        F = MonotonicCubicInterpolator(xmin, xmax, self.n, func)
        tol = 1.0e-10
        for x in xgen(self.ntests, xmin, xmax):
            passing = fuzzyEqual(F(x), func(x), tol)
            if not passing:
                print F.vals
                self.plotem(x, xmin, xmax, func, F)
            self.failUnless(passing,
                            "Error interpolating F(x): %g != %g" % (F(x), func(x)))

    #===========================================================================
    # Interpolate a quadratic function (func + gradient)
    #===========================================================================
    def test_quad_interp_with_grad(self):
        xmin = -10.0
        xmax =  40.0
        A = rangen.uniform(-100.0, 100.0)
        B = rangen.uniform(-100.0, 100.0)
        C = rangen.uniform(-100.0, 100.0)
        func = Fquad(A, B, C)
        F = MonotonicCubicInterpolator(xmin, xmax, self.n, func, Fgrad(func))
        tol = 1.0e-10
        for x in xgen(self.ntests, xmin, xmax):
            passing = fuzzyEqual(F(x), func(x), tol)
            if not passing:
                print F.vals
                self.plotem(x, xmin, xmax, func, F)
            self.failUnless(passing,
                            "Error interpolating F(x): %g != %g" % (F(x), func(x)))

    #===========================================================================
    # Interpolate a cubic function (func only)
    #===========================================================================
    def test_cubic_interp(self):
        xmin = -10.0
        xmax =  40.0
        A = rangen.uniform(-100.0, 100.0)
        B = rangen.uniform(-100.0, 100.0)
        C = rangen.uniform(-100.0, 100.0)
        D = rangen.uniform(-100.0, 100.0)
        func = Fcubic(A, B, C, D)
        F = MonotonicCubicInterpolator(xmin, xmax, self.n, func)
        tol = 1.0e-10
        for x in xgen(self.ntests, xmin, xmax):
            passing = fuzzyEqual(F(x), func(x), tol)
            if not passing:
                print F.vals
                self.plotem(x, xmin, xmax, func, F)
            self.failUnless(passing,
                            "Error interpolating F(x): %g != %g, err = %g" % (F(x), func(x), abs(F(x) - func(x))/(1e-30 + abs(F(x)) + abs(func(x)))))

    #===========================================================================
    # Interpolate a cubic function (func + gradient)
    #===========================================================================
    def test_cubic_interp_with_grad(self):
        xmin = -10.0
        xmax =  40.0
        A = rangen.uniform(-100.0, 100.0)
        B = rangen.uniform(-100.0, 100.0)
        C = rangen.uniform(-100.0, 100.0)
        D = rangen.uniform(-100.0, 100.0)
        func = Fcubic(A, B, C, D)
        F = MonotonicCubicInterpolator(xmin, xmax, self.n, func, Fgrad(func))
        tol = 1.0e-10
        for x in xgen(self.ntests, xmin, xmax):
            passing = fuzzyEqual(F(x), func(x), tol)
            if not passing:
                print F.vals
                self.plotem(x, xmin, xmax, func, F)
            self.failUnless(passing,
                            "Error interpolating F(x): %g != %g, err = %g" % (F(x), func(x), abs(F(x) - func(x))/(1e-30 + abs(F(x)) + abs(func(x)))))

if __name__ == "__main__":
    unittest.main()
