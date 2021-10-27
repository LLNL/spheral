#ATS:test(SELF, label="QuadraticInterpolator unit tests")

from Spheral import *
from SpheralTestUtilities import *

import unittest

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
# TestQuadraticInterpolator
#===============================================================================
class TestQuadraticInterpolator(unittest.TestCase):

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
        yvals = vector_of_double([self.y(self.xmin + i*self.dx) for i in xrange(self.n)])
        self.F = QuadraticInterpolator(self.xmin, self.xmax, yvals)
        self.fuzz = 1.0e-10
        return

    #===========================================================================
    # Compute y(x) answer
    #===========================================================================
    def y(self, x):
        return self.A + self.B*x + self.C*x*x

    #===========================================================================
    # Compute dy/dx(x) answer
    #===========================================================================
    def yprime(self, x):
        return self.B + 2.0*self.C*x

    #===========================================================================
    # Compute d^2y/dx^2(x) answer
    #===========================================================================
    def yprime2(self, x):
        return 2.0*self.C

    #===========================================================================
    # Interpolate y
    #===========================================================================
    def test_yinterp(self):
        for x in xgen(self.ntests, self.xmin, self.xmax):
            self.failUnless(fuzzyEqual(self.F(x), self.y(x), self.fuzz),
                            "Error interpolating F(x): %g != %g" % (self.F(x), self.y(x)))

    #===========================================================================
    # Interpolate dy/dx
    #===========================================================================
    def test_dyinterp(self):
        for x in xgen(self.ntests, self.xmin, self.xmax):
            self.failUnless(fuzzyEqual(self.yprime(x), self.F.prime(x), self.fuzz),
                            "Error interpolating F'(x): %g != %g" % (self.F.prime(x), self.yprime(x)))

    #===========================================================================
    # Interpolate d^2y/dx^2
    #===========================================================================
    def test_dy2interp(self):
        for x in xgen(self.ntests, self.xmin, self.xmax):
            self.failUnless(fuzzyEqual(self.yprime2(x), self.F.prime2(x), self.fuzz),
                            "Error interpolating F''(x): %g != %g" % (self.F.prime2(x), self.yprime2(x)))

if __name__ == "__main__":
    unittest.main()
