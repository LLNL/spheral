#ATS:test(SELF, label="BiQuadraticInterpolator unit tests")

from Spheral2d import *
from SpheralTestUtilities import *

import unittest
import numpy as np
import numpy.polynomial.polynomial.polyval2d as npP2D

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
        self.ntests = 100
        self.n = 100
        return

    #===========================================================================
    # A quadratic answer object
    #===========================================================================
    class QuadraticFunctor(VectorScalarFunctor2d):
        def __init__(self,
                     c00, c01, c02,
                     c10, c11, c12
                     c20, c21, c22):
            VectorScalarFunctor2d.__init__(self)
            self.c = np.array([c00, c01, c02],
                              [c10, c11, c12],
                              [c20, c21, c22])
            return

        def __call__(self, pos):
            return npP2D(pos.x, pos.y, self.c)

    #===========================================================================
    # Interpolate a quadratic function with linear spacing (should be exact)
    #===========================================================================
    def test_interp_quadratic_linear(self):
        for (nx, ny) in ((3, 3), (10, 10), (3, 20)):
            for itest in xrange(self.ntests):
                F = self.QuadraticFunctor(rangen.uniform(-10.0, 10.0),
                                          rangen.uniform(-10.0, 10.0),
                                          rangen.uniform(-10.0, 10.0),
                                          rangen.uniform(-10.0, 10.0),
                                          rangen.uniform(-10.0, 10.0),
                                          rangen.uniform(-10.0, 10.0),
                                          rangen.uniform(-10.0, 10.0),
                                          rangen.uniform(-10.0, 10.0),
                                          rangen.uniform(-10.0, 10.0))
                xmin = Vector(-100.0, -100.0)
                xmax = Vector( 100.0,  100.0)
                Finterp = BiQuadraticInterpolator(xmin, xmax, nx, ny, F)
                assert fuzzyEqual(Finterp.coeffs[0], F.c[0][0])
                assert fuzzyEqual(Finterp.coeffs[1], F.c[0][1])
                assert fuzzyEqual(Finterp.coeffs[2], F.c[0][2])
                assert fuzzyEqual(Finterp.coeffs[3], F.c[1][0])
                assert fuzzyEqual(Finterp.coeffs[4], F.c[1][1])
                assert fuzzyEqual(Finterp.coeffs[5], F.c[1][2])
                assert fuzzyEqual(Finterp.coeffs[6], F.c[2][0])
                assert fuzzyEqual(Finterp.coeffs[7], F.c[2][1])
                assert fuzzyEqual(Finterp.coeffs[8], F.c[2][2])
                for i in xrange(self.n):
                    pos = Vector(rangen.uniform(xmin.x, xmax.x),
                                 rangen.uniform(xmin.y, xmax.y))
                    self.failUnless(fuzzyEqual(Finterp(pos), F(pos)),
                                    "Interpolation off: %g != %g" % (Finterp(pos), F(pos)))


if __name__ == "__main__":
    unittest.main()
