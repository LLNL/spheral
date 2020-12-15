#ATS:test(SELF, label="BiQuadraticInterpolator unit tests")

from Spheral2d import *
from SpheralTestUtilities import *

import unittest

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
        self.ntests = 1000
        self.n = 10
        return

    #===========================================================================
    # A quadratic answer object
    #===========================================================================
    class BiQuadraticFunctor(VectorScalarFunctor2d):
        def __init__(self, c0, c1, c2, c3, c4, c5):
            VectorScalarFunctor2d.__init__(self)
            self.c0 = c0
            self.c1 = c1
            self.c2 = c2
            self.c3 = c3
            self.c4 = c4
            self.c5 = c5
            return

        def __call__(self, x):
            return (self.c0 +
                    self.c1 * x.x +
                    self.c2 * x.y +
                    self.c3 * x.x*x.y +
                    self.c4 * x.x*x.x +
                    self.c5 * x.y*x.y)

    #===========================================================================
    # Interpolate a quadratic function with linear spacing (should be exact)
    #===========================================================================
    def test_interp_quadratic_linear(self):
        for (nx, ny) in ((3, 3), (10, 10)):
            for itest in xrange(self.ntests):
                F = self.BiQuadraticFunctor(0.0, # rangen.uniform(-10.0, 10.0),
                                            rangen.uniform(-10.0, 10.0),
                                            rangen.uniform(-10.0, 10.0),
                                            rangen.uniform(-10.0, 10.0),
                                            rangen.uniform(-10.0, 10.0),
                                            rangen.uniform(-10.0, 10.0))
                xmin = Vector(-100.0, -100.0)
                xmax = Vector( 100.0,  100.0)
                Finterp = BiQuadraticInterpolator(xmin, xmax, nx, ny, False, False, F)
                print "--------------------------------------------------------------------------------"
                print "Input : ", [F.c0, F.c1, F.c2, F.c3, F.c4, F.c5]
                print "Interp: ", list(Finterp.coeffs)
                assert fuzzyEqual(Finterp.coeffs[0], F.c0)
                assert fuzzyEqual(Finterp.coeffs[1], F.c1)
                assert fuzzyEqual(Finterp.coeffs[2], F.c2)
                assert fuzzyEqual(Finterp.coeffs[3], F.c3)
                assert fuzzyEqual(Finterp.coeffs[4], F.c4)
                assert fuzzyEqual(Finterp.coeffs[5], F.c5)
                for i in xrange(self.n):
                    pos = Vector(rangen.uniform(xmin.x, xmax.x),
                                 rangen.uniform(xmin.y, xmax.y))
                    assert fuzzyEqual(Finterp(pos), F(pos))

    #===========================================================================
    # Interpolate a quadratic function with log spacing (should be exact)
    #===========================================================================
    # def test_interp_quadratic_log(self):
    #     for (nx, ny) in ((2, 2), (10, 10)):
    #         for itest in xrange(self.ntests):
    #             F = self.BiQuadraticFunctor(rangen.uniform(-10.0, 10.0),
    #                                         rangen.uniform(-10.0, 10.0),
    #                                         rangen.uniform(-10.0, 10.0),
    #                                         rangen.uniform(-10.0, 10.0),
    #                                         rangen.uniform(-10.0, 10.0),
    #                                         rangen.uniform(-10.0, 10.0))
    #             xmin = Vector(1e-10, 1e-10)
    #             xmax = Vector(100.0, 100.0)
    #             Finterp = BiQuadraticInterpolator(xmin, xmax, nx, ny, True, True, F)
    #             print "--------------------------------------------------------------------------------"
    #             print "Input : ", [F.c0, F.c1, F.c2, F.c3, F.c4, F.c5]
    #             print "Interp: ", list(Finterp.coeffs)
    #             assert fuzzyEqual(Finterp.coeffs[0], F.c0)
    #             assert fuzzyEqual(Finterp.coeffs[1], F.c1)
    #             assert fuzzyEqual(Finterp.coeffs[2], F.c2)
    #             assert fuzzyEqual(Finterp.coeffs[3], F.c3)
    #             assert fuzzyEqual(Finterp.coeffs[4], F.c4)
    #             assert fuzzyEqual(Finterp.coeffs[5], F.c5)
    #             for i in xrange(self.n):
    #                 pos = Vector(rangen.uniform(xmin.x, xmax.x),
    #                              rangen.uniform(xmin.y, xmax.y))
    #                 assert fuzzyEqual(Finterp(pos), F(pos))

if __name__ == "__main__":
    unittest.main()
