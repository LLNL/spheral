#ATS:test(SELF, label="BiCubicInterpolator unit tests")

from Spheral import *
from SpheralTestUtilities import *
from XYInterpolatorTestingBase import *
import unittest

#===============================================================================
# TestCubicInterpolator
#===============================================================================
class TestBiCubicInterpolator(XYInterpolatorTestingBase,
                                  unittest.TestCase):

    #===========================================================================
    # setUp
    #===========================================================================
    def setUp(self):
        self.order = 3
        self.tol = {1 : 1e-10,
                    2 : 1e-6,
                    3 : 1e-6}
        self.ntests = 20
        self.n = 100
        self.xmin, self.ymin = -100.0, -100.0
        self.xmax, self.ymax =  100.0,  100.0
        return

    #===========================================================================
    # Generate a BiCubicInterpolator
    #===========================================================================
    def generateInterpolator(self,
                             nx,
                             ny,
                             xlog,
                             ylog,
                             func):
        gradFunc = GradPolynomialFunctor(func)
        return BiCubicInterpolator(xmin = self.xmin,
                                   xmax = self.xmax,
                                   ymin = self.ymin,
                                   ymax = self.ymax,
                                   nx = nx,
                                   ny = ny,
                                   F = func,
                                   #gradF = gradFunc,
                                   xlog = xlog,
                                   ylog = ylog)

#===============================================================================
# Run the tests
#===============================================================================
if __name__ == "__main__":
    unittest.main()

