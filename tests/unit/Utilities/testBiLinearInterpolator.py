#ATS:test(SELF, label="BiLinearInterpolator unit tests")

from Spheral import *
from SpheralTestUtilities import *
from XYInterpolatorTestingBase import *
import unittest

#===============================================================================
# TestLinearInterpolator
#===============================================================================
class TestBiLinearInterpolator(XYInterpolatorTestingBase,
                               unittest.TestCase):

    #===========================================================================
    # setUp
    #===========================================================================
    def setUp(self):
        self.order = 1
        self.tol = {1 : 1e-9,
                    2 : 100.0,
                    3 : 100.0}
        self.ntests = 20
        self.n = 100
        self.xmin, self.ymin = -100.0, -100.0
        self.xmax, self.ymax =  100.0,  100.0
        return

    #===========================================================================
    # Generate a BiLinearInterpolator
    #===========================================================================
    def generateInterpolator(self,
                             nx,
                             ny,
                             xlog,
                             ylog,
                             func):
        return BiLinearInterpolator(xmin = self.xmin,
                                    xmax = self.xmax,
                                    ymin = self.ymin,
                                    ymax = self.ymax,
                                    nx = nx,
                                    ny = ny,
                                    F = func,
                                    xlog = xlog,
                                    ylog = ylog)

#===============================================================================
# Run the tests
#===============================================================================
if __name__ == "__main__":
    unittest.main()
