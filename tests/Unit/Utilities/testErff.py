from Spheral import *
from SpheralTestUtilities import *

import unittest
import Gnuplot

#===============================================================================
# Test the erff and erffc methods.
#===============================================================================
class TestErrorFunctions(unittest.TestCase):

    #===========================================================================
    # Set up, create arrays of the function values.
    #===========================================================================
    def setUp(self):
        self.minx = -3.0
        self.maxx = 3.0
        self.n = 1000
        self.x = [self.minx + (self.maxx - self.minx)/(self.n - 1)*i
                  for i in xrange(self.n)]
        self.erf = [erff(x) for x in self.x]
        self.erfc = [erffc(x) for x in self.x]

        return

    #===========================================================================
    # Plot the values for good measure.
    #===========================================================================
    def testPlotValues(self):
        self.erfData = Gnuplot.Data(self.x, self.erf,
                                    title = "Erf(x)")
        self.erfcData = Gnuplot.Data(self.x, self.erfc,
                                     title = "Erfc(x)")
        self.plot = Gnuplot.Gnuplot(persist = True)
        self.plot.replot(self.erfData)
        self.plot.replot(self.erfcData)

        return

    #===========================================================================
    # Check that every value is positive.
    #===========================================================================
    def testSign(self):
        for (x, erf, erfc) in zip(self.x, self.erf, self.erfc):
            self.failUnless(x*erf >= 0.0,
                            "Wrong sign for Erf:  %f %f" % (x, erf))
        return

    #===========================================================================
    # Check that Erf + Erfc == 1.
    #===========================================================================
    def testSumUnity(self):
        for (x, erf, erfc) in zip(self.x, self.erf, self.erfc):
            self.failUnless(fuzzyEqual(erf + erfc, 1.0),
                            "Erf(x) + Erfc(x) != 1.0: %f %f %f" % (erf, erfc, x))
        return

if __name__ == "__main__":
    unittest.main()

