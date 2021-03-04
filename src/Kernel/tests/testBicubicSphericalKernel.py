from math import *
import unittest
import numpy as np

from Spheral1d import *
from SpheralTestUtilities import fuzzyEqual

# Build the SphericalKernel
WT = TableKernel3d(BSplineKernel3d(), 500)
t0 = time.time()
W = SphericalTableKernel(WT)
t1 = time.time()
etamax = W.etamax
print("Required %0.4f sec to construct SphericalTableKernel"% (t1 - t0))

# The r/h distance from the origin for point i (the central point we're probing around)
etavals_i = ((i + 1)*0.1 for i in xrange(100))

# The values of h we'll consider
hvals_i = (0.1, 0.5, 1.0, 2.0, 5.0)

#-------------------------------------------------------------------------------
# The analytic form of the quadratic bi-cubic spline from Omang et al.
#-------------------------------------------------------------------------------
def W3S1(rj, ri, h):
    def C(q):
        return q*q - 0.75*q**4 + 0.3*q**5
    def D(q):
        return 2.0*(q*q - q**3) + 0.75*q**4 - 0.1*q**5
    sigj = rj/h
    sigi = ri/h
    sigdiff = abs(sigj - sigi)
    sigplus = sigj + sigi
    result = 0.0
    if sigplus <= 1.0:
        result = C(sigplus) - C(sigdiff)
    elif sigplus <= 2.0:
        if sigdiff < 1.0:
            result = -0.1 + D(sigplus) - C(sigdiff)
        elif sigdiff < 2.0:
            result = D(sigplus) - D(sigdiff)
    else:
        if sigdiff < 1.0:
            result = 0.7 - C(sigdiff)
        elif sigdiff < 2.0:
            result = 0.8 - D(sigdiff)
    return result/(h*rj*ri)

#-------------------------------------------------------------------------------
# The analytic gradient of the quadratic bi-cubic spline from Omang et al.
#-------------------------------------------------------------------------------
def gradW3S1(rj, ri, h):
    def C(q):
        return q*q - 0.75*q**4 + 0.3*q**5
    def D(q):
        return 2.0*(q*q - q**3) + 0.75*q**4 - 0.1*q**5
    def gradC(q):
        return 2.0*q - 3.0*q**3 + 1.5*q**4
    def gradD(q):
        return 4.0*q - 6.0*q**2 + 3.0*q**3 - 0.5*q**4
    sigj = rj/h
    sigi = ri/h
    sigdiff = abs(sigj - sigi)
    sigplus = sigj + sigi

    # \partial_rj
    if sigj > sigi:
        sgnfac = -1.0
    else:
        sgnfac = 1.0
    #sgnfac = 1.0

    if sigplus <= 1.0:
        return -W3S1(rj, ri, h)/rj + (gradC(sigplus) - sgnfac*gradC(sigdiff))/(h*h*ri*rj)

    elif sigplus <= 2.0:
        if sigdiff < 1.0:
            return -W3S1(rj, ri, h)/rj + (gradD(sigplus) - sgnfac*gradC(sigdiff))/(h*h*ri*rj)
        else:
            return -W3S1(rj, ri, h)/rj + (gradD(sigplus) - sgnfac*gradD(sigdiff))/(h*h*ri*rj)

    else:
        if sigdiff < 1.0:
            return -W3S1(rj, ri, h)/rj - sgnfac*gradC(sigdiff)/(h*h*ri*rj)
        else:
            return -W3S1(rj, ri, h)/rj - sgnfac*gradD(sigdiff)/(h*h*ri*rj)

#-------------------------------------------------------------------------------
# Return a useful r_j range for a given r_i
#-------------------------------------------------------------------------------
def rprange(r, h, etastep=0.05):
    return h*np.arange(max(0.01, r/h - 2.0), r/h + 2.0, etastep)

#-------------------------------------------------------------------------------
# Test tolerances
#-------------------------------------------------------------------------------
def Wtol(etai, etaj, etamax):
    delta = abs(etai - etaj)/etamax
    if delta < 0.5:
        return 1e-6
    if delta < 0.75:
        return 1e-5
    else:
        return 0.1

#-------------------------------------------------------------------------------
# Error estimator
#-------------------------------------------------------------------------------
def error(val0, val1, fuzz=1.0e-8):
    return abs(val1 - val0)/max(fuzz, abs(val0))

#-------------------------------------------------------------------------------
# Unit tests for the SphericalKernel class using an ordinary bicubic spline base
# kernel
#-------------------------------------------------------------------------------
class TestSphericalKernel(unittest.TestCase):

    def test_kernel_vals(self):
        for hi in hvals_i:
            for etai in etavals_i:
                ri = hi*etai
                for rj in rprange(ri, hi):
                    W0 = W3S1(rj, ri, hi)
                    Wij = W(Vector(ri/hi), Vector(rj/hi), 1.0/hi)
                    self.failUnless(error(W0, Wij) < Wtol(ri/hi, rj/hi, etamax),
                                    "Kernel value outside tolerance @ (%g, %g, %g): %g != %g, %g !<= %g" % (ri, rj, hi, W0, Wij, error(W0, Wij), Wtol(ri/hi, rj/hi, etamax)))


if __name__ == "__main__":
    unittest.main()
