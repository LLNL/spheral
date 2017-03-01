from math import *
import unittest

from Spheral import *
from SpheralTestUtilities import fuzzyEqual

#===============================================================================
# Main testing class.
#===============================================================================
class TestTableKernel(unittest.TestCase):

    #===========================================================================
    # Initialization for each test.
    #===========================================================================
    def setUp(self):

        self.W1 = BSplineKernel1d()
        self.WT1 = TableKernel1d(self.W1)

        self.W2 = BSplineKernel2d()
        self.WT2 = TableKernel2d(self.W2)

        self.W3 = BSplineKernel3d()
        self.WT3 = TableKernel3d(self.W3)

        self.realKernels = [self.W1, self.W2, self.W3]
        self.tableKernels = [self.WT1, self.WT2, self.WT3]
        self.kernelPairs = zip(self.realKernels, self.tableKernels)

        self.nsamples = 1000
        self.W0tol = 1.0e-3
        self.W1tol = 1.0e-2
        self.W2tol = 1.0e-2
        self.Wsumtol = 1.0e-10
        
        return

    #===========================================================================
    # Check that the table kernel accurately reflects the real kernels values.
    #===========================================================================
    def testWlookup(self):
        for W, WT in self.kernelPairs:
            deta = W.kernelExtent()/(self.nsamples - 1)
            for i in xrange(self.nsamples):
                eta = i*deta
                self.failUnless(fuzzyEqual(WT.kernelValue(eta, 1.0),
                                           W.kernelValue(eta, 1.0), self.W0tol),
                                "TableKernel value does match original within tolerance:  %g %g" %
                                (WT.kernelValue(eta, 1.0), W.kernelValue(eta, 1.0)))
                self.failUnless(fuzzyEqual(WT.gradValue(eta, 1.0),
                                           W.gradValue(eta, 1.0), self.W1tol),
                                "TableKernel grad value does match original within tolerance:  %g %g" %
                                (WT.gradValue(eta, 1.0), W.gradValue(eta, 1.0)))
                self.failUnless(fuzzyEqual(WT.grad2Value(eta, 1.0),
                                           W.grad2Value(eta, 1.0), self.W2tol),
                                "TableKernel grad2 value does match original within tolerance:  %g %g" %
                                (WT.grad2Value(eta, 1.0), W.grad2Value(eta, 1.0)))

    #===========================================================================
    # Confirm that the nperh and Wsum values are monotonic.
    #===========================================================================
    def testMonotonicity(self):
        for W in self.tableKernels:
            nperh = W.nperh()
            Wsum = W.Wsum()
            assert len(nperh) == len(Wsum)
            for i in xrange(len(nperh) - 1):
                self.failUnless(nperh[i] < nperh[i + 1],
                                "Failed monotonicity test in nperh table: %f %f" %
                                (nperh[i], nperh[i + 1]))
                self.failUnless(Wsum[i] < Wsum[i + 1],
                                "Failed monotonicity test in Wsum table: %f %f" %
                                (Wsum[i], Wsum[i + 1]))
        return

    #===========================================================================
    # Check that the W sum values are reasonable for the 1D kernel.
    #===========================================================================
    def testWsumValues1d(self):
        W = self.WT1
        assert len(W.nperh()) == len(W.Wsum())
        for nperh, Wsum in zip(W.nperh(), W.Wsum()):
            deta = 1.0/nperh
            nx = int(W.kernelExtent()*nperh) + 1
            testSum = 0.0
            for ix in xrange(nx):
                eta = ix*deta
                delta = W.kernelValue(eta, 1.0)
                if ix > 0:
                    delta *= 2.0
                testSum += delta
            self.failUnless(fuzzyEqual(Wsum, testSum, self.Wsumtol),
                            "Wsum failure: %g != %g: " %
                            (Wsum, testSum))
            self.failUnless(fuzzyEqual(W.equivalentNodesPerSmoothingScale(testSum),
                                       nperh,
                                       self.Wsumtol),
                            "Lookup n per h failure: %g %g %g" %
                            (testSum, W.equivalentNodesPerSmoothingScale(testSum), nperh))
        return

    #===========================================================================
    # Check that the W sum values are reasonable for the 2D kernel.
    #===========================================================================
    def testWsumValues2d(self):
        W = self.WT2
        assert len(W.nperh()) == len(W.Wsum())
        for nperh, Wsum in zip(W.nperh(), W.Wsum()):
            deta = 1.0/nperh
            nx = int(W.kernelExtent()*nperh) + 1
            testSum = 0.0
            for iy in xrange(nx):
                etay = iy*deta
                for ix in xrange(nx):
                    etax = ix*deta
                    eta = sqrt(etax*etax + etay*etay)
                    delta = W.kernelValue(eta, 1.0)
                    if ix > 0:
                        delta *= 2.0
                    if iy > 0:
                        delta *= 2.0
                    testSum += delta
            self.failUnless(fuzzyEqual(Wsum, testSum, self.Wsumtol),
                            "Wsum failure: %g != %g: " %
                            (Wsum, testSum))
            self.failUnless(fuzzyEqual(W.equivalentNodesPerSmoothingScale(testSum),
                                       nperh,
                                       self.Wsumtol),
                            "Lookup n per h failure: %g %g %g" %
                            (testSum, W.equivalentNodesPerSmoothingScale(testSum), nperh))
                             
        return

    #===========================================================================
    # Check that the W sum values are reasonable for the 3D kernel.
    #===========================================================================
    def testWsumValues3d(self):
        W = self.WT3
        assert len(W.nperh()) == len(W.Wsum())
        for nperh, Wsum in zip(W.nperh(), W.Wsum()):
            deta = 1.0/nperh
            nx = int(W.kernelExtent()*nperh) + 1
            testSum = 0.0
            for iz in xrange(nx):
                etaz = iz*deta
                for iy in xrange(nx):
                    etay = iy*deta
                    for ix in xrange(nx):
                        etax = ix*deta
                        eta = sqrt(etax*etax + etay*etay + etaz*etaz)
                        delta = W.kernelValue(eta, 1.0)
                        if ix > 0:
                            delta *= 2.0
                        if iy > 0:
                            delta *= 2.0
                        if iz > 0:
                            delta *= 2.0
                        testSum += delta
            self.failUnless(fuzzyEqual(Wsum, testSum, self.Wsumtol),
                            "Wsum failure: %g != %g: " %
                            (Wsum, testSum))
            self.failUnless(fuzzyEqual(W.equivalentNodesPerSmoothingScale(testSum),
                                       nperh,
                                       self.Wsumtol),
                            "Lookup n per h failure: %g %g %g" %
                            (testSum, W.equivalentNodesPerSmoothingScale(testSum), nperh))
                             
        return

if __name__ == "__main__":
    unittest.main()

