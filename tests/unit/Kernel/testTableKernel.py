#ATS:test(SELF, label="TableKernel unit tests")
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
    # kernelExtent
    #===========================================================================
    def testWextent(self):
        for W, WT in self.kernelPairs:
            assert W.kernelExtent == WT.kernelExtent

    #===========================================================================
    # inflectionPoint
    #===========================================================================
    def testWinflectionPoint(self):
        for W, WT in self.kernelPairs:
            assert W.inflectionPoint == WT.inflectionPoint

    #===========================================================================
    # Check that the table kernel accurately reflects the real kernels values.
    #===========================================================================
    def testW_simultaneous_vals(self):
        for WT, Vector, SymTensor in ((self.WT1, Vector1d, SymTensor1d),
                                      (self.WT2, Vector2d, SymTensor2d),
                                      (self.WT3, Vector3d, SymTensor3d)):
            deta = WT.kernelExtent/(self.nsamples - 1)
            for i in xrange(self.nsamples):
                eta = i*deta
                H = 0.5*SymTensor.one
                Hdet = H.Determinant()
                Wval, gradWval, deltaWsum = WT.kernelAndGrad(Vector.zero, Vector(eta), H)
                self.failUnless(fuzzyEqual(Wval, WT.kernelValue(eta, Hdet), self.W0tol), 
                                "TableKernel value does match single lookup within tolerance:  %g %g" %
                                (Wval, WT.kernelValue(eta, Hdet)))
                self.failUnless(fuzzyEqual((gradWval - H*Vector(eta).unitVector()*WT.gradValue(eta, Hdet)).magnitude(), 0.0, self.W1tol),
                                "TableKernel grad value does match single lookup within tolerance:  %s %s" %
                                (gradWval, H*Vector(eta).unitVector()*WT.gradValue(eta, Hdet)))

    #===========================================================================
    # Check that the table kernel accurately reflects the real kernels values.
    #===========================================================================
    def testWlookup(self):
        for W, WT in self.kernelPairs:
            deta = W.kernelExtent/(self.nsamples - 1)
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

    #===========================================================================
    # Confirm that the nperh and Wsum values are monotonic.
    #===========================================================================
    def testMonotonicity(self):
        for W in self.tableKernels:
            nperh = W.nperhValues
            Wsum = W.WsumValues
            assert len(nperh) == len(Wsum)
            for i in xrange(len(nperh) - 1):
                self.failUnless(nperh[i] < nperh[i + 1],
                                "Failed monotonicity test in nperh table: %f %f" %
                                (nperh[i], nperh[i + 1]))
                self.failUnless(Wsum[i] <= Wsum[i + 1],
                                "Failed monotonicity test in Wsum table: %f %f" %
                                (Wsum[i], Wsum[i + 1]))
        return

    #===========================================================================
    # Check that the W sum values are reasonable for the 1D kernel.
    #===========================================================================
    def testWsumValues1d(self):
        W = self.WT1
        assert len(W.nperhValues) == len(W.WsumValues)
        for nperh, Wsum in zip(W.nperhValues, W.WsumValues):
            if Wsum > 0.0:
                deta = 1.0/nperh
                etax = deta
                testSum = 0.0
                while etax < W.kernelExtent:
                    delta = 2.0*abs(W.gradValue(etax, 1.0))
                    testSum += delta
                    etax += deta
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
        assert len(W.nperhValues) == len(W.WsumValues)
        for nperh, Wsum in zip(W.nperhValues, W.WsumValues):
            if Wsum > 0.0:
                deta = 1.0/nperh
                testSum = 0.0
                etay = 0.0
                while etay < W.kernelExtent:
                    etax = 0.0
                    while etax < W.kernelExtent:
                        eta = Vector2d(etax, etay)
                        delta = abs(W.gradValue(eta.magnitude(), 1.0))
                        if etax > 0.0:
                            delta *= 2.0
                        if etay > 0.0:
                            delta *= 2.0
                        testSum += delta
                        etax += deta
                    etay += deta
                testSum = sqrt(testSum)
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
        assert len(W.nperhValues) == len(W.WsumValues)
        for nperh, Wsum in zip(W.nperhValues, W.WsumValues):
            if Wsum > 0.0:
                deta = 1.0/nperh
                testSum = 0.0
                etaz = 0.0
                while etaz < W.kernelExtent:
                    etay = 0.0
                    while etay < W.kernelExtent:
                        etax = 0.0
                        while etax < W.kernelExtent:
                            eta = Vector3d(etax, etay, etaz)
                            delta = abs(W.gradValue(eta.magnitude(), 1.0))
                            if etax > 0.0:
                                delta *= 2.0
                            if etay > 0.0:
                                delta *= 2.0
                            if etaz > 0.0:
                                delta *= 2.0
                            testSum += delta
                            etax += deta
                        etay += deta
                    etaz += deta
                testSum = testSum**(1.0/3.0)
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

