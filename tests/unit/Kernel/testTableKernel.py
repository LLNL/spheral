#ATS:test(SELF, label="TableKernel unit tests")
from Spheral import *
from SpheralTestUtilities import fuzzyEqual

from math import *
import numpy as np
import unittest
import random

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
        self.kernelPairs = list(zip(self.realKernels, self.tableKernels))

        self.nsamples = 1000
        self.W0tol = 1.0e-3
        self.W1tol = 1.0e-2
        self.W2tol = 1.0e-2
        self.Wsumtol = 1.0e-1
        
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
            for i in range(self.nsamples):
                eta = i*deta
                H = 0.5*SymTensor.one
                Hdet = H.Determinant()
                Wval, gradWval, deltaWsum = WT.kernelAndGrad(Vector.zero, Vector(eta), H)
                self.assertTrue(fuzzyEqual(Wval, WT.kernelValue(eta, Hdet), self.W0tol), 
                                "TableKernel value does match single lookup within tolerance:  %g %g" %
                                (Wval, WT.kernelValue(eta, Hdet)))
                self.assertTrue(fuzzyEqual((gradWval - H*Vector(eta).unitVector()*WT.gradValue(eta, Hdet)).magnitude(), 0.0, self.W1tol),
                                "TableKernel grad value does match single lookup within tolerance:  %s %s" %
                                (gradWval, H*Vector(eta).unitVector()*WT.gradValue(eta, Hdet)))

    #===========================================================================
    # Check that the table kernel accurately reflects the real kernels values.
    #===========================================================================
    def testWlookup(self):
        for W, WT in self.kernelPairs:
            deta = W.kernelExtent/(self.nsamples - 1)
            for i in range(self.nsamples):
                eta = i*deta
                self.assertTrue(fuzzyEqual(WT.kernelValue(eta, 1.0),
                                           W.kernelValue(eta, 1.0), self.W0tol),
                                "TableKernel value does match original within tolerance:  %g %g" %
                                (WT.kernelValue(eta, 1.0), W.kernelValue(eta, 1.0)))
                self.assertTrue(fuzzyEqual(WT.gradValue(eta, 1.0),
                                           W.gradValue(eta, 1.0), self.W1tol),
                                "TableKernel grad value does match original within tolerance:  %g %g" %
                                (WT.gradValue(eta, 1.0), W.gradValue(eta, 1.0)))

    #===========================================================================
    # Confirm that the nperh and Wsum values are monotonic.
    #===========================================================================
    def testMonotonicity(self):
        for W in self.tableKernels:
            WsumMin = W.equivalentWsum(W.minNperhLookup)
            WsumMax = W.equivalentWsum(W.maxNperhLookup)
            n = 2*W.numPoints

            nperh = np.array([W.equivalentNodesPerSmoothingScale(x) for x in np.linspace(WsumMin, WsumMax, n)])
            self.assertTrue(np.all(np.diff(nperh)) > 0.0, "nperh lookup values not monotonic")

            Wsum = np.array([W.equivalentWsum(x) for x in np.linspace(W.minNperhLookup, W.maxNperhLookup, n)])
            self.assertTrue(np.all(np.diff(Wsum)) > 0.0, "Wsum lookup values not monotonic")
        return

    #===========================================================================
    # Check that the W sum values are reasonable for the 1D kernel.
    #===========================================================================
    def testWsumValues1d(self):
        W = self.WT1
        n = 2*W.numPoints
        minNperh = max(W.minNperhLookup, 0.5*W.kernelExtent)
        for nperh in np.linspace(minNperh, W.maxNperhLookup, n):
            deta = 1.0/nperh
            etax = deta
            testSum = 0.0
            while etax < W.kernelExtent:
                testSum += 2.0*abs(W.gradValue(etax, 1.0))
                etax += deta
            tol = self.Wsumtol / (W.kernelExtent/deta)
            self.assertTrue(fuzzyEqual(W.equivalentWsum(nperh), testSum, tol),
                            "Wsum failure: %g != %g @ %g: " %
                            (W.equivalentWsum(nperh), testSum, nperh))
            self.assertTrue(fuzzyEqual(W.equivalentNodesPerSmoothingScale(testSum),
                                       nperh,
                                       tol),
                            "Lookup n per h failure: %g %g %g" %
                            (testSum, W.equivalentNodesPerSmoothingScale(testSum), nperh))
        return

    #===========================================================================
    # Check that the W sum values are reasonable for the 2D kernel.
    #===========================================================================
    def testWsumValues2d(self):
        W = self.WT2
        minNperh = max(W.minNperhLookup, 0.5*W.kernelExtent)
        for itest in range(10):
            nperh = random.uniform(minNperh, W.maxNperhLookup)
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
            tol = self.Wsumtol / (W.kernelExtent/deta)**2
            self.assertTrue(fuzzyEqual(W.equivalentWsum(nperh), testSum, tol),
                            "Wsum failure: %g != %g @ %g: " %
                            (W.equivalentWsum(nperh), testSum, nperh))
            self.assertTrue(fuzzyEqual(W.equivalentNodesPerSmoothingScale(testSum),
                                       nperh,
                                       tol),
                            "Lookup n per h failure: %g %g %g" %
                            (testSum, W.equivalentNodesPerSmoothingScale(testSum), nperh))
                             
        return

    #===========================================================================
    # Check that the W sum values are reasonable for the 3D kernel.
    #===========================================================================
    def testWsumValues3d(self):
        W = self.WT3
        minNperh = max(W.minNperhLookup, 0.5*W.kernelExtent)
        for itest in range(10):
            nperh = random.uniform(minNperh, W.maxNperhLookup)
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
            tol = 5.0*self.Wsumtol / (W.kernelExtent/deta)**3
            self.assertTrue(fuzzyEqual(W.equivalentWsum(nperh), testSum, tol),
                            "Wsum failure: %g != %g @ %g: " %
                            (W.equivalentWsum(nperh), testSum, nperh))
            self.assertTrue(fuzzyEqual(W.equivalentNodesPerSmoothingScale(testSum),
                                       nperh,
                                       tol),
                            "Lookup n per h failure: %g %g %g" %
                            (testSum, W.equivalentNodesPerSmoothingScale(testSum), nperh))
                             
        return

if __name__ == "__main__":
    unittest.main()

