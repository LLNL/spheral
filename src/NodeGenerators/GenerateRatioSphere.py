from math import *
import mpi

from NodeGeneratorBase import *

from Spheral import SimpsonsIntegrationDoubleFunction, simpsonsIntegrationDouble
from Spheral import Vector2d, Tensor2d, SymTensor2d, CylindricalBoundary, rotationMatrix2d
from Spheral import Vector3d, Tensor3d, SymTensor3d, CylindricalBoundary, rotationMatrix3d
from Spheral import CylindricalBoundary, generateCylDistributionFromRZ
from Spheral import vector_of_int, vector_of_double, vector_of_vector_of_double, vector_of_SymTensor3d

#-------------------------------------------------------------------------------
# This version ratios from the center out in 2D.  Kind of a misnomer with the
# sphere thing...
#-------------------------------------------------------------------------------
class GenerateRatioSphere2d(NodeGeneratorBase):
    
    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def __init__(self,
                 drCenter, drRatio,
                 rho,
                 rmin,
                 rmax,
                 thetamin = 0.0,
                 thetamax = 0.5*pi,
                 ntheta = 1,
                 center = (0.0, 0.0),
                 distributionType = "constantDTheta",   # one of (constantDTheta, constantNTheta)
                 nNodePerh = 2.01,
                 SPH = False,
                 rejecter = None):

        assert drCenter > 0.0
        assert drRatio > 0.0
        assert nNodePerh > 0.0
        assert rmin >= 0.0
        assert rmax > rmin
        assert thetamax > thetamin
        assert distributionType.lower() in ("constantdtheta", "constantntheta")
        
        self.center = center

        # Did we get passed a function or a constant for the density?
        if type(rho) == type(1.0):
            def rhofunc(ri):
                return rho
        else:
            rhofunc = rho
        self.rhofunc = rhofunc

        self.x, self.y, self.m, self.H = [], [], [], []

        constantN = (distributionType.lower() == "constantntheta")
        Dtheta = thetamax - thetamin

        # Decide the actual drCenter we're going to use to arrive at an integer number of radial bins.
        if abs(drRatio - 1.0) > 1e-4:
            neff = max(1, int(log(1.0 - (rmax - rmin)*(1.0 - drRatio)/drCenter)/log(drRatio) + 0.5))
            drCenter = (rmax - rmin)*(1.0 - drRatio)/(1.0 - drRatio**neff)
        else:
            neff = max(1, int((rmax - rmin)/drCenter + 0.5))
            drCenter = (rmax - rmin)/neff
        print "Adjusting initial radial spacing to %g in order to create an integer radial number of bins %i." % (drCenter, neff)

        # Work our way out from the center.
        r0 = rmin
        dr = drCenter
        for i in xrange(neff):
            if abs(drRatio - 1.0) > 1e-4:
                r0 = rmin + drCenter*(1.0 - drRatio**i)/(1.0 - drRatio)
                r1 = rmin + drCenter*(1.0 - drRatio**(i + 1))/(1.0 - drRatio)
            else:
                r0 = rmin + i*drCenter
                r1 = rmin + (i + 1)*drCenter
            dr = r1 - r0
            ri = 0.5*(r0 + r1)
            li = Dtheta*ri
            if constantN:
                ntheta = ntheta
            else:
                ntheta = max(1, int(li/dr))
            dtheta = Dtheta/ntheta
            mring = (r1**2 - r0**2) * 0.5*Dtheta * rhofunc(ri)
            mi = mring/ntheta
            hr = nNodePerh * dr
            ha = nNodePerh * ri*dtheta
            for j in xrange(ntheta):
                thetai = thetamin + (j + 0.5)*dtheta
                xi = ri*cos(thetai)
                yi = ri*sin(thetai)
                self.x.append(xi + center[0])
                self.y.append(yi + center[1])
                self.m.append(mi)
                if SPH:
                    hi = sqrt(hr*ha)
                    self.H.append(SymTensor2d(1.0/hi, 0.0, 0.0, 1.0/hi))
                else:
                    self.H.append(SymTensor2d(1.0/hr, 0.0, 0.0, 1.0/ha))
                    runit = Vector2d(xi, yi).unitVector()
                    T = rotationMatrix2d(runit).Transpose()
                    self.H[-1].rotationalTransform(T)

        # Do a numerical integral to get the expected total mass.
        class integfunc(SimpsonsIntegrationDoubleFunction):
            def __init__(self, rho, Dtheta):
                SimpsonsIntegrationDoubleFunction.__init__(self)
                self.rho = rho
                self.Dtheta = Dtheta
                return
            def __call__(self, ri):
                return Dtheta*ri*self.rho(ri)
        M1 = simpsonsIntegrationDouble(integfunc(rhofunc, Dtheta), rmin, rmax, 10000)

        # Make sure the total mass is what we intend it to be, by applying
        # a multiplier to the particle masses.
        M0 = sum(self.m)
        assert M0 > 0.0
        massCorrection = M1/M0
        for i in xrange(len(self.m)):
            self.m[i] *= massCorrection
        print "Applied a mass correction of %f to ensure total mass is %f." % (massCorrection, M1)

        # If the user provided a "rejecter", give it a pass
        # at the nodes.
        if rejecter:
            self.x, self.y, self.z, self.m, self.H = rejecter(self.x,
                                                              self.y,
                                                              self.m,
                                                              self.H)

        # Have the base class break up the serial node distribution
        # for parallel cases.
        NodeGeneratorBase.__init__(self, True,
                                   self.x, self.y, self.m, self.H)
        return

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localPosition(self, i):
        assert i >= 0 and i < len(self.x)
        assert len(self.x) == len(self.y)
        return Vector2d(self.x[i], self.y[i])
    
    #---------------------------------------------------------------------------
    # Get the mass for the given node index.
    #---------------------------------------------------------------------------
    def localMass(self, i):
        assert i >= 0 and i < len(self.m)
        return self.m[i]
    
    #---------------------------------------------------------------------------
    # Get the mass density for the given node index.
    #---------------------------------------------------------------------------
    def localMassDensity(self, i):
        ri = sqrt((self.x[i] - self.center[0])**2 + (self.y[i] - self.center[1])**2)
        return self.rhofunc(ri)
    
    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        assert i >= 0 and i < len(self.H)
        return self.H[i]
    
#--------------------------------------------------------------------------------
# 3D version, actual sphere.
# Based on spinning the case above.
#--------------------------------------------------------------------------------
class GenerateRatioSphere3d(NodeGeneratorBase):

    def __init__(self, 
                 drCenter, drRatio,
                 rho,
                 rmin,
                 rmax,
                 thetamin = 0.0,
                 thetamax = 0.5*pi,
                 phi = pi,
                 ntheta = 1,
                 center = (0.0, 0.0, 0.0),
                 distributionType = "constantDTheta",   # one of (constantDTheta, constantNTheta)
                 nNodePerh = 2.01,
                 SPH = False,
                 rejecter = None):

        assert thetamax <= pi

        self.gen2d = GenerateRatioSphere2d(drCenter, drRatio, rho, rmin, rmax, thetamin, thetamax, ntheta, 
                                           (0.0, 0.0), distributionType, nNodePerh, SPH)

        # The 2D class already split the nodes up between processors, but
        # we want to handle that ourselves.  Distribute the full set of RZ
        # nodes to every process, then redecompose them below.
        self.x = mpi.allreduce(self.gen2d.x[:], mpi.SUM)
        self.y = mpi.allreduce(self.gen2d.y[:], mpi.SUM)
        self.m = mpi.allreduce(self.gen2d.m[:], mpi.SUM)
        self.H = mpi.allreduce(self.gen2d.H[:], mpi.SUM)
        n = len(self.x)
        self.z = [0.0]*n
        self.globalIDs = [0]*n

        # Convert the 2-D H tensors to 3-D, and correct the masses.
        for i in xrange(n):
            xi = self.x[i]
            yi = self.y[i]
            H2d = SymTensor2d(self.H[i])
            H2dinv = H2d.Inverse()

            hxy0 = 0.5*(H2dinv.Trace())
            dphi = CylindricalBoundary.angularSpacing(yi, hxy0, nNodePerh, 2.0)
            assert dphi > 0.0
            nsegment = max(1, int(phi/dphi + 0.5))
            dphi = phi/nsegment

            hz = dphi*yi*nNodePerh
            self.H[i] = SymTensor3d(H2d.xx, H2d.xy, 0.0,
                                    H2d.yx, H2d.yy, 0.0,
                                    0.0,    0.0,    1.0/hz)
            if SPH:
                h0 = self.H[i].Determinant()**(1.0/3.0)
                self.H[-1] = SymTensor3d(h0, 0.0, 0.0,
                                         0.0, h0, 0.0,
                                         0.0, 0.0, h0)

            # Convert the mass to the full hoop mass, which will then be used in
            # generateCylDistributionFromRZ to compute the actual nodal masses.
            mi = self.m[i]
            circ = 2.0*pi*yi
            mhoop = mi*circ
            self.m[i] = mhoop

        assert len(self.m) == n
        assert len(self.H) == n

        # Duplicate the nodes from the xy-plane, creating rings of nodes about
        # the x-axis.  We use a C++ helper method for the sake of speed.
        kernelExtent = 2.0
        extras = []
        xvec = self.vectorFromList(self.x, vector_of_double)
        yvec = self.vectorFromList(self.y, vector_of_double)
        zvec = self.vectorFromList(self.z, vector_of_double)
        mvec = self.vectorFromList(self.m, vector_of_double)
        Hvec = self.vectorFromList(self.H, vector_of_SymTensor3d)
        globalIDsvec = self.vectorFromList(self.globalIDs, vector_of_int)
        extrasVec = vector_of_vector_of_double()
        for extra in extras:
            extrasVec.append(self.vectorFromList(extra, vector_of_double))
        generateCylDistributionFromRZ(xvec, yvec, zvec, mvec, Hvec, globalIDsvec,
                                      extrasVec,
                                      nNodePerh, kernelExtent, phi,
                                      mpi.rank, mpi.procs)
        self.x = [x + center[0] for x in xvec]
        self.y = [x + center[1] for x in yvec]
        self.z = [x + center[2] for x in zvec]
        self.m = list(mvec)
        self.H = [SymTensor3d(x) for x in Hvec]
        self.globalIDs = list(globalIDsvec)
        for i in xrange(len(extras)):
            extras[i] = list(extrasVec[i])

        self.center = Vector3d(*center)

        # If the user provided a "rejecter", give it a pass
        # at the nodes.
        if rejecter:
            self.x, self.y, self.z, self.m, self.H = rejecter(self.x,
                                                              self.y,
                                                              self.z,
                                                              self.m,
                                                              self.H)
        # Initialize the base class.
        NodeGeneratorBase.__init__(self, False,
                                   self.x, self.y, self.z, self.m, self.H)
        return

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localPosition(self, i):
        return Vector3d(self.x[i], self.y[i], self.z[i])

    #---------------------------------------------------------------------------
    # Get the mass for the given node index.
    #---------------------------------------------------------------------------
    def localMass(self, i):
        return self.m[i]

    #---------------------------------------------------------------------------
    # Get the mass density for the given node index.
    #---------------------------------------------------------------------------
    def localMassDensity(self, i):
        return self.gen2d.rhofunc((self.localPosition(i) - self.center).magnitude())

    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        return self.H[i]

