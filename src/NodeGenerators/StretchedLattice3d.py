from math import *

from NodeGeneratorBase import *

from Spheral import Vector3d
from Spheral import Tensor3d
from Spheral import SymTensor3d

from Spheral import vector_of_int, vector_of_double, vector_of_SymTensor3d, vector_of_vector_of_double

import mpi
procID = mpi.rank
nProcs = mpi.procs

#-------------------------------------------------------------------------------
# Class to generate 3-D node positions in a stretched lattice
#-------------------------------------------------------------------------------
class GenerateStretchedLattice3d(NodeGeneratorBase):

    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def __init__(self, nr, densityProfileMethod,
                 rmin = 0.0,
                 rmax = 1.0,
                 thetaMin = 0.0,
                 thetaMax = pi,
                 phiMin = 0.0,
                 phiMax = 2.0*pi,
                 nNodePerh = 2.01,
                 offset=None):
        
        assert nr > 0
        assert rmin >= 0
        assert rmin < rmax
        assert thetaMin < thetaMax
        assert thetaMin >= 0.0 and thetaMin <= 2.0*pi
        assert thetaMax >= 0.0 and thetaMax <= 2.0*pi
        assert phiMin < phiMax
        assert phiMin >= 0.0 and phiMin <= 2.0*pi
        assert phiMax >= 0.0 and phiMax <= 2.0*pi
        assert nNodePerh > 0.0
        assert offset is None or len(offset)==3
        
        if offset is None:
            self.offset = Vector3d(0,0,0)
        else:
            self.offset = Vector3d(offset[0],offset[1],offset[2])
        
        self.nr         = nr
        self.rmin       = rmin
        self.rmax       = rmax
        self.thetaMin   = thetaMin
        self.thetaMax   = thetaMax
        self.phiMin     = phiMin
        self.phiMax     = phiMax
        self.nNodePerh  = nNodePerh
        
        self.xmin       = Vector3d(-rmax,-rmax,-rmax)
        self.xmax       = Vector3d(rmax,rmax,rmax)
        
        # no reason to support a constant density method here, just use a regular lattice for that
        self.densityProfileMethod = densityProfileMethod
        
        # Determine how much total mass there is in the system.
        self.totalMass = self.integrateTotalMass(self.densityProfileMethod,
                                                 rmax)
        
        self.ntot       = 4.0/3.0*pi*(nr**3)
        self.m0         = self.totalMass/self.ntot
        self.vol        = 4.0/3.0*pi*(rmax**3)
        # what this means is this currently only supports creating a full sphere and then
        # cutting out the middle to rmin if rmin > 0
        self.rho0       = self.totalMass/self.vol
    
        # compute kappa first
        k = 3/(self.rho0*rmax**3) * self.totalMass/(4.0*pi)
        print "Found kappa={0:3.3f}. Was that what you expected?".format(k)
            
        # create the unstretched lattice
        self.x, self.y, self.z, self.m, self.H = \
            self.latticeDistribution(self.nr,
                                     self.rho0,
                                     self.m0,
                                     self.xmin,    # (xmin, ymin, zmin)
                                     self.xmax,    # (xmax, ymax, zmax)
                                     self.rmax,
                                     self.nNodePerh)

        nx  = 2*nr+1
        eta = (self.xmax[0] - self.xmin[0])/nx
        
        # Initialize the base class.  If "serialInitialization" is True, this
        # is where the points are broken up between processors as well.
        serialInitialization = True
        NodeGeneratorBase.__init__(self, serialInitialization,
                                   self.x, self.y, self.z, self.m, self.H)
    
        return
    
    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localPosition(self, i):
        assert i >= 0 and i < len(self.x)
        assert len(self.x) == len(self.y) == len(self.z)
        return Vector3d(self.x[i], self.y[i], self.z[i])
    
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
        loc = Vector3d(0,0,0)
        loc = self.localPosition(i) - self.offset
        return self.densityProfileMethod(loc.magnitude())
    
    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        assert i >= 0 and i < len(self.H)
        return self.H[i]


    #---------------------------------------------------------------------------
    # Numerically integrate the given density profile to determine the total
    # enclosed mass.
    #---------------------------------------------------------------------------
    def integrateTotalMass(self, densityProfileMethod,
                           rmax,
                           nbins = 10000):
        assert nbins > 0
        assert nbins % 2 == 0
        
        result = 0
        dr = rmax/nbins
        for i in xrange(1,nbins):
            r1 = (i-1)*dr
            r2 = i*dr
            result += 0.5*dr*(r2*r2*densityProfileMethod(r2)+r1*r1*densityProfileMethod(r1))
        result = result * 4.0*pi
        return result

    #---------------------------------------------------------------------------
    # Helper functions for newton raphson
    #---------------------------------------------------------------------------
    def func(self,x,k,rho0,r0,eta,rp):
        const = k*r0**(2)*rho0*eta
        return ((x-rp)*x**(2)*rho(x)-const)

    def dfunc(self,x,eta,rp,dr):
        fst = (rho(x) - rho(x-dr))/dr * (x-rp)*x*x
        scd = rho(x) * (x*x + (x-rp)*(2)*x*x)
        return fst + scd

    def rootFind(self,x,k,rho0,r0,eta,rp,dr,guess,nsteps=100):
        for i in xrange(ns):
            f = self.func(x = guess,
                          k = k,
                          rho0 = rho0,
                          r0 = r0,
                          eta = eta,
                          rp = rp)
            df = self.dfunc(x = guess,
                            eta = eta,
                            rp = rp,
                            dr = dr)
            guess = guess - f/df

    #-------------------------------------------------------------------------------
    # Seed positions/masses on the unstretched lattice
    #-------------------------------------------------------------------------------
    def latticeDistribution(self, nr, rho0, m0,
                            xmin,
                            xmax,
                            rmax,
                            nNodePerh = 2.01):
        
        assert nr > 0
        assert rho0 > 0
        
        nx = 2*nr+1
        ny = 2*nr+1
        nz = 2*nr+1
        
        dx = (xmax[0] - xmin[0])/nx
        dy = (xmax[1] - xmin[1])/ny
        dz = (xmax[2] - xmin[2])/nz
        
        n = nx*ny*nz
        print xmax,xmin
        print "nx={0:f} dx={1:3.3e}".format(nx,dx)

        hx = 1.0/(nNodePerh*dx)
        hy = 1.0/(nNodePerh*dy)
        hz = 1.0/(nNodePerh*dz)
        H0 = SymTensor3d(hx, 0.0, 0.0,
                         0.0, hy, 0.0,
                         0.0, 0.0, hz)
                       
        x = []
        y = []
        z = []
        m = []
        H = []

        imin, imax = self.globalIDRange(n)
        for iglobal in xrange(imin, imax):
            i = iglobal % nx
            j = (iglobal // nx) % ny
            k = iglobal // (nx*ny)
            xx = xmin[0] + (i + 0.5)*dx
            yy = xmin[1] + (j + 0.5)*dy
            zz = xmin[2] + (k + 0.5)*dz
            r2 = (xx)**2 + (yy)**2 + (zz)**2
            if (rmax is None or r2 <= rmax**2):
               x.append(xx)
               y.append(yy)
               z.append(zz)
               m.append(m0)
               H.append(H0)
       
        return x, y, z, m, H