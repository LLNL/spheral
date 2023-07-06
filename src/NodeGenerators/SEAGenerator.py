from math import *

from NodeGeneratorBase import *

from Spheral import Vector3d
from Spheral import Tensor3d
from Spheral import SymTensor3d

from Spheral import vector_of_int, vector_of_double, vector_of_SymTensor3d, vector_of_vector_of_double

import mpi
rank = mpi.rank
procs = mpi.procs

class SEAGenerator3d(NodeGeneratorBase):
    def __init__(self,n,densityProfileMethod,
                 rmin = 0.0,
                 rmax = 0.0,
                 thetaMin = 0.0,
                 thetaMax = pi,
                 phiMin = 0.0,
                 phiMax = 2.0*pi,
                 nNodePerh = 2.01,
                 offset = None,
                 rejecter = None,
                 m0 = 0.0):

        assert n > 0
        assert rmin < rmax
        assert thetaMin < thetaMax
        assert thetaMin >= 0.0 and thetaMin <= 2.0*pi
        assert thetaMax >= 0.0 and thetaMax <= 2.0*pi
        assert phiMin < phiMax
        assert phiMin >= 0.0 and phiMin <= 2.0*pi
        assert phiMax >= 0.0 and phiMax <= 2.0*pi
        assert nNodePerh > 0.0
        assert offset is None or len(offset)==3
        
        self.rejecter = None
        if rejecter:
            self.rejecter = rejecter
        
        import random
        
        if offset is None:
            self.offset = Vector3d(0,0,0)
        else:
            self.offset = Vector3d(offset[0],offset[1],offset[2])
        
        self.n = n
        self.rmin = rmin
        self.rmax = rmax
        self.thetaMin = thetaMin
        self.thetaMax = thetaMax
        self.phiMin = phiMin
        self.phiMax = phiMax
        self.nNodePerh = nNodePerh        
        
        # If the user provided a constant for rho, then use the constantRho
        # class to provide this value.
        if type(densityProfileMethod) == type(1.0):
            self.densityProfileMethod = ConstantRho(densityProfileMethod)
        else:
            self.densityProfileMethod = densityProfileMethod
        
        # Determine how much total mass there is in the system.
        self.totalMass = self.integrateTotalMass(self.densityProfileMethod,
                                                     rmin, rmax,
                                                     thetaMin, thetaMax,
                                                     phiMin, phiMax)

        print("Total mass of %g in the range r = (%g, %g), theta = (%g, %g), phi = (%g, %g)" % \
            (self.totalMass, rmin, rmax, thetaMin, thetaMax, phiMin, phiMax))

        # Now set the nominal mass per node.
        if (m0 == 0.0):
            self.m0 = float(self.totalMass/n)
        else:
            self.m0 = float(m0)
            n = int(self.totalMass/self.m0)
        assert self.m0 > 0.0
        print("Nominal mass per node of %f for %d nodes." % (self.m0,n))
        
        from Spheral import SymTensor3d
        self.x = []
        self.y = []
        self.z = []
        self.m = []
        self.H = []

        # first shell is a tetrahedron
        self.positions = []
        nshell = 4
        rhoc = self.densityProfileMethod(0.0)
        print(rhoc,self.m0)
        mi   = self.m0
        drc  = pow(0.333*mi/rhoc,1.0/3.0)
        print(drc)
        hi   = nNodePerh*(drc)

        random.seed(nshell)
        dt = random.random()*pi
        dt2 = random.random()*pi
            
        rot = [[1.0,0.0,0.0],[0.0,cos(dt),-sin(dt)],[0.0,sin(dt),cos(dt)]]
        rot2 = [[cos(dt2),0.0,sin(dt2)],[0.0,1.0,0.0],[-sin(dt2),0.0,cos(dt2)]]

        t = sqrt(3.0)/3.0
        p1 = self.rotater([t,t,t],rot,rot2)
        p2 = self.rotater([t,-t,-t],rot,rot2)
        p3 = self.rotater([-t,-t,t],rot,rot2)
        p4 = self.rotater([-t,t,-t],rot,rot2)

        ps = [p1,p2,p3,p4]
        
        for pn in ps:
            self.positions.append(pn)
        mi = self.densityProfileMethod(drc*0.5)*(4.0/3.0*3.14159*drc**3)/4.0

        for pn in ps:
            x = pn[0]*0.5*drc
            y = pn[1]*0.5*drc
            z = pn[2]*0.5*drc

            if rejecter:
                if rejecter.accpet(x,y,z):
                    self.x.append(x)
                    self.y.append(y)
                    self.z.append(z)
                    self.m.append(mi)
                    self.H.append(SymTensor3d.one*(1.0/hi))
            else:
                self.x.append(x)
                self.y.append(y)
                self.z.append(z)
                self.m.append(mi)
                self.H.append(SymTensor3d.one*(1.0/hi))

        # now march up to the surface
        r = 0
        dr = drc
        while r <= rmax:
            r += dr
                
    
        # If requested, shift the nodes.
        if offset:
            for i in range(len(self.x)):
                self.x[i] += offset[0]
                self.y[i] += offset[1]
                self.z[i] += offset[2]
            
        print("Generated a total of %i nodes." % mpi.allreduce(len(self.x),mpi.SUM))
        NodeGeneratorBase.__init__(self, False,
                                   self.x, self.y, self.z, self.m, self.H)
        return

    def rotater(self,pos,rot1,rot2):
        posp = [0,0,0]
        for k in range(3):
            for j in range(3):
                posp[k] += pos[j]*rot1[k][j]
                
        x = posp[0]
        y = posp[1]
        z = posp[2]
        
        pos = [x,y,z]
        posp= [0,0,0]
        for k in range(3):
            for j in range(3):
                posp[k] += pos[j]*rot2[k][j]
        return posp

    #---------------------------------------------------------------------------
    # Compute the number of vertices for a given shape at a specific refinement
    # level.
    #  new formula for calculating number of points for a given subdivision level
    #  (Nf * Np(n) - Ne * Npe(n) + Nc)
    #  Nf = Number of faces of primitive shape
    #  Np(n) = Number of points in a triangle subdivided n times
    #       2^(2n-1) + 3*2^(n-1) + 1
    #  Ne = Number of edges of primitive shape
    #  Npe(n) = Number of points along an edge of primitive shape subdivided n times
    #       2^n + 1
    #  Nc = Number of corners
    #---------------------------------------------------------------------------
    def shapeCount(self, refinement, shape):
        Nf  = shape[0]
        Ne  = shape[1]
        Nc  = shape[2]
        n   = refinement
    
        Npe = 2**n + 1
        Np  = 2**(2*n-1) + 3*(2**(n-1)) + 1
        return (Nf * Np - Ne * Npe + Nc)
    
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
                           rmin, rmax,
                           thetaMin, thetaMax,
                           phiMin, phiMax,
                           nbins = 10000):
        assert nbins > 0
        assert nbins % 2 == 0
        
        result = 0
        dr = (rmax-rmin)/nbins
        nbp = nbins/procs
        binmin = nbp*rank if (rank!=0) else (1)
        binmax = nbp*(rank+1) if (rank!=procs-1) else (nbins)
        for i in range(binmin,binmax):
            r1 = rmin + (i-1)*dr
            r2 = rmin + i*dr
            result += 0.5*dr*(r2*r2*densityProfileMethod(r2)+r1*r1*densityProfileMethod(r1))
        result = result * (phiMax-phiMin) * (cos(thetaMin)-cos(thetaMax))
        result = mpi.allreduce(result,mpi.SUM)
        return result
