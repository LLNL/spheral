from math import *

from NodeGeneratorBase import *
#from Spheral import generateCylDistributionFromRZ

from Spheral import Vector2d, Vector3d
from Spheral import Tensor2d, Tensor3d
from Spheral import SymTensor2d, SymTensor3d
from Spheral import CylindricalBoundary
from Spheral import vector_of_double, vector_of_int, vector_of_SymTensor3d, vector_of_vector_of_double

import mpi
rank = mpi.rank
procs = mpi.procs

delimiter = ' '

#-------------------------------------------------------------------------------
# Read in a 2-D setup.
#-------------------------------------------------------------------------------
class AsciiFileNodeGenerator2D(NodeGeneratorBase):

    #---------------------------------------------------------------------------
    # Constructor.
    #---------------------------------------------------------------------------
    def __init__(self,
                 filename,
                 materialName,
                 nNodePerh = 2.01,
                 SPH = False,
                 precision = 20,
                 Hscalefactor = 1.0,
                 extraFields = [],
                 initializeBase = True,
                 readFileToMemory = False,
                 refineNodes = 0):
        
        
        self.filename = filename
        self.nPerh = nNodePerh
        self.SPH = SPH
        self.extraFields = extraFields
        
        # For now we restrict to reading from a single (serial) file.
        allfiles = mpi.allreduce([filename], mpi.SUM)
        assert min([x == filename for x in allfiles])
        self.serialfile = True
        
        self.H = []
        self.fieldNames = []
        vals = []
        self.nf = 0
        self.nv = 0
        
        print "using " + str(mpi.procs) +  " procs"
        
        if mpi.rank == 0:
            f = open(filename,'r')
            self.f = f
            
            gotFieldNames = 0
            
            for line in self.f:
                line.rstrip('\r\n')
                data = line.split(delimiter)
                if data[0][0] != "#" and gotFieldNames == 1:
                    vals.append(data)
                if data[0][0] != "#" and gotFieldNames == 0:
                    self.fieldNames.append(data)
                    gotFieldNames = 1
            
            self.f.close()
            
            for i in xrange(len(self.fieldNames[0])):
                print self.fieldNames[0][i]
                self.nf = self.nf + 1
        
        
        self.nf = mpi.bcast(self.nf, root=0)
        self.nv = mpi.bcast(len(vals), root=0)
        self.fieldNames = mpi.bcast(self.fieldNames, root=0)
        mpi.barrier()
        
        for i in xrange(len(self.fieldNames[0])):
            self.__dict__[self.fieldNames[0][i]] = []
        
        assert self.nf == len(self.fieldNames[0])
        
        if mpi.rank == 0:
            n = len(vals)
            for i in xrange(n):
                for j in xrange(len(self.fieldNames[0])):
                    self.__dict__[self.fieldNames[0][j]].append(float(vals[i][j]))
                self.H.append((1.0/self.h[i]) * SymTensor2d.one)
        
        for j in xrange(len(self.fieldNames[0])):
            self.__dict__[self.fieldNames[0][j]] = mpi.bcast(self.__dict__[self.fieldNames[0][j]], root=0)
        
        self.H = mpi.bcast(self.H, root=0)

        
        # Initialize the base class.
        if initializeBase:
            fields = tuple([self.x, self.y, self.m, self.rho, self.vx, self.vy, self.eps, self.H] +
                           [self.__dict__[x] for x in extraFields])
            NodeGeneratorBase.__init__(self, self.serialfile, *fields)

        
        # Apply the requested number of refinements.
        for i in xrange(refineNodes):
            refineNodes2d(self)
        
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
        assert i >= 0 and i < len(self.rho)
        return self.rho[i]

    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        assert i >= 0 and i < len(self.H)
        return self.H[i]



#-------------------------------------------------------------------------------
# Read in a 3-D setup.
#-------------------------------------------------------------------------------
class AsciiFileNodeGenerator3D(NodeGeneratorBase):


    #---------------------------------------------------------------------------
    # Constructor.
    #---------------------------------------------------------------------------
    def __init__(self,
                 filename,
                 materialName,
                 nNodePerh = 2.01,
                 SPH = False,
                 precision = 20,
                 Hscalefactor = 1.0,
                 extraFields = [],
                 initializeBase = True,
                 readFileToMemory = False,
                 refineNodes = 0):
                 
                 
        self.filename = filename
        self.nPerh = nNodePerh
        self.SPH = SPH
        self.extraFields = extraFields
        
        # For now we restrict to reading from a single (serial) file.
        allfiles = mpi.allreduce([filename], mpi.SUM)
        assert min([x == filename for x in allfiles])
        self.serialfile = True

        self.H = []
        self.fieldNames = []
        vals = []
        self.nf = 0
        self.nv = 0
        
        print "using " + str(mpi.procs) +  " procs"

        if mpi.rank == 0:
            f = open(filename,'r')
            self.f = f
            
            gotFieldNames = 0
            
            for line in self.f:
                line.rstrip('\r\n')
                data = line.split(delimiter)
                if data[0][0] != "#" and gotFieldNames == 1:
                    vals.append(data)
                if data[0][0] != "#" and gotFieldNames == 0:
                    self.fieldNames.append(data)
                    gotFieldNames = 1
            
            self.f.close()
                
            for i in xrange(len(self.fieldNames[0])):
                print self.fieldNames[0][i]
                self.nf = self.nf + 1
                
        
        self.nf = mpi.bcast(self.nf, root=0)
        self.nv = mpi.bcast(len(vals), root=0)
        self.fieldNames = mpi.bcast(self.fieldNames, root=0)
        mpi.barrier()

        for i in xrange(len(self.fieldNames[0])):
            self.__dict__[self.fieldNames[0][i]] = []

        assert self.nf == len(self.fieldNames[0])

        if mpi.rank == 0:
            n = len(vals)
            for i in xrange(n):
                for j in xrange(len(self.fieldNames[0])):
                    self.__dict__[self.fieldNames[0][j]].append(float(vals[i][j]))
                self.H.append((1.0/self.h[i]) * SymTensor3d.one)
        
        for j in xrange(len(self.fieldNames[0])):
            self.__dict__[self.fieldNames[0][j]] = mpi.bcast(self.__dict__[self.fieldNames[0][j]], root=0)
        
        self.H = mpi.bcast(self.H, root=0)
        
    
        # Initialize the base class.
        if initializeBase:
            fields = tuple([self.x, self.y, self.z, self.m, self.rho, self.vx, self.vy, self.vz, self.eps, self.H] +
                           [self.__dict__[x] for x in extraFields])
            NodeGeneratorBase.__init__(self, self.serialfile, *fields)


        # Apply the requested number of refinements.
        for i in xrange(refineNodes):
            refineNodes3d(self)

        return


    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localPosition(self, i):
        assert i >= 0 and i < len(self.x)
        assert len(self.x) == len(self.y)
        return Vector3d(self.x[i], self.y[i], self.z[i])
    
    #---------------------------------------------------------------------------
    # Get the velocity for the given node index.
    #---------------------------------------------------------------------------
    def localVelocity(self, i):
        assert i >= 0 and i < len(self.vx)
        assert len(self.vx) == len(self.vy)
        return Vector3d(self.vx[i], self.vy[i], self.vz[i])

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
        assert i >= 0 and i < len(self.rho)
        return self.rho[i]

    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        assert i >= 0 and i < len(self.H)
        return self.H[i]

