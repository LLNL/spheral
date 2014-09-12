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
        
        if mpi.rank == 0:
            f = open(filename,'r')
            self.f = f
    
            # create the field arrays
            vals = []
            self.H = []
            
            fieldNames = []
            gotFieldNames = 0
            
            for line in self.f:
                data = line.split(delimiter)
                if data[0][0] != "#" and gotFieldNames == 1:
                    vals.append(data)
                if data[0][0] != "#" and gotFieldNames == 0:
                    fieldNames.append(data)
                    gotFieldNames = 1
            
            self.f.close()
    
            print "in " + filename + " found " + str(len(fieldNames[0])) + " fields:"
            for i in xrange(len(fieldNames[0])):
                print fieldNames[0][i]
                self.__dict__[fieldNames[0][i]] = []
            
            
            n = len(vals)
            for i in xrange(n):
                for j in xrange(len(fieldNames[0])):
                    self.__dict__[fieldNames[0][j]].append(float(vals[i][j]))
                self.H.append((1.0/self.h[i]) * SymTensor2d.one)
        
        
        
        
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

        if mpi.rank == 0:
            f = open(filename,'r')
            self.f = f
    
            # create the field arrays
            vals = []
            self.H = []
            
            fieldNames = []
            gotFieldNames = 0
            
            for line in self.f:
                data = line.split(delimiter)
                if data[0][0] != "#" and gotFieldNames == 1:
                    vals.append(data)
                if data[0][0] != "#" and gotFieldNames == 0:
                    fieldNames.append(data)
                    gotFieldNames = 1
            
            self.f.close()
                
            print "in " + filename + " found " + str(len(fieldNames[0])) + " fields:"
            for i in xrange(len(fieldNames[0])):
                print fieldNames[0][i]
                self.__dict__[fieldNames[0][i]] = []


            n = len(vals)
            for i in xrange(n):
                for j in xrange(len(fieldNames[0])):
                    self.__dict__[fieldNames[0][j]].append(float(vals[i][j]))
                self.H.append((1.0/self.h[i]) * SymTensor3d.one)
    
        self.x = mpi.bcast(self.x, root=0)
        self.y = mpi.bcast(self.y, root=0)
        self.z = mpi.bcast(self.z, root=0)
        self.m = mpi.bcast(self.m, root=0)
        self.H = mpi.bcast(self.H, root=0)
        self.vx = mpi.bcast(self.vx, root=0)
        self.vy = mpi.bcast(self.vy, root=0)
        self.vz = mpi.bcast(self.vz, root=0)
        self.eps = mpi.bcast(self.eps, root=0)
        self.rho = mpi.bcast(self.rho, root=0)

    
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
