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
        #if readFileToMemory:
        #self.f = f.readlines()
        #f.close()
        # i don't want to build in this functionality right now
        #else:
        # self.f = f
        else:
            self.f = None
        
        # create the field arrays
        vals = []
        
        self.x = []
        self.y = []
        self.m = []
        self.rho = []
        self.eps = []
        self.vx = []
        self.vy = []
        self.vz = []
        self.H = []
        
        for line in f:
            data = line.split()
            vals.append(data)
        
        n = len(vals) - 1
        for i in range(n):
            j = i + 1
            self.x.append(float(vals[j][0]))
            self.y.append(float(vals[j][1]))
            # (1/h) * [[1 0 ][0 1]]
            H = (1.0/float(vals[j][2])) * SymTensor2d.one
            self.H.append(H)
            
            self.m.append(float(vals[j][3]))
            self.rho.append(float(vals[j][4]))
            #pressure
            self.eps.append(float(vals[j][6]))
            self.vx.append(float(vals[j][7]))
            self.vy.append(float(vals[j][8]))
            #temperature
            #abund array
        
        
        
        
        # Initialize the base class.
        if initializeBase:
            fields = tuple([self.x, self.y, self.m, self.rho, self.vx, self.vy, self.eps, self.H] +
                           [self.__dict__[x] for x in extraFields])
            NodeGeneratorBase.__init__(self, self.serialfile, *fields)
        
        if mpi.rank == 0:
            self.f.close()
        
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
            #if readFileToMemory:
                #self.f = f.readlines()
                #f.close()
                # i don't want to build in this functionality right now
                #else:
                # self.f = f
        else:
            self.f = None
            
        # create the field arrays
        vals = []
        
        self.x = []
        self.y = []
        self.z = []
        self.m = []
        self.rho = []
        self.eps = []
        self.vx = []
        self.vy = []
        self.vz = []
        self.H = []
        
        for line in f:
            data = line.split()
            vals.append(data)
            
        n = len(vals) - 1
        for i in range(n):
            j = i + 1
            self.x.append(float(vals[j][0]))
            self.y.append(float(vals[j][1]))
            self.z.append(float(vals[j][2]))
            # (1/h) * [[1 0 0][0 1 0][0 0 1]]
            H = (1.0/float(vals[j][3])) * SymTensor3d.one
            self.H.append(H)
            
            self.m.append(float(vals[j][4]))
            self.rho.append(float(vals[j][5]))
            #pressure
            self.eps.append(float(vals[j][7]))
            self.vx.append(float(vals[j][8]))
            self.vy.append(float(vals[j][9]))
            self.vz.append(float(vals[j][10]))
            #temperature
            #abund array
            
        # Read in extra fields
        # maybe later...

        
        # Initialize the base class.
        if initializeBase:
            fields = tuple([self.x, self.y, self.z, self.m, self.rho, self.vx, self.vy, self.vz, self.eps, self.H] +
                           [self.__dict__[x] for x in extraFields])
            NodeGeneratorBase.__init__(self, self.serialfile, *fields)

        if mpi.rank == 0:
            self.f.close()

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
