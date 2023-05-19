from math import *

from NodeGeneratorBase import *
#from Spheral import generateCylDistributionFromRZ

from Spheral import Vector2d, Vector3d
from Spheral import Tensor2d, Tensor3d
from Spheral import ScalarField2d, ScalarField3d
from Spheral import SymTensor2d, SymTensor3d, SymTensorField2d, SymTensorField3d
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
                 refineNodes = 0,
                 delimiter = ' ',
                 offset = None,
                 nodes = None):
        
        
        self.filename = filename
        self.nPerh = nNodePerh
        self.SPH = SPH
        self.extraFields = extraFields
        self.nodes = nodes
        
        # For now we restrict to reading from a single (serial) file.
        allfiles = mpi.allreduce([filename], mpi.SUM)
        assert min([x == filename for x in allfiles])
        self.serialfile = True
        
        self.H = []
        self.fieldNames = []
        vals = []
        self.nf = 0
        self.nv = 0
        
        print("using " + str(mpi.procs) +  " procs")
        
        if mpi.rank == 0:
            f = open(filename,'r')
            self.f = f
            
            gotFieldNames = 0
            
            for line in self.f:
                #line.rstrip('\r\n')
                #data = line.split(delimiter)
                data = (line.strip()).split(delimiter)
                if data[0][0] != "#" and gotFieldNames == 1:
                    vals.append(data)
                if data[0][0] != "#" and gotFieldNames == 0:
                    self.fieldNames.append(data)
                    gotFieldNames = 1
            
            self.f.close()
            
            for i in range(len(self.fieldNames[0])):
                print(self.fieldNames[0][i])
                self.nf = self.nf + 1
        
        
        self.nf = mpi.bcast(self.nf, root=0)
        self.nv = mpi.bcast(len(vals), root=0)
        self.fieldNames = mpi.bcast(self.fieldNames, root=0)
        mpi.barrier()
        
        for i in range(len(self.fieldNames[0])):
            self.__dict__[self.fieldNames[0][i]] = []
        
        assert self.nf == len(self.fieldNames[0])
        
        if mpi.rank == 0:
            n = len(vals)
            for i in range(n):
                for j in range(len(self.fieldNames[0])):
                    self.__dict__[self.fieldNames[0][j]].append(float(vals[i][j]))
                self.H.append((1.0/self.h[i]) * SymTensor2d.one)
        
        for j in range(len(self.fieldNames[0])):
            self.__dict__[self.fieldNames[0][j]] = mpi.bcast(self.__dict__[self.fieldNames[0][j]], root=0)
        
        self.H = mpi.bcast(self.H, root=0)

        if offset:
            for i in range(n):
                self.x[i] += offset[0]
                self.y[i] += offset[1]
        
        try:
            test = self.vx
        except AttributeError:
            print("Ascii file did not supply velocities (intentional?)")
            self.vx = []
            self.vy = []
            for i in range(len(self.x)):
                self.vx.append(0.0)
                self.vy.append(0.0)
    
        # Initialize the base class.
        if initializeBase:
            fields = tuple([self.x, self.y, self.m, self.rho, self.vx, self.vy, self.eps, self.H] +
                           [self.__dict__[x] for x in extraFields])
            NodeGeneratorBase.__init__(self, self.serialfile, *fields)
        
        # Apply the requested number of refinements.
        for i in range(refineNodes):
            refineNodes2d(self)
        
        # Finally, if the user provided a NodeList convert all our stored data to Fields
        # in order to have them follow the nodes around during redistribution.
        if nodes:
            n0 = len(self.x)
            nodes.numInternalNodes = n0
            for name in ['x', 'y', 'z', 'm', 'rho', 'vx', 'vy', 'vz', 'eps'] + extraFields:
                stuff = self.__dict__[name]
                assert len(stuff) == n0
                field = ScalarField2d("generator_" + name, nodes)
                for i in range(n0):
                    field[i] = stuff[i]
                self.__dict__[name] = field
            # H is a SymTensor, so do it separately.
            field = SymTensorField2d("generator_H", nodes)
            for i in range(n0):
                field[i] = self.H[i]
            self.H = field

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
                 refineNodes = 0,
                 rejecter=None,
                 delimiter = ' ',
                 offset=None,
                 nodes=None):
                 
                 
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
        
        print("using " + str(mpi.procs) +  " procs")

        if mpi.rank == 0:
            f = open(filename,'r')
            self.f = f
            
            gotFieldNames = 0
            
            for line in self.f:
                #line.rstrip('\r\n')
                #data = line.split(delimiter)
                data = (line.strip()).split(delimiter)
                if data[0][0] != "#" and gotFieldNames == 1:
                    vals.append(data)
                if data[0][0] != "#" and gotFieldNames == 0:
                    self.fieldNames.append(data)
                    gotFieldNames = 1
            
            self.f.close()
                
            for i in range(len(self.fieldNames[0])):
                print(self.fieldNames[0][i])
                self.nf = self.nf + 1
                
        
        self.nf = mpi.bcast(self.nf, root=0)
        self.nv = mpi.bcast(len(vals), root=0)
        self.fieldNames = mpi.bcast(self.fieldNames, root=0)
        mpi.barrier()

        for i in range(len(self.fieldNames[0])):
            self.__dict__[self.fieldNames[0][i]] = []

        assert self.nf == len(self.fieldNames[0])

        n = 0

        if mpi.rank == 0:
            n = len(vals)
            for i in range(n):
                for j in range(len(self.fieldNames[0])):
                    self.__dict__[self.fieldNames[0][j]].append(float(vals[i][j]))
                self.H.append((1.0/self.h[i]) * SymTensor3d.one)
        
        for j in range(len(self.fieldNames[0])):
            self.__dict__[self.fieldNames[0][j]] = mpi.bcast(self.__dict__[self.fieldNames[0][j]], root=0)
        
        self.H = mpi.bcast(self.H, root=0)
        
        n = mpi.bcast(n,root=0)
        
        if offset:
            for i in range(n):
                self.x[i] += offset[0]
                self.y[i] += offset[1]
                self.z[i] += offset[2]

        if rejecter:
            self.newH = []
                
            for i in range(len(self.fieldNames[0])):
                name = self.fieldNames[0][i] + 'new'
                self.__dict__[name] = []

            for i in range(n):
                #print rejecter.xmin, rejecter.xmax
                #print self.x
                if ((self.x[i] < rejecter.xmax and self.x[i] > rejecter.xmin) and (self.y[i] < rejecter.ymax and self.y[i] > rejecter.ymin) and (self.z[i] < rejecter.zmax and self.z[i] > rejecter.zmin)):
                    #print "condition green"
                    self.newH.append(self.H[i])
                    for j in range(len(self.fieldNames[0])):
                        name = self.fieldNames[0][j] + 'new'
                        self.__dict__[name].append(self.__dict__[self.fieldNames[0][j]][i])

            self.H = self.newH

            for j in range(len(self.fieldNames[0])):
                name = self.fieldNames[0][j] + 'new'
                self.__dict__[self.fieldNames[0][j]] = self.__dict__[name]

        try:
            test = self.vx
        except AttributeError:
            print("Ascii file did not supply velocities (intentional?)")
            self.vx = []
            self.vy = []
            self.vz = []
            for i in range(len(self.x)):
                self.vx.append(0.0)
                self.vy.append(0.0)
                self.vz.append(0.0)


        # Initialize the base class.
        if initializeBase:
            fields = tuple([self.x, self.y, self.z, self.m, self.rho, self.vx, self.vy, self.vz, self.eps, self.H] +
                           [self.__dict__[x] for x in extraFields])
            NodeGeneratorBase.__init__(self, self.serialfile, *fields)


        # Apply the requested number of refinements.
        for i in range(refineNodes):
            refineNodes3d(self)

        # Finally, if the user provided a NodeList convert all our stored data to Fields
        # in order to have them follow the nodes around during redistribution.
        if nodes:
            n0 = len(self.x)
            nodes.numInternalNodes = n0
            for name in ['x', 'y', 'z', 'm', 'rho', 'vx', 'vy', 'vz', 'eps'] + extraFields:
                stuff = self.__dict__[name]
                assert len(stuff) == n0
                field = ScalarField3d("generator_" + name, nodes)
                for i in range(n0):
                    field[i] = stuff[i]
                self.__dict__[name] = field
            # H is a SymTensor, so do it separately.
            field = SymTensorField3d("generator_H", nodes)
            for i in range(n0):
                field[i] = self.H[i]
            self.H = field

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

