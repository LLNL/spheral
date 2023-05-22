from math import *
import mpi

from NodeGeneratorBase import *

from Spheral import *

#-------------------------------------------------------------------------------
# Read a polymesh from a silo file, turning each zone into an SPH point.
# Assumed 3D here.
#-------------------------------------------------------------------------------
class SiloPolyMeshGenerator(NodeGeneratorBase):

    #---------------------------------------------------------------------------
    # Ordinary constructor.
    #---------------------------------------------------------------------------
    def __init__(self, 
                 fileName,               # Name of the file
                 meshName,               # Name of mesh variable in file
                 rho0,                   # Initial mass density
                 serialFile = True,      # Should this be treated as a serial file or broken up in parallel?
                 nNodePerh = 2.01,       # number of nodes per smoothing scale
                 SPH = False,            # Force round H tensors
                 scale = 1.0):           # Optionally scale the coordinates by some factor

        self.x, self.y, self.z, self.m, self.H = [], [], [], [], []
        self.rho0 = rho0
        if rank == 0 or (not serialFile):

            # Read the file to a set of positions, volumes, and H's.
            pos = vector_of_Vector3d()
            vol = vector_of_double()
            H = vector_of_SymTensor3d()
            readSiloPolyMesh(fileName, meshName, pos, vol, H)
            print("Read %i points from %s." % (len(pos), fileName))
            assert len(pos) == len(vol) == len(H)

            self.x = [scale * x.x for x in pos]
            self.y = [scale * x.y for x in pos]
            self.z = [scale * x.z for x in pos]
            self.m = [scale**3 * x*rho0 for x in vol]
            self.H = [x/scale for x in H]

        self.x = mpi.bcast(self.x, root=0)
        self.y = mpi.bcast(self.y, root=0)
        self.z = mpi.bcast(self.z, root=0)
        self.m = mpi.bcast(self.m, root=0)
        self.H = mpi.bcast(self.H, root=0)

        # Initialize the base class, which will break up the serial node distribution
        # for parallel cases if required.
        NodeGeneratorBase.__init__(self, serialFile,
                                   self.x, self.y, self.z, self.m, self.H)

        # If SPH has been specified, make sure the H tensors are round.
        if SPH:
            self.makeHround()

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
        assert i >= 0 and i < len(self.x)
        return self.rho0

    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        assert i >= 0 and i < len(self.H)
        return self.H[i]

