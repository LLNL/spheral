from math import *
import gzip

from NodeGeneratorBase import *
from Spheral import generateCylDistributionFromRZ

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
class GzipFileNodeGenerator2D(NodeGeneratorBase):

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

        self.precision = "%" + "%i.%ie" % (precision + 3, precision)

        # For now we restrict to reading from a single (serial) file.
        allfiles = mpi.allreduce([filename], mpi.SUM)
        assert min([x == filename for x in allfiles])
        self.serialfile = True

        # Open the file.
        if mpi.rank == 0:
            f = gzip.open(filename, "r")
            if readFileToMemory:
                self.f = f.readlines()
                f.close()
            else:
                self.f = f
        else:
            self.f = None

        # Read the positions.
        vals = readField2String(materialName, "positions", self.f)
        n = len(vals)
        self.x = []
        self.y = []
        for val in vals:
            x, y = tuple([float(x) for x in val.split()])
            self.x.append(x)
            self.y.append(y)
        assert len(self.x) == n
        assert len(self.y) == n

        # Read the masses.
        vals = readField2String(materialName, "mass", self.f)
        assert len(vals) == n
        self.m = [float(x) for x in vals]
        assert len(self.m) == n

        # Read the mass densities.
        vals = readField2String(materialName, "density", self.f)
        assert len(vals) == n
        self.rho = [float(x) for x in vals]
        assert len(self.rho) == n

        # Read the velocities.
        vals = readField2String(materialName, "velocity", self.f)
        assert len(vals) == n
        self.vx = []
        self.vy = []
        for val in vals:
            vx, vy = tuple([float(x) for x in val.split()])
            self.vx.append(vx)
            self.vy.append(vy)
        assert len(self.vx) == n
        assert len(self.vy) == n

        # Read the specific thermal energies.
        vals = readField2String(materialName, "specificThermalEnergy", self.f)
        assert len(vals) == n
        self.eps = [float(x) for x in vals]
        assert len(self.eps) == n

        # Read the H tensors.
        vals = readField2String(materialName, "Hinverse2", self.f)
        assert len(vals) == n
        self.H = []
        for val in vals:
            Hi2 = SymTensor2d(*tuple([float(x) for x in val.split()])) * Hscalefactor
            H = Hi2.Inverse().sqrt()
            if SPH:
                hxy = sqrt(H.Determinant())
                H = SymTensor2d.one * hxy
            self.H.append(H)
        assert len(self.H) == n

        # Read in any extra fields the user requested.
        # Note we make the assumption here that any extra fields are in fact scalar fields.
        for fname in extraFields:
            vals = readField2String(materialName, fname, self.f)
            assert len(vals) == n
            self.__dict__[fname] = [float(x) for x in vals]

        # Initialize the base class.
        if initializeBase:
            fields = tuple([self.x, self.y, self.m, self.rho, self.vx, self.vy, self.eps, self.H] +
                           [self.__dict__[x] for x in extraFields])
            NodeGeneratorBase.__init__(self, self.serialfile, *fields)

        if mpi.rank == 0:
            self.f.close()

        # Apply the requested number of refinements.
        for i in range(refineNodes):
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
class GzipFileNodeGenerator3D(NodeGeneratorBase):

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

        self.precision = "%" + "%i.%ie" % (precision + 3, precision)

        # For now we restrict to reading from a single (serial) file.
        allfiles = mpi.allreduce([filename], mpi.SUM)
        assert min([x == filename for x in allfiles])
        self.serialfile = True

        # Open the file.
        if mpi.rank == 0:
            f = gzip.open(filename, "r")
            if readFileToMemory:
                self.f = f.readlines()
                f.close()
            else:
                self.f = f
        else:
            self.f = None

        # Read the positions.
        vals = readField2String(materialName, "positions", self.f)
        n = len(vals)
        self.x = []
        self.y = []
        self.z = []
        for val in vals:
            x, y, z = tuple([float(x) for x in val.split()])
            self.x.append(x)
            self.y.append(y)
            self.z.append(z)
        assert len(self.x) == n
        assert len(self.y) == n
        assert len(self.z) == n

        # Read the masses.
        vals = readField2String(materialName, "mass", self.f)
        assert len(vals) == n
        self.m = [float(x) for x in vals]
        assert len(self.m) == n

        # Read the mass densities.
        vals = readField2String(materialName, "density", self.f)
        assert len(vals) == n
        self.rho = [float(x) for x in vals]
        assert len(self.rho) == n

        # Read the velocities.
        vals = readField2String(materialName, "velocity", self.f)
        assert len(vals) == n
        self.vx = []
        self.vy = []
        self.vz = []
        for val in vals:
            vx, vy, vz = tuple([float(x) for x in val.split()])
            self.vx.append(vx)
            self.vy.append(vy)
            self.vz.append(vz)
        assert len(self.vx) == n
        assert len(self.vy) == n
        assert len(self.vz) == n

        # Read the specific thermal energies.
        vals = readField2String(materialName, "specificThermalEnergy", self.f)
        assert len(vals) == n
        self.eps = [float(x) for x in vals]
        assert len(self.eps) == n

        # Read the H tensors.
        vals = readField2String(materialName, "Hinverse2", self.f)
        assert len(vals) == n
        self.H = []
        for val in vals:
            Hi2 = SymTensor3d(*tuple([float(x) for x in val.split()])) * Hscalefactor
            H = Hi2.Inverse().sqrt()
            if SPH:
                hxyz = (H.Determinant())**(1.0/3.0)
                H = SymTensor3d.one * hxyz
            self.H.append(H)
        assert len(self.H) == n

        # Read in any extra fields the user requested.
        # Note we make the assumption here that any extra fields are in fact scalar fields.
        for fname in extraFields:
            vals = readField2String(materialName, fname, self.f)
            assert len(vals) == n
            self.__dict__[fname] = [float(x) for x in vals]

        # Initialize the base class.
        if initializeBase:
            fields = tuple([self.x, self.y, self.z, self.m, self.rho, self.vx, self.vy, self.vz, self.eps, self.H] +
                           [self.__dict__[x] for x in extraFields])
            NodeGeneratorBase.__init__(self, self.serialfile, *fields)

        if mpi.rank == 0:
            self.f.close()

        # Apply the requested number of refinements.
        for i in range(refineNodes):
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

#-------------------------------------------------------------------------------
# Read in an RZ generated setup as 2-D.
#-------------------------------------------------------------------------------
class GzipFileNodeGeneratorRZto2D(GzipFileNodeGenerator2D):

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
                 refineNodes = 0):

        # Invoke the base class construction, but don't initialize the
        # NodeGeneratorBase yet.
        GzipFileNodeGenerator2D.__init__(self,
                                         filename,
                                         materialName,
                                         nNodePerh,
                                         SPH,
                                         precision,
                                         Hscalefactor,
                                         extraFields,
                                         False)

        # Flip the components of the positions and velocities, since
        # x->z, y->r.
        n = len(self.x)
        for i in range(n):
            self.x[i], self.y[i] = self.y[i], self.x[i]
            self.vx[i], self.vy[i] = self.vy[i], self.vx[i]

        # Scale the masses from RZ to 2-D planar geometry.
        for i in range(n):
            ri = self.y[i]
            c = 2.0*pi*ri
            self.m[i] *= c/(c*c + 1.0e-50)

        # Initialize the base class.
        if initializeBase:
            fields = tuple([self.x, self.y, self.m, self.rho, self.vx, self.vy, self.eps, self.H] +
                           [self.__dict__[x] for x in extraFields])
            NodeGeneratorBase.__init__(self, self.serialfile, *fields)

        # Apply the requested number of refinements.
        for i in range(refineNodes):
            refineNodes2d(self)

        return

#-------------------------------------------------------------------------------
# Read in an RZ generated setup as 3-D.
#-------------------------------------------------------------------------------
class GzipFileNodeGeneratorRZto3D(GzipFileNodeGeneratorRZto2D):

    #---------------------------------------------------------------------------
    # Constructor.
    #---------------------------------------------------------------------------
    def __init__(self,
                 filename,
                 materialName,
                 nNodePerh = 2.01,
                 phi = 2.0*pi,
                 SPH = False,
                 precision = 20,
                 Hscalefactor = 1.0,
                 extraFields = [],
                 refineNodes = 0):

        # Invoke the base class construction, but don't initialize the
        # NodeGeneratorBase yet.
        GzipFileNodeGeneratorRZto2D.__init__(self,
                                             filename,
                                             materialName,
                                             nNodePerh,
                                             SPH,
                                             precision,
                                             Hscalefactor,
                                             extraFields,
                                             False,
                                             refineNodes)

        # Convert the 2-D H tensors to 3-D, and correct the masses.
        n = len(self.x)
        H2dlist = self.H[:]
        assert len(H2dlist) == n
        self.H = []
        for i in range(n):
            xi = self.x[i]
            yi = self.y[i]
            H2d = H2dlist[i]

            # hxy0 = 0.5*(H2d.Inverse().Trace())
            # dphi = CylindricalBoundary.angularSpacing(yi, hxy0, nNodePerh, 2.0)
            # nhoop = int(2.0*pi/dphi + 0.5)
            # dphi = 2.0*pi/nhoop

            hmax = H2d.Inverse().eigenValues().maxElement()
            circ = phi*yi
            nhoop = max(1, int(circ/(hmax/nNodePerh) + 0.5))
            dphi = phi/nhoop

            hz = dphi*yi*nNodePerh
            self.H.append(SymTensor3d(H2d.xx, H2d.xy, 0.0,
                                      H2d.yx, H2d.yy, 0.0,
                                      0.0,    0.0,    1.0/hz))
            if SPH:
                h0 = self.H[-1].Determinant()**(1.0/3.0)
                self.H[-1] = SymTensor3d.one * h0

            circ = 2.0*pi*yi
            mhoop = self.m[i]*circ
            self.m[i] = mhoop

        assert len(self.m) == n
        assert len(self.H) == n

        # Duplicate the nodes from the xy-plane, creating rings of nodes about
        # the x-axis.  We use a C++ helper method for the sake of speed.
        kernelExtent = 2.0
        self.z = [0.0]*len(self.x)
        self.globalIDs = [0]*len(self.x)
        extras = [self.rho, self.vx, self.vy, self.eps]
        for fname in extraFields:
            extras.append(self.__dict__[fname])

        # Duplicate the nodes from the xy-plane, creating rings of nodes about
        # the x-axis.  We use a C++ helper method for the sake of speed.
        kernelExtent = 2.0
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
                                      rank, procs)
        self.x = [x for x in xvec]
        self.y = [x for x in yvec]
        self.z = [x for x in zvec]
        self.m = [x for x in mvec]
        self.H = [SymTensor3d(x) for x in Hvec]
        self.globalIDs = [x for x in globalIDsvec]
        for i in range(len(extras)):
            extras[i] = [x for x in extrasVec[i]]
        self.rho = extras[0][:]
        self.vx = extras[1][:]
        self.vy = extras[2][:]
        self.eps = extras[3][:]
        assert len(extraFields) == (len(extras) - 4)
        for i in range(len(extraFields)):
            self.__dict__[extraFields[i]] = extras[i + 4]

        # For consistency also generate the z velocities.
        n = len(self.x)
        self.vz = [0.0]*n

        # Post conditions.
        for fname in (["y", "z", "m", "rho", "vx", "vy", "vz", "eps", "H", "globalIDs"] + extraFields):
            if len(self.__dict__[fname]) != n:
                raise ValueError("GzipFileNodeGenerator3d: wrong field size for %s: %i != %i." % (fname,
                                                                                                   len(self.__dict__[fname]),
                                                                                                   n))
        
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

#-------------------------------------------------------------------------------
# Various conversion routines from strings.
#-------------------------------------------------------------------------------
def _string2int(x):
    return int(x)

def _string2bool(x):
    return x.replace(" ", "") == "True"

def _string2float(x):
    return float(x)

def _string2str(x):
    return x

def _string2Vector2d(x):
    return Vector2d(*tuple(x.split()))

def _string2Vector3d(x):
    return Vector3d(*tuple(x.split()))

def _string2Tensor2d(x):
    return Tensor2d(*tuple(x.split()))

def _string2Tensor3d(x):
    return Tensor3d(*tuple(x.split()))

def _string2SymTensor2d(x):
    return SymTensor2d(*tuple(x.split()))

def _string2SymTensor3d(x):
    return SymTensor3d(*tuple(x.split()))

_fromString = {int               : _string2int,
               bool              : _string2bool,
               float             : _string2float,
               str               : _string2str,
               type(Vector2d)    : _string2Vector2d,
               type(Vector3d)    : _string2Vector3d,
               type(Tensor2d)    : _string2Tensor2d,
               type(Tensor3d)    : _string2Tensor3d,
               type(SymTensor2d) : _string2Tensor2d,
               type(SymTensor3d) : _string2Tensor3d,
               }

#-------------------------------------------------------------------------------
# Read the given field values form the file.
# Returns the individual elements as strings.
#-------------------------------------------------------------------------------
def readField(field,
              materialName,
              fieldName,
              file):

    # First get the string values.
    strings = readField2string(materialName, fieldName, file)

    # Size the NodeList appropriately.
    n = len(strings)
    field.nodeList.numInternalNodes = n

    # Convert the field elements.
    for i in range(n):
        field[i] = _fromString[type(field[i])](strings[i])

    return

#-------------------------------------------------------------------------------
# Read the given field values form the file.
# Returns the individual elements as strings.
#-------------------------------------------------------------------------------
def readField2String(materialName,
                     fieldName,
                     file):
    delimiter = "$"
    result = ""
    if mpi.rank == 0:
        found = False
        for line in file:
            vals = line.split(delimiter)
            if vals[0][0] != "#":
                if (vals[0] == materialName and
                    vals[1] == fieldName):
                    n = int(vals[2])
                    result = vals[3:-1]
                    assert len(result) == n
                    found = True
                    break
        if not found:
            raise ValueError("Unable to find %s %s" % (materialName, fieldName))

    result = mpi.bcast(result, 0)
    return result

