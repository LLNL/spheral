from math import *

from NodeGeneratorBase import *

from Spheral import Vector3d as Vector
from Spheral import Tensor3d as Tensor
from Spheral import SymTensor3d as SymTensor

import mpi
procID = mpi.rank
nProcs = mpi.procs

#-------------------------------------------------------------------------------
# Node generator to read ASCI Abaqus formatted files.
# Assumed 3D here.
#-------------------------------------------------------------------------------
class AbaqusNodeGenerator(NodeGeneratorBase):

    def __init__(self, 
                 fileName,               # Name of Abaqus file.
                 materialLabel,          # Name of the material we're generating from in the file
                 elsetLabel,             # Name of the element set we're generating from in the file
                 serialFile = True,      # Should this be treated as a serial file or broken up in parallel?
                 nNodePerh = 2.01,       # number of nodes per smoothing scale
                 SPH = False):           # Force round H tensors

        if rank == 0 or (not serialFile):

            # Read the whole darn file.
            f = open(fileName, "r")
            lines = f.readlines()
            f.close()

            # Find and read the vertex coordinates.
            vertices = {}
            iline = 0
            while (iline < len(lines) and
                   lines[iline][:5] != "*NODE"):
                iline += 1
            if iline == len(lines):
                raise RuntimeError, "Unable to find *NODE specification in %s" % fileName
            iline += 1
            while lines[iline][0] != "*":
                stuff = lines[iline].split(",")
                assert len(stuff) == 4
                i = int(stuff[0])
                xi, yi, zi = (float(stuff[j]) for j in xrange(1,4))
                vertices[i] = Vector(xi, yi, zi)
                iline += 1
            print "AbaqusNodeGenerator : Read %i vertices from file %s" % (len(vertices), fileName)

            # Read the element types for each element, and build our info based on that.
            self.x, self.y, self.z, self.m, self.H = [], [], [], [], []
            elset = "ELSET=%s" % elsetLabel
            iline = 0
            while (iline < len(lines) and
                   lines[iline][:5] != "*ELEMENT" and 
                   not elset in lines[iline].replace(" ", "")):
                iline += 1
            if iline == len(lines):
                raise RuntimeError, "Unable to find *ELEMENT section for %s in %s" % (elset, fileName)
            iline += 1
            while lines[iline][0] != "*":
                stuff = lines[iline].split(",")
                i = int(stuff[0])
                indices = [int(x) for x in stuff[1:]]
                assert len(indices) >= 4
                vol, pos, psi = self.computeCellMoments([vertices[i] for i in indices])
                self.m.append(vol)                     # Will rescale with density once we read it in.
                self.x.append(pos.x)
                self.y.append(pos.y)
                self.z.append(pos.z)
                assert psi.Determinant() > 0.0
                psi = psi.sqrt()
                psi *= (vol/psi.Determinant())**(1.0/3.0) * nNodePerh
                psi = psi.Inverse()
                self.H.append(psi)
                iline += 1
            print "AbaqusNodeGenerator : Read %i cells for element set %s from file %s" % (len(self.x), elsetLabel, fileName)

        # Find the mass density and turn the masses into actual masses rather than volumes.
        material = "*MATERIAL,NAME=%s" % materialLabel
        iline = 0
        while (iline < len(lines) and 
               not material in lines[iline].replace(" ", "")):
            iline += 1
        if iline == len(lines):
            raise RuntimeError, "Unable to find *MATERIAL section for %s in %s" % (materialLabel, fileName)
        while (iline < len(lines) and
               lines[iline][:8] != "*DENSITY"):
            iline += 1
        if iline == len(lines):
            raise RuntimeError, "Unable to find *DENSITY section for %s in %s" % (materialLabel, fileName)
        iline += 1
        assert iline < len(lines)
        self.rho0 = float(lines[iline])
        self.m = [x*self.rho0 for x in self.m]

        # Initialize the base class, which will break up the serial node distribution
        # for parallel cases if required.
        NodeGeneratorBase.__init__(self, serialFile,
                                   self.x, self.y, self.z, self.m, self.H)

        # If SPH has been specified, make sure the H tensors are round.
        if SPH:
            self.makeHround()

        return

    #-------------------------------------------------------------------------------
    # Compute the zeroth (volume), first, and second moments of the given point
    # distribution.
    #
    # Assumed node orderings:
    # Tetrahedron:                  3
    #                              /|\
    #                             / | \
    #                            /  |  \
    #                           /   |   \
    #                          /    2    \
    #                         /   *   *   \
    #                        /  *       *  \
    #                       / *           * \
    #                      0-----------------1
    # 
    # Hexahedron:                7--------------6
    #                           /*             /|
    #                          / *   f6       / |    face 1: (0,1,2,3)
    #                         4--------------5  |    face 2: (0,3,7,4)
    #                         |f2*    f5     |f4|    face 3: (0,4,5,1)
    #                         |  *   f3      |  |    face 4: (1,5,6,2)
    #                         |  3***********|**2    face 5: (2,6,7,3)
    #                         | *    f1      | /     face 6: (4,7,6,5)
    #                         0--------------1
    #
    #-------------------------------------------------------------------------------
    def computeCellMoments(self, vertices):
        assert len(vertices) in (4, 8)

        # Compute the zeroth and first moments based on decomposing into tets.
        if len(vertices) == 4:
            # Tetrahedron.
            vol, cent = self.tetMoments(vertices)
        elif len(vertices) == 8:
            # Hexahedron.  Facet the faces.
            f1 = 0.25*(vertices[0] + vertices[1] + vertices[2] + vertices[3])
            f2 = 0.25*(vertices[0] + vertices[3] + vertices[7] + vertices[4])
            f3 = 0.25*(vertices[0] + vertices[4] + vertices[5] + vertices[1])
            f4 = 0.25*(vertices[1] + vertices[5] + vertices[6] + vertices[2])
            f5 = 0.25*(vertices[2] + vertices[6] + vertices[7] + vertices[3])
            f6 = 0.25*(vertices[4] + vertices[7] + vertices[6] + vertices[5])
            cell = (f1 + f2 + f3 + f4 + f5 + f6)/6.0
            vol, cent = 0.0, Vector()
            for tetSet in ((vertices[0], vertices[1], f1, cell),
                           (vertices[1], vertices[2], f1, cell),
                           (vertices[2], vertices[3], f1, cell),
                           (vertices[3], vertices[0], f1, cell),
                           (vertices[4], vertices[0], f2, cell),
                           (vertices[0], vertices[3], f2, cell),
                           (vertices[3], vertices[7], f2, cell),
                           (vertices[7], vertices[4], f2, cell),
                           (vertices[0], vertices[4], f3, cell),
                           (vertices[4], vertices[5], f3, cell),
                           (vertices[5], vertices[1], f3, cell),
                           (vertices[1], vertices[0], f3, cell),
                           (vertices[1], vertices[5], f4, cell),
                           (vertices[5], vertices[6], f4, cell),
                           (vertices[6], vertices[2], f4, cell),
                           (vertices[2], vertices[1], f4, cell),
                           (vertices[3], vertices[2], f5, cell),
                           (vertices[2], vertices[6], f5, cell),
                           (vertices[6], vertices[7], f5, cell),
                           (vertices[7], vertices[3], f5, cell),
                           (vertices[5], vertices[4], f6, cell),
                           (vertices[4], vertices[7], f6, cell),
                           (vertices[7], vertices[6], f6, cell),
                           (vertices[6], vertices[5], f6, cell)):
                vol_tet, cent_tet = self.tetMoments(tetSet)
                vol += vol_tet
                cent += vol_tet * cent_tet
            assert vol > 0.0
            cent /= vol

        # The second moment is just the dyadic sum of the vertices.
        psi = SymTensor()
        for vi in vertices:
            psi += (vi - cent).selfdyad()
        psi /= len(vertices)

        return vol, cent, psi

    #-------------------------------------------------------------------------------
    # Compute the volume and centroid of a tet.
    # Assumes ordering shown above.
    #-------------------------------------------------------------------------------
    def tetMoments(self, vertices):
        assert len(vertices) == 4
        vol = ((vertices[1] - vertices[0]).cross(vertices[2] - vertices[0])).dot(vertices[3] - vertices[0])
        cent = 0.25*(vertices[0] + vertices[1] + vertices[2] + vertices[3])
        assert vol > 0.0
        return vol, cent
            
    #-------------------------------------------------------------------------------
    # Get the position for the given node index.
    #-------------------------------------------------------------------------------
    def localPosition(self, i):
        assert i >= 0 and i < len(self.x)
        assert len(self.x) == len(self.y)
        return Vector(self.x[i], self.y[i], self.z[i])

    #-------------------------------------------------------------------------------
    # Get the mass for the given node index.
    #-------------------------------------------------------------------------------
    def localMass(self, i):
        assert i >= 0 and i < len(self.m)
        return self.m[i]

    #-------------------------------------------------------------------------------
    # Get the mass density for the given node index.
    #-------------------------------------------------------------------------------
    def localMassDensity(self, i):
        assert i >= 0 and i < len(self.x)
        return self.rho0

    #-------------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #-------------------------------------------------------------------------------
    def localHtensor(self, i):
        assert i >= 0 and i < len(self.H)
        return self.H[i]

