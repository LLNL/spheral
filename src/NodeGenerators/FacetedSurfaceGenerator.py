#-------------------------------------------------------------------------------
# Convenience wrappers for 3D NodeGenerators that fill a faceted surface with
# points.
#-------------------------------------------------------------------------------
import mpi
from GenerateNodeDistribution3d import GenerateNodeDistribution3d
from Spheral3d import Vector, Tensor, SymTensor, Polyhedron, \
    vector_of_Vector, vector_of_unsigned, vector_of_vector_of_unsigned

#-------------------------------------------------------------------------------
# General case where you hand in the surface polyhedron.
#-------------------------------------------------------------------------------
class PolyhedralSurfaceGenerator(GenerateNodeDistribution3d):

    def __init__(self,
                 surface,
                 rho,
                 nx,
                 ny = None,
                 nz = None,
                 seed = "hcp",
                 xmin = None,
                 xmax = None,
                 nNodePerh = 2.01,
                 SPH = False):
        self.surface = surface
        if xmin is None:
            xmin = surface.xmin
        if xmax is None:
            xmax = surface.xmax
        box = xmax - xmin
        assert box.minElement > 0.0
        dx = box.x/nx
        if ny is None:
            ny = int(box.y/dx + 0.5)
        if nz is None:
            nz = int(box.z/dx + 0.5)
        GenerateNodeDistribution3d.__init__(self,
                                            n1 = nx,
                                            n2 = ny,
                                            n3 = nz,
                                            rho = rho,
                                            distributionType = seed,
                                            xmin = xmin,
                                            xmax = xmax,
                                            nNodePerh = nNodePerh,
                                            SPH = SPH,
                                            rejecter = PolyhedralSurfaceRejecter(surface))
        return

#-------------------------------------------------------------------------------
# Read a VF (vertex-facet labeled) shape file to create the polyhedral surface.
# This is the format used by the NASA radar asteroid shape models:
#  http://echo.jpl.nasa.gov/asteroids/shapes/shapes.html
#-------------------------------------------------------------------------------
def VFSurfaceGenerator(filename,
                       rho,
                       nx,
                       ny = None,
                       nz = None,
                       seed = "hcp",
                       xmin = None,
                       xmax = None,
                       nNodePerh = 2.01,
                       SPH = False):
    surface = None
    if mpi.rank == 0:
        f = open(filename, "r")
        verts = vector_of_Vector()
        facets = vector_of_vector_of_unsigned()
        for line in f:
            stuff = line.split()
            assert stuff[0] in ("v", "f")
            if stuff[0] == "v":
                assert len(stuff) == 4
                verts.append(Vector(float(stuff[1]), float(stuff[2]), float(stuff[3])))
            else:
                assert len(stuff) >= 4
                facets.append(vector_of_unsigned())
                for x in stuff[1:]:
                    facets[-1].append(int(x) - 1)
        f.close()
        nverts = len(verts)
        for i in xrange(len(facets)):
            for j in xrange(len(facets[i])):
                assert facets[i][j] < nverts
        surface = Polyhedron(verts, facets)
    surface = mpi.bcast(surface)
    return PolyhedralSurfaceGenerator(surface, rho, nx, ny, nz, seed, xmin, xmax,
                                      nNodePerh, SPH)

#-------------------------------------------------------------------------------
# Helper rejecter method given a polyhedral surface.
#-------------------------------------------------------------------------------
class PolyhedralSurfaceRejecter:

    def __init__(self, surface):
        self.surface = surface
        return

    def __call__(self, x, y, z, m, H):
        n = len(x)
        assert len(y) == n
        assert len(z) == n
        assert len(m) == n
        assert len(H) == n
        xnew, ynew, znew, mnew, Hnew = [], [], [], [], []
        for i in xrange(n):
            j = 0
            if self.surface.contains(Vector(x[i], y[i], z[i])):
                xnew.append(x[i])
                ynew.append(y[i])
                znew.append(z[i])
                mnew.append(m[i])
                Hnew.append(H[i])
        return xnew, ynew, znew, mnew, Hnew
