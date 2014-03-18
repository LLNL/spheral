#-------------------------------------------------------------------------------
# Convenience wrappers for 3D NodeGenerators that fill a faceted surface with
# points.
#-------------------------------------------------------------------------------
import mpi
from NodeGeneratorBase import NodeGeneratorBase
from Spheral3d import Vector, Tensor, SymTensor, Polyhedron, \
    vector_of_Vector, vector_of_unsigned, vector_of_vector_of_unsigned

#-------------------------------------------------------------------------------
# General case where you hand in the surface polyhedron.
#-------------------------------------------------------------------------------
class PolyhedralSurfaceGenerator(NodeGeneratorBase):

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
        self.rho0 = rho

        # Figure out bounds and numbers of nodes to scan the volume with.
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

        # Some local geometry.
        ntot0 = nx*ny*nz
        dy = box.y/ny
        dz = box.z/nz
        volume = box.x * box.y * box.z
        self.m0 = rho*volume/ntot0
        hx = 1.0/(nNodePerh*dx)
        hy = 1.0/(nNodePerh*dy)
        hz = 1.0/(nNodePerh*dz)
        self.H0 = SymTensor(hx, 0.0, 0.0,
                            0.0, hy, 0.0,
                            0.0, 0.0, hz)

        # Determine the range of indices this domain will check.
        imin0, imax0 = self.globalIDRange(ntot0)

        # Check our local range of IDs and find what points are inside the specified surface.
        self.x, self.y, self.z = [], [], []
        for i in xrange(imin0, imax0):
            xi, yi, zi = self.hcpPosition(i, nx, ny, nz, dx, dy, dz, xmin, xmax)
            if self.surface.contains(Vector(xi, yi, zi)):
                self.x.append(xi)
                self.y.append(yi)
                self.z.append(zi)

        # At this point we have a less than optimal domain decomposition, but this will
        # be redistributed later anyway so take it and run.
        n = len(self.x)
        self.m = [self.m0]*n
        self.H = [self.H0]*n
        self.rho = [self.rho0]*n
        NodeGeneratorBase.__init__(self, False)
        return

    #-------------------------------------------------------------------------------
    # Seed positions/masses on an hcp lattice.
    #-------------------------------------------------------------------------------
    def hcpPosition(self, i, nx, ny, nz, dx, dy, dz, xmin, xmax):
        nxy = nx*ny
        ix = i % nx
        iy = (i / nx) % ny
        iz = i / nxy
        xx = xmin[0] + (ix + 0.5*((iy % 2) + (iz % 2)))*dx
        yy = xmin[1] + (iy + 0.5*(iz % 2))*dy
        zz = xmin[2] + (iz + 0.5)*dz
        return xx, yy, zz

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localPosition(self, i):
        assert i < len(self.x)
        return Vector(self.x[i], self.y[i], self.z[i])

    #---------------------------------------------------------------------------
    # Get the mass for the given node index.
    #---------------------------------------------------------------------------
    def localMass(self, i):
        assert i < len(self.m)
        return self.m[i]

    #---------------------------------------------------------------------------
    # Get the mass density for the given node index.
    #---------------------------------------------------------------------------
    def localMassDensity(self, i):
        assert i < len(self.rho)
        return self.rho[i]

    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        assert i < len(self.H)
        return self.H[i]

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
                       SPH = False,
                       scaleFactor = 1.0):
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
                verts.append(Vector(float(stuff[1]), float(stuff[2]), float(stuff[3]))*scaleFactor)
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

# #-------------------------------------------------------------------------------
# # Helper rejecter method given a polyhedral surface.
# #-------------------------------------------------------------------------------
# class PolyhedralSurfaceRejecter:

#     def __init__(self, surface):
#         self.surface = surface
#         return

#     def __call__(self, x, y, z, m, H):
#         n = len(x)
#         assert len(y) == n
#         assert len(z) == n
#         assert len(m) == n
#         assert len(H) == n

#         # We'll take advantage of any available parallelism to split
#         # up the containment testing.  The following algorithm is borrowed 
#         # from NodeGeneratorBase to divvy up the ID range.
#         ndomain0 = n/mpi.procs
#         remainder = n % mpi.procs
#         assert remainder < mpi.procs
#         ndomain = ndomain0
#         if mpi.rank < remainder:
#             ndomain += 1
#         imin = mpi.rank*ndomain0 + min(mpi.rank, remainder)
#         imax = imin + ndomain

#         # Check our local range of IDs.
#         xloc, yloc, zloc, mloc, Hloc = [], [], [], [], []
#         localIndices = [i for i in xrange(imin, imax)
#                            if self.surface.contains(Vector(x[i], y[i], z[i]))]

#         # Now cull to the interior values.
#         xnew, ynew, znew, mnew, Hnew = [], [], [], [], []
#         for iproc in xrange(mpi.procs):
#             otherIndices = mpi.bcast(localIndices, iproc)
#             for i in otherIndices:
#                 xnew.append(x[i])
#                 ynew.append(y[i])
#                 znew.append(z[i])
#                 mnew.append(m[i])
#                 Hnew.append(H[i])

#         # That's it.
#         return xnew, ynew, znew, mnew, Hnew
