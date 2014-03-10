#-------------------------------------------------------------------------------
# Convenience wrappers for 3D NodeGenerators that fill a faceted surface with
# points.
#-------------------------------------------------------------------------------
import mpi
from math import *
from numpy.polynomial import Polynomial as P
from NodeGeneratorBase import NodeGeneratorBase
from Spheral3d import Vector, Tensor, SymTensor, Polyhedron, \
    vector_of_Vector, vector_of_unsigned, vector_of_vector_of_unsigned, \
    rotationMatrix

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
def VFSurfacePolyhedron(filename,
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
    return surface

#-------------------------------------------------------------------------------
# Create a FacetedSurfaceGenerator based on a polyhedron in VF format in a file.
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
    surface = VFSurfacePolyhedron(filename, scaleFactor)
    return PolyhedralSurfaceGenerator(surface, rho, nx, ny, nz, seed, xmin, xmax,
                                      nNodePerh, SPH)

#-------------------------------------------------------------------------------
# Generate nodes inside a surface by extruding inward from the facets of a 
# surface.
#-------------------------------------------------------------------------------
class ExtrudedSurfaceGenerator(NodeGeneratorBase):

    def __init__(self,
                 surface,
                 lextrude,
                 nextrude,
                 dltarget,
                 dstarget,
                 rho,
                 flags = None,
                 nNodePerh = 2.01,
                 SPH = False):
        surfaceFacets = surface.facets()
        surfaceVertices = surface.vertices()
        
        # Find the facets surrounding each vertex.
        verts2facets = {}
        for fi in xrange(len(surfaceFacets)):
            f = surfaceFacets[fi]
            for ip in f.ipoints:
                if ip in verts2facets:
                    verts2facets[ip].append(fi)
                else:
                    verts2facets[ip] = [fi]
        
        # Find the effective normal for each vertex.
        vertnorms = {}
        for i in verts2facets:
            vertnorms[i] = Vector()
            for fi in verts2facets[i]:
                nhat = surfaceFacets[fi].normal
                wi = 1.0/surfaceFacets[fi].area
                vertnorms[i] += nhat * wi
            vertnorms[i] = vertnorms[i].unitVector()
        
        # Get the facets we need to extrude.
        if flags is None:
            facets = surfaceFacets
        else:
            assert len(flags) == len(surfaceFacets)
            facets = [surfaceFacets[i] for i in xrange(len(surfaceFacets)) if (flags[i] == 1)]
        
        # Find the maximum extent we need to use to cover the volumes of all extruded
        # facets.
        ymin, ymax, zmin, zmax = 1e200, -1e200, 1e200, -1e200
        for f in facets:
            nhat = f.normal
            T = rotationMatrix(nhat)
            p = f.position
            verts = [T*(surfaceVertices[i] - p) for i in f.ipoints]
            ymin = min(ymin, min([v.y for v in verts]))
            ymax = max(ymax, max([v.y for v in verts]))
            zmin = min(zmin, min([v.z for v in verts]))
            zmax = max(zmax, max([v.z for v in verts]))
        
        # Find the ratio needed for the spacing in the x direction.
        # We have to check for ratio=1 explicitly, since the series sum
        # doesn't work in that case.
        if abs(lextrude - nextrude*dltarget) < 1e-5*lextrude:
            ratio = 1.0
            l = lextrude
            dx = dltarget
        else:
            stuff = [0.0]*(nextrude + 1)
            stuff[0:1] = [dltarget - lextrude, -lextrude]
            stuff[-1] = dltarget
            p = P(stuff, [1.0e-3, 1.0e10], [1.0e-3, 1.0e10])
            complex_ratio = p.roots()[-1]
            assert complex_ratio.imag == 0.0
            ratio = complex_ratio.real
        
            # Unfortnately the root finding above isn't 100% accurate, so we 
            # adjust the initial step size to get the correct total length.
            l = dltarget*(1.0 - ratio**nextrude)/(1.0 - ratio)
            dx = dltarget * lextrude/l
        if mpi.rank == 0:
            print "FacetedSurfaceGenerator: selected ratio=%g, dxfirst=%g, l=%g." % (ratio, dx, l)
        
        # Build the template values we'll use to stamp into each facet volume.
        rt, mt, Ht = [], [], []
        xi = 0.0
        dxi = dx
        for ix in xrange(nextrude):
            dxi = dx*ratio**ix
            if ratio == 1.0:
                xi = -(ix + 0.5)*dx
            else:
                xi = -dx*(1.0 - ratio**ix)/(1.0 - ratio) + 0.5*dxi
            ds = min(4.0*dstarget, max(dstarget, dxi))
            ny = max(1, int((ymax - ymin)/ds + 0.5))
            nz = max(1, int((zmax - zmin)/ds + 0.5))
            dy = (ymax - ymin)/ny
            dz = (zmax - zmin)/nz
            for iy in xrange(ny):
                yi = ymin + (iy + 0.5)*dy
                for iz in xrange(nz):
                    zi = zmin + (iz + 0.5)*dz
                    rt.append(Vector(xi, yi, zi))
                    mt.append(rho*dxi*dy*dz)
                    Ht.append(SymTensor(1.0/(nNodePerh*dxi), 0.0, 0.0,
                                        0.0, 1.0/(nNodePerh*dy), 0.0,
                                        0.0, 0.0, 1.0/(nNodePerh*dz)))
        if mpi.rank == 0:
            print "FacetedSurfaceGenerator: built template block of %i points per facet." % len(rt)
        
        # Figure out a crude partitioning of the facets between processors.
        if len(facets) > mpi.procs:
            ndomain0 = len(facets)/mpi.procs
            remainder = len(facets) % mpi.procs
            assert remainder < mpi.procs
            ndomain = ndomain0
            if mpi.rank < remainder:
                ndomain += 1
            imin = mpi.rank*ndomain0 + min(mpi.rank, remainder)
            imax = imin + ndomain
        else:
            if mpi.rank < len(facets):
                imin = mpi.rank
                imax = imin + 1
            else:
                imin = 0
                imax = 0
                    
        # Now walk the facets and build our values.
        self.x, self.y, self.z, self.m, self.H = [], [], [], [], []
        for fi in xrange(imin, imax):
            if mpi.rank == 0:
                print "%i..." % fi
            f = facets[fi]
            nhat = f.normal
            p = f.position
            T = rotationMatrix(nhat)
            Ti = T.Transpose()
            verts = vector_of_Vector()
            for ip in f.ipoints:
                fphat = (p - surfaceVertices[ip]).unitVector()
                vi = surfaceVertices[ip] + fphat*0.05*dstarget
                verts.append(vi)
                verts.append(vi - vertnorms[ip]*lextrude)
            poly = Polyhedron(verts)   # Better be convex!
            for i in xrange(len(rt)):
                ri = Ti*rt[i] + p
                if poly.contains(ri):
                    self.x.append(ri.x)
                    self.y.append(ri.y)
                    self.z.append(ri.z)
                    self.m.append(mt[i])
                    self.H.append(SymTensor(Ht[i]))
                    self.H[-1].rotationalTransform(Ti)
        self.rho = [rho] * len(self.x)
        
        # Invoke the base class to finish up.
        NodeGeneratorBase.__init__(self, False, self.x, self.y, self.z, self.m, self.H, self.rho)

        return

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
