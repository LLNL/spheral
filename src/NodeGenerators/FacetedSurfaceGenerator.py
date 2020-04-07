#-------------------------------------------------------------------------------
# Convenience wrappers for 3D NodeGenerators that fill a faceted surface with
# points.
#-------------------------------------------------------------------------------
import mpi
from math import *
import string, random
from numpy.polynomial import Polynomial as P
from NodeGeneratorBase import NodeGeneratorBase
from PolyhedronFileUtilities import readPolyhedronOBJ
from Spheral3d import Vector, Tensor, SymTensor, Polyhedron, \
    vector_of_Vector, vector_of_unsigned, vector_of_vector_of_unsigned, \
    rotationMatrix, refinePolyhedron, fillFacetedVolume2, \
    LinearOrder, TableKernel, BSplineKernel, GammaLawGasMKS, makeFluidNodeList
from centroidalRelaxNodes import *

#-------------------------------------------------------------------------------
# General case where you hand in the surface polyhedron.
#-------------------------------------------------------------------------------
class PolyhedralSurfaceGenerator(NodeGeneratorBase):

    def __init__(self,
                 surface,
                 rho,
                 resolution,
                 nNodePerh = 2.01,
                 SPH = False,
                 rejecter = None):
        self.surface = surface
        self.rho0 = rho

        # Figure out bounds and numbers of nodes to scan the volume with.
        xmin = surface.xmin
        xmax = surface.xmax
        box = xmax - xmin
        assert box.minElement > 0.0
        nx = max(1, int(box.x/resolution + 0.5))
        ny = max(1, int(box.y/resolution + 0.5))
        nz = max(1, int(box.z/resolution + 0.5))

        # Some local geometry.
        ntot0 = nx*ny*nz
        volume = box.x * box.y * box.z
        self.m0 = rho*volume/ntot0
        hx = 1.0/(nNodePerh*resolution)
        hy = 1.0/(nNodePerh*resolution)
        hz = 1.0/(nNodePerh*resolution)
        self.H0 = SymTensor(hx, 0.0, 0.0,
                            0.0, hy, 0.0,
                            0.0, 0.0, hz)

        # Build the intial positions.
        pos = fillFacetedVolume2(surface, resolution, mpi.rank, mpi.procs)
        nsurface = mpi.allreduce(len(pos), mpi.SUM)

        # Apply any rejecter.
        if rejecter:
            print "Applying rejection..."
            mask = [rejecter.accept(ri.x, ri.y, ri.z) for ri in pos]
        else:
            mask = [True]*len(pos)
        n0 = len(pos)
        self.x = [pos[i].x for i in xrange(n0) if mask[i]]
        self.y = [pos[i].y for i in xrange(n0) if mask[i]]
        self.z = [pos[i].z for i in xrange(n0) if mask[i]]
        n = len(self.x)

        # Pick a mass per point so we get exactly the correct total mass inside the surface
        # before any rejection.
        M0 = surface.volume * self.rho0
        self.m0 = M0/nsurface
        self.m = [self.m0]*n

        # At this point we have a less than optimal domain decomposition, but this will
        # be redistributed later anyway so take it and run.
        self.rho = [self.rho0]*n
        self.H = [self.H0]*n
        NodeGeneratorBase.__init__(self, False, self.m)
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
# This version includes centroidal optimization of the points.
#-------------------------------------------------------------------------------
class CentroidalPolyhedralSurfaceGenerator(PolyhedralSurfaceGenerator):

    def __init__(self,
                 surface,
                 rho,
                 resolution,
                 nNodePerh = 2.01,
                 SPH = False,
                 rejecter = None,
                 boundaries = [],
                 maxIterations = 100,
                 fracTol = 1.0e-3,
                 correctionOrder = LinearOrder,
                 centroidFrac = 0.25,
                 tessellationBaseDir = ".",
                 tessellationFileName = None):
        PolyhedralSurfaceGenerator.__init__(self,
                                            surface,
                                            rho,
                                            resolution,
                                            nNodePerh,
                                            SPH,
                                            rejecter)

        # Make temporary data structures so we can relax these points.
        W = TableKernel(BSplineKernel(), 1000)
        eos = GammaLawGasMKS(2.0, 2.0)
        n = self.localNumNodes()
        nodes = makeFluidNodeList("nodes" + "".join([random.choice(string.ascii_letters + string.digits) for nnn in xrange(32)]),
                                  eos,
                                  numInternal = n,
                                  kernelExtent = W.kernelExtent,
                                  xmin = surface.xmin,
                                  xmax = surface.xmax,
                                  nPerh = nNodePerh)
        pos = nodes.positions()
        H = nodes.Hfield()
        for i in xrange(n):
            pos[i] = self.localPosition(i)
            H[i] = self.localHtensor(i)
        nodes.neighbor().updateNodes()

        # Call the relaxer
        centroidalRelaxNodes(nodeListsAndBounds = [(nodes, surface)],
                             W = W,
                             rho = rho,
                             boundaries = boundaries,
                             maxIterations = maxIterations,
                             fracTol = fracTol,
                             correctionOrder = correctionOrder,
                             centroidFrac = centroidFrac,
                             tessellationBaseDir = tessellationBaseDir,
                             tessellationFileName = tessellationFileName)


#-------------------------------------------------------------------------------
# Create a FacetedSurfaceGenerator based on a polyhedron in VF format in a file.
#-------------------------------------------------------------------------------
def VFSurfaceGenerator(filename,
                       rho,
                       nx,
                       nNodePerh = 2.01,
                       SPH = False,
                       scaleFactor = 1.0,
                       refineFactor = 0,
                       rejecter = None):
    surface = None
    if mpi.rank == 0:
        surface = readPolyhedronOBJ(filename)
        if refineFactor != 0:
            surface = refinePolyhedron(surface, refineFactor)
        if scaleFactor != 1.0:
            surface *= scaleFactor
    surface = mpi.bcast(surface, 0)
    return PolyhedralSurfaceGenerator(surface, rho, nx, nNodePerh, SPH, rejecter)

#-------------------------------------------------------------------------------
# Generate nodes inside a surface by extruding inward from the facets of a 
# surface.
#-------------------------------------------------------------------------------
class ExtrudedSurfaceGenerator(NodeGeneratorBase):

    def __init__(self,
                 surface,
                 lconstant,
                 lextrude,
                 nextrude,
                 dltarget,
                 dstarget,
                 rho,
                 flags = None,
                 nNodePerh = 2.01,
                 SPH = False):
        self.surface = surface
        surfaceFacets = surface.facets
        surfaceVertices = surface.vertices
        vertexNorms = surface.vertexUnitNorms
        facetNeighbors = surface.facetFacetConnectivity

        assert lconstant <= lextrude

        # Figure out the extent and resolution of the ratioed region.
        nconstant = int(lconstant/dltarget + 0.5)
        lratioed = lextrude - lconstant
        nratioed = nextrude - nconstant

        # Get the facets we need to extrude.
        if flags is None:
            facets = surfaceFacets
            realFacetIDs = range(len(facets))
            flags = [1]*len(facets)
        else:
            assert len(flags) == len(surfaceFacets)
            realFacetIDs = [fi for fi in xrange(len(surfaceFacets)) if (flags[fi] == 1)]
            facets = [surfaceFacets[i] for i in realFacetIDs]
        
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
        if abs(lratioed - nratioed*dltarget) < 1e-5*lratioed:
            ratio = 1.0
            l = lratioed
            dx = dltarget
        else:
            stuff = [0.0]*(nratioed + 1)
            stuff[0:1] = [dltarget - lratioed, -lratioed]
            stuff[-1] = dltarget
            p = P(stuff, [1.0e-3, 1.0e10], [1.0e-3, 1.0e10])
            complex_ratio = p.roots()[-1]
            assert complex_ratio.imag == 0.0
            ratio = complex_ratio.real
        
            # Unfortunately the root finding above isn't 100% accurate, so we 
            # adjust the initial step size to get the correct total length.
            l = dltarget*(1.0 - ratio**nratioed)/(1.0 - ratio)
            dx = dltarget # * lratioed/l
        if mpi.rank == 0:
            print "FacetedSurfaceGenerator: selected ratio=%g, dxfirst=%g, l=%g." % (ratio, dx, l)
        
        # Build the template values we'll use to stamp into each facet volume.
        rt, mt, Ht = [], [], []
        dxi = dx
        xi = 0.5*dx
        ix = -1
        while xi > -lextrude:
            ix += 1
            if ix < nconstant:
                dxi = dx # lconstant/nconstant
            else:
                dxi = dx*ratio**(ix - nconstant)
            xi -= dxi
            ds = min(2.0*dstarget, max(dstarget, dxi))
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
                    
        # We need all the extruded facet polyhedra.  In general these may be intersecting!
        print "FacetedSurfaceGenerator: building extruded facets..."
        self.extrudedFacets, localExtrudedFacets = [], []
        for i in xrange(imin, imax):
            f = facets[i]
            p = f.position
            verts = vector_of_Vector()
            for ip in f.ipoints:
                fphat = (p - surfaceVertices[ip]).unitVector()
                vi = surfaceVertices[ip] + fphat*0.05*dstarget
                verts.append(vi)
                verts.append(vi - vertexNorms[ip]*lextrude)
            localExtrudedFacets.append(Polyhedron(verts))   # Better be convex!
        for i in xrange(mpi.procs):
            self.extrudedFacets.extend(mpi.bcast(localExtrudedFacets, i))
        assert len(self.extrudedFacets) == len(facets)

        # Now walk the facets and build our values.
        print "FacetedSurfaceGenerator: seeding points..."
        self.x, self.y, self.z, self.m, self.H, self.nodes2facets = [], [], [], [], [], []
        for i in xrange(imin, imax):
            f = facets[i]
            fi = realFacetIDs[i]
            neighbors = facetNeighbors[fi]
            p = f.position
            nhat = f.normal
            T = rotationMatrix(nhat)
            Ti = T.Transpose()
            for j in xrange(len(rt)):
                rj = Ti*rt[j] + p
                if self.extrudedFacets[i].contains(rj): # and surface.contains(rj):
                    stuff = [(surfaceFacets[k].distance(rj), k) for k in neighbors if flags[k] == 1]
                    stuff.sort()
                    if stuff[0][1] == fi:
                        self.x.append(rj.x)
                        self.y.append(rj.y)
                        self.z.append(rj.z)
                        self.m.append(mt[j])
                        self.H.append(SymTensor(Ht[j]))
                        self.H[-1].rotationalTransform(Ti)
                        self.nodes2facets.append(fi)
        self.rho = [rho] * len(self.x)
        
        # Invoke the base class to finish up.
        NodeGeneratorBase.__init__(self, False, self.x, self.y, self.z, self.m, self.H, self.rho, self.nodes2facets)

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
