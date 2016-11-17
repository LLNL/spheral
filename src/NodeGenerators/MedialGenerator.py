from math import *
import mpi
import random

from NodeGeneratorBase import *
from Spheral import Vector2d, Tensor2d, SymTensor2d, \
     rotationMatrix2d, testPointInBox2d

from sobol import i4_sobol
from centroidalRelaxNodes import centroidalRelaxNodes

#-------------------------------------------------------------------------------
# 2D Generator.  Seeds positions within the given boundary using the Sobol 
# sequence to create the seeds.
#-------------------------------------------------------------------------------
class MedialGenerator2d(NodeGeneratorBase):

    #---------------------------------------------------------------------------
    # Constructor.
    #---------------------------------------------------------------------------
    def __init__(self,
                 n,
                 rho,
                 boundary,
                 gradrho = None,
                 holes = [],
                 maxIterations = 100,
                 fracTol = 1.0e-3,
                 tessellationFileName = None,
                 nNodePerh = 2.01,
                 offset = (0.0, 0.0),
                 rejecter = None,
                 randomseed = 492739149274):

        assert n > 0
        #assert len(holes) == 0   # Not supported yet, but we'll get there.

        # Load our handy 2D aliases.
        import Spheral2d as sph

        # Did we get passed a function or a constant for the density?
        if type(rho) == type(1.0):
            def rhofunc(posi):
                return rho
        else:
            rhofunc = rho
        self.rhofunc = rhofunc

        # Create a temporary NodeList we'll use store and update positions.
        eos = sph.GammaLawGasMKS(2.0, 2.0)
        WT = sph.TableKernel(sph.NBSplineKernel(7), 1000)
        nodes = sph.makeFluidNodeList("tmp generator nodes", 
                                      eos,
                                      hmin = 1e-10,
                                      hmax = 1e10,
                                      kernelExtent = WT.kernelExtent,
                                      hminratio = 1.0,
                                      nPerh = nNodePerh)

        # Make a first pass looking for the maximum density (roughly).
        pos = nodes.positions()
        mass = nodes.mass()
        rhof = nodes.massDensity()
        H = nodes.Hfield()
        nodes.numInternalNodes = n
        length = max(boundary.xmax.x - boundary.xmin.x,
                     boundary.xmax.y - boundary.xmin.y)
        rhomax = 0.0
        seed = 0
        i = 0
        while i < n:
            [coords, seed] = i4_sobol(2, seed)
            p = boundary.xmin + length*Vector2d(coords[0], coords[1])
            ihole = 0
            use = boundary.contains(p, False)
            if use:
                while use and ihole < len(holes):
                    use = not holes[ihole].contains(p, True)
                    ihole += 1
            if use:
                rhomax = max(rhomax, rhofunc(p))
                i += 1
        rhomax = mpi.allreduce(rhomax, mpi.MAX)
        print "MedialGenerator: selected a maximum density of ", rhomax

        # Initialize the desired number of generators in the boundary using the Sobol sequence.
        rangen = random.Random(randomseed)
        area = boundary.volume
        seed = 0
        i = 0
        while i < n:
            [coords, seed] = i4_sobol(2, seed)
            p = boundary.xmin + length*Vector2d(coords[0], coords[1])
            ihole = 0
            use = boundary.contains(p, False)
            if use:
                while use and ihole < len(holes):
                    use = not holes[ihole].contains(p, True)
                    ihole += 1
            if use:
                rhoi = rhofunc(p)
                if rangen.uniform(0.0, 1.0) < rhoi/rhomax:
                    pos[i] = p
                    rhof[i] = rhoi
                    mass[i] = rhoi * area/n  # Not actually correct, but mass will be updated in centroidalRelaxNodes
                    hi = 2.0 * nNodePerh * sqrt(area/n)
                    assert hi > 0.0
                    H[i] = SymTensor2d(1.0/hi, 0.0, 0.0, 1.0/hi)
                    i += 1

        # Add the holes to the boundary.
        points = sph.vector_of_Vector(boundary.vertices())
        facets = sph.vector_of_vector_of_unsigned(boundary.facetVertices)
        for hole in holes:
            ps = hole.vertices()
            fs = hole.facetVertices
            nold = points.size()
            nnew = nold + ps.size()
            for p in ps:
                points.append(p)
            for f in fs:
                assert len(f) == 2
                facets.append(sph.vector_of_unsigned(2))
                facets[-1][0] = nold + f[0]
                facets[-1][1] = nold + f[1]
                # facets[-1][0] = nold + f[1]
                # facets[-1][1] = nold + f[0]
        bound = sph.Polygon(points, facets)

        # Iterate the points toward centroidal relaxation.
        vol, surfacePoint = centroidalRelaxNodes([(nodes, bound)],
                                                 W = WT,
                                                 rho = rhofunc,
                                                 gradrho = gradrho,
                                                 maxIterations = maxIterations,
                                                 tessellationFileName = tessellationFileName)

        # Now we can fill out the usual Spheral generator info.
        self.x, self.y, self.m, self.H = [], [], [], []
        for i in xrange(n):
            self.x.append(pos[i].x)
            self.y.append(pos[i].y)
            self.m.append(vol(0,i) * rhofunc(pos[i]))
            hi = nNodePerh * sqrt(vol(0,i)/pi)
            assert hi > 0.0
            self.H.append(SymTensor2d(1.0/hi, 0.0, 0.0, 1.0/hi))
        assert len(self.x) == n
        assert len(self.y) == n
        assert len(self.m) == n
        assert len(self.H) == n

        # If the user provided a "rejecter", give it a pass
        # at the nodes.
        if rejecter:
            self.x, self.y, self.m, self.H = rejecter(self.x,
                                                      self.y,
                                                      self.m,
                                                      self.H)

        # Have the base class break up the serial node distribution
        # for parallel cases.
        NodeGeneratorBase.__init__(self, True,
                                   self.x, self.y, self.m, self.H)
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
        assert i >= 0 and i < len(self.H)
        return self.rhofunc(Vector2d(self.x[i], self.y[i]))
    
    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        assert i >= 0 and i < len(self.H)
        return self.H[i]
