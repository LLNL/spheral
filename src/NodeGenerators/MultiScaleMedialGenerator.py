from math import *
import mpi

from NodeGeneratorBase import *
from MedialGenerator import *

#-------------------------------------------------------------------------------
# 
#-------------------------------------------------------------------------------
class MultiScaleMedialGenerator2d(NodeGeneratorBase):

    #---------------------------------------------------------------------------
    # Constructor.
    #---------------------------------------------------------------------------
    def __init__(self,
                 n,
                 rho,
                 boundary,
                 nstart = 100,
                 gradrho = None,
                 holes = [],
                 centroidFrac = 1.0,
                 maxIterationsPerStage = 100,
                 fracTol = 1.0e-3,
                 tessellationFileName = None,
                 nNodePerh = 2.01,
                 offset = (0.0, 0.0),
                 rejecter = None,
                 randomseed = 492739149274,
                 maxNodesPerDomain = 1000):

        # Kick things off by generating the coarsest level.
        ntarget = min(n, nstart)
        gen = MedialGeneratorBase(ndim = 2,
                                  n = ntarget,
                                  rho = rho,
                                  boundary = boundary,
                                  gradrho = gradrho,
                                  holes = holes,
                                  centroidFrac = centroidFrac,
                                  maxIterations = maxIterationsPerStage,
                                  fracTol = fracTol,
                                  tessellationFileName = tessellationFileName + "_gen000",
                                  nNodePerh = nNodePerh,
                                  randomseed = randomseed,
                                  maxNodesPerDomain = maxNodesPerDomain,
                                  seedPositions = None)

        # Iterate from coarse generators to fine until we hit the target number of
        # points.
        rangen = random.Random(randomseed)
        igeneration = 0
        while gen.globalNumNodes() != n:
            igeneration += 1

            # Split the generators we've created from the last stage.
            ntarget = min(n, 4*ntarget)
            seedPositions = []
            while mpi.allreduce(len(seedPositions), mpi.SUM) < ntarget:
                zeta = Vector2d(rangen.uniform(-0.5, 0.5), rangen.uniform(-0.5, 0.5))
                for i in xrange(gen.localNumNodes()):
                    hscale = sqrt(gen.vol[i]/pi)
                    xi = gen.pos[i] + hscale*zeta
                    if boundary.contains(xi, False):
                        use = True
                        ihole = 0
                        while use and ihole < len(holes):
                            use = not holes[ihole].contains(xi, False)
                            ihole += 1
                        if use:
                            seedPositions.append(xi)

            # Remove any excess points.
            nglobal = mpi.allreduce(len(seedPositions), mpi.SUM)
            while nglobal > ntarget:
                if mpi.rank < nglobal - ntarget and len(seedPositions) > 0:
                    seedPositions = seedPositions[:-1]
                nglobal = mpi.allreduce(len(seedPositions), mpi.SUM)
            assert nglobal == ntarget

            # Now let the MedialGenerator do its thing.
            gen = MedialGeneratorBase(ndim = 2,
                                      n = ntarget,
                                      seedPositions = seedPositions,
                                      rho = rho,
                                      boundary = boundary,
                                      gradrho = gradrho,
                                      holes = holes,
                                      centroidFrac = centroidFrac,
                                      maxIterations = maxIterationsPerStage,
                                      fracTol = fracTol,
                                      tessellationFileName = tessellationFileName + "_gen%03i" % igeneration,
                                      nNodePerh = nNodePerh,
                                      randomseed = randomseed,
                                      maxNodesPerDomain = maxNodesPerDomain)

        # Convert to our now regrettable standard coordinate storage for generators.
        self.x = [x.x + offset[0] for x in gen.pos]
        self.y = [x.y + offset[1] for x in gen.pos]

        # Copy the rest of the state.
        self.rhofunc = gen.rhofunc
        self.m = list(gen.m)
        self.H = list(gen.H)
        self.vol = list(gen.vol)
        self.surface = list(gen.surface)

        # If the user provided a "rejecter", give it a pass
        # at the nodes.
        if rejecter:
            self.x, self.y, self.m, self.H, self.vol, self.surface = rejecter(self.x,
                                                                              self.y,
                                                                              self.m,
                                                                              self.H,
                                                                              self.vol,
                                                                              self.surface)

        # Initialize the base class, taking into account the fact we've already broken
        # up the nodes between domains.
        NodeGeneratorBase.__init__(self, False,
                                   self.x, self.y, self.m, self.H, self.vol, self.surface)
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

        
