from math import *
import mpi

from NodeGeneratorBase import *
from MedialGenerator import *

#-------------------------------------------------------------------------------
# Dimension agnostic base version.
#-------------------------------------------------------------------------------
class MultiScaleMedialGeneratorBase(NodeGeneratorBase):

    #---------------------------------------------------------------------------
    # Constructor.
    #---------------------------------------------------------------------------
    def __init__(self,
                 ndim,
                 n,
                 rho,
                 boundary,
                 nstart,
                 gradrho,
                 holes,
                 centroidFrac,
                 maxIterationsPerStage,
                 fracTol,
                 tessellationFileName,
                 nNodePerh,
                 offset,
                 rejecter,
                 randomseed,
                 maxNodesPerDomain,
                 enforceConstantMassPoints,
                 cacheFileName):

        # Load our handy aliases.
        if ndim == 2:
            import Spheral2d as sph
        else:
            import Spheral3d as sph

        # Kick things off by generating the coarsest level.
        ntarget = min(n, nstart)
        tfname = tessellationFileName
        if tfname:
            tfname += "_gen000"
        gen = MedialGeneratorBase(ndim = ndim,
                                  n = ntarget,
                                  rho = rho,
                                  boundary = boundary,
                                  gradrho = gradrho,
                                  holes = holes,
                                  centroidFrac = centroidFrac,
                                  maxIterations = maxIterationsPerStage,
                                  fracTol = fracTol,
                                  tessellationFileName = tfname,
                                  nNodePerh = nNodePerh,
                                  randomseed = randomseed,
                                  maxNodesPerDomain = maxNodesPerDomain,
                                  enforceConstantMassPoints = enforceConstantMassPoints,
                                  seedPositions = None,
                                  cacheFileName = None)

        # If there is an pre-existing cache file, load it instead of doing all the work.
        if not gen.restoreState(cacheFileName):

            # Iterate from coarse generators to fine until we hit the target number of
            # points.
            rangen = random.Random(randomseed)
            igeneration = 0
            while gen.globalNumNodes() != n:
                igeneration += 1
        
                # Split the generators we've created from the last stage.
                ntarget = min(n, ntarget * 2**ndim)
                seedPositions = []
                while mpi.allreduce(len(seedPositions), mpi.SUM) < ntarget:
                    zeta = sph.Vector(rangen.uniform(-0.5, 0.5), rangen.uniform(-0.5, 0.5), rangen.uniform(-0.5, 0.5))
                    for i in range(gen.localNumNodes()):
                        hscale = (gen.vol[i]/pi)**(1.0/ndim)
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
                tfname = tessellationFileName
                if tfname:
                    tfname += "_gen%03i" % igeneration
                gen = MedialGeneratorBase(ndim = ndim,
                                          n = ntarget,
                                          seedPositions = seedPositions,
                                          rho = rho,
                                          boundary = boundary,
                                          gradrho = gradrho,
                                          holes = holes,
                                          centroidFrac = centroidFrac,
                                          maxIterations = max(10, maxIterationsPerStage//(2**(ndim*igeneration))),
                                          fracTol = fracTol,
                                          tessellationFileName = tfname,
                                          nNodePerh = nNodePerh,
                                          randomseed = randomseed,
                                          maxNodesPerDomain = maxNodesPerDomain,
                                          enforceConstantMassPoints = enforceConstantMassPoints,
                                          cacheFileName = None)
        
            # If requested, we can store the state of the generator such that it can be
            # later restored without going through all that work.
            if cacheFileName:
                gen.dumpState(cacheFileName)

        # Convert to our now regrettable standard coordinate storage for generators.
        self.x = [x.x + offset[0] for x in gen.pos]
        self.y = [x.y + offset[1] for x in gen.pos]
        if ndim == 3:
            self.z = [x.z + offset[2] for x in gen.pos]

        # Copy the rest of the state.
        self.ndim = ndim
        self.rhofunc = gen.rhofunc
        self.m = list(gen.m)
        self.H = list(gen.H)
        self.vol = list(gen.vol)
        self.surface = list(gen.surface)

        # Some dimension trickery.
        if ndim == 2:
            args = (self.x, self.y, self.m, self.H, self.vol, self.surface)
        else:
            args = (self.x, self.y, self.z, self.m, self.H, self.vol, self.surface)

        # If the user provided a "rejecter", give it a pass
        # at the nodes.
        if rejecter:
            newargs = rejecter(*args)
            if ndim == 2:
                self.x, self.y, self.m, self.H, self.vol, self.surface = newargs
                args = (self.x, self.y, self.m, self.H, self.vol, self.surface)
            else:
                self.x, self.y, self.z, self.m, self.H, self.vol, self.surface = newargs
                args = (self.x, self.y, self.z, self.m, self.H, self.vol, self.surface)

        # Initialize the base class, taking into account the fact we've already broken
        # up the nodes between domains.
        NodeGeneratorBase.__init__(self, False, *args)
        return

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localPosition(self, i):
        assert i >= 0 and i < len(self.x)
        assert len(self.x) == len(self.y)
        if self.ndim == 2:
            return Vector2d(self.x[i], self.y[i])
        else:
            assert len(self.x) == len(self.z)
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
        assert i >= 0 and i < len(self.H)
        return self.rhofunc(self.localPosition(i))
    
    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        assert i >= 0 and i < len(self.H)
        return self.H[i]
        
#-------------------------------------------------------------------------------
# 2D
#-------------------------------------------------------------------------------
class MultiScaleMedialGenerator2d(MultiScaleMedialGeneratorBase):

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
                 maxIterationsPerStage = 1000,
                 fracTol = 1.0e-3,
                 tessellationFileName = None,
                 nNodePerh = 2.01,
                 offset = (0.0, 0.0),
                 rejecter = None,
                 randomseed = 492739149274,
                 maxNodesPerDomain = 1000,
                 enforceConstantMassPoints = False,
                 cacheFileName = None):

        MultiScaleMedialGeneratorBase.__init__(self,
                                               ndim = 2,
                                               n = n,
                                               rho = rho,
                                               boundary = boundary,
                                               nstart = nstart,
                                               gradrho = gradrho,
                                               holes = holes,
                                               centroidFrac = centroidFrac,
                                               maxIterationsPerStage = maxIterationsPerStage,
                                               fracTol = fracTol,
                                               tessellationFileName = tessellationFileName,
                                               nNodePerh = nNodePerh,
                                               offset = offset,
                                               rejecter = rejecter,
                                               randomseed = randomseed,
                                               maxNodesPerDomain = maxNodesPerDomain,
                                               enforceConstantMassPoints = enforceConstantMassPoints,
                                               cacheFileName = cacheFileName)

        return

#-------------------------------------------------------------------------------
# 3D
#-------------------------------------------------------------------------------
class MultiScaleMedialGenerator3d(MultiScaleMedialGeneratorBase):

    #---------------------------------------------------------------------------
    # Constructor.
    #---------------------------------------------------------------------------
    def __init__(self,
                 n,
                 rho,
                 boundary,
                 nstart = 2000,
                 gradrho = None,
                 holes = [],
                 centroidFrac = 1.0,
                 maxIterationsPerStage = 1000,
                 fracTol = 1.0e-3,
                 tessellationFileName = None,
                 nNodePerh = 2.01,
                 offset = (0.0, 0.0, 0.0),
                 rejecter = None,
                 randomseed = 492739149274,
                 maxNodesPerDomain = 1000,
                 enforceConstantMassPoints = False,
                 cacheFileName = None):

        MultiScaleMedialGeneratorBase.__init__(self,
                                               ndim = 3,
                                               n = n,
                                               rho = rho,
                                               boundary = boundary,
                                               nstart = nstart,
                                               gradrho = gradrho,
                                               holes = holes,
                                               centroidFrac = centroidFrac,
                                               maxIterationsPerStage = maxIterationsPerStage,
                                               fracTol = fracTol,
                                               tessellationFileName = tessellationFileName,
                                               nNodePerh = nNodePerh,
                                               offset = offset,
                                               rejecter = rejecter,
                                               randomseed = randomseed,
                                               maxNodesPerDomain = maxNodesPerDomain,
                                               enforceConstantMassPoints = enforceConstantMassPoints,
                                               cacheFileName = cacheFileName)

        return

