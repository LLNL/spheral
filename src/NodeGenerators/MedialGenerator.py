from math import *
import mpi
import random

from NodeGeneratorBase import *
from Spheral import Vector2d, Tensor2d, SymTensor2d, \
     rotationMatrix2d, testPointInBox2d

from sobol import i4_sobol
from centroidalRelaxNodes import centroidalRelaxNodes

#-------------------------------------------------------------------------------
# Base MedialGenerator.  This implements the generic 2D/3D algorithm.
#-------------------------------------------------------------------------------
class MedialGeneratorBase(NodeGeneratorBase):

    #---------------------------------------------------------------------------
    # Constructor.
    #---------------------------------------------------------------------------
    def __init__(self,
                 ndim,
                 n,
                 rho,
                 boundary,
                 gradrho,
                 holes,
                 maxIterations,
                 fracTol,
                 tessellationFileName,
                 nNodePerh,
                 randomseed,
                 maxNodesPerDomain):

        assert ndim in (2,3)
        assert n > 0

        # Load our handy 2D aliases.
        if ndim == 2:
            import Spheral2d as sph
        else:
            import Spheral3d as sph

        # Did we get passed a function or a constant for the density?
        if type(rho) in (float, int):
            def rhofunc(posi):
                return rho
            rhomax = rho
        else:
            rhofunc = rho
            rhomax = None
        self.rhofunc = rhofunc

        # Some useful geometry.
        box = boundary.xmax - boundary.xmin
        length = box.maxElement()
        boundvol = boundary.volume
        for hole in holes:
            boundvol -= hole.volume
        boxvol = 1.0
        for idim in xrange(ndim):
            boxvol *= box[idim]
        fracOccupied = min(1.0, boxvol/boundvol)
        assert fracOccupied > 0.0 and fracOccupied <= 1.0

        # Create a temporary NodeList we'll use store and update positions.
        eos = sph.GammaLawGasMKS(2.0, 2.0)
        WT = sph.TableKernel(sph.NBSplineKernel(7), 1000)
        hmax = 2.0*length
        nodes = sph.makeFluidNodeList("tmp generator nodes", 
                                      eos,
                                      hmin = 1e-10,
                                      hmax = hmax,
                                      kernelExtent = WT.kernelExtent,
                                      hminratio = 1.0,
                                      nPerh = nNodePerh,
                                      topGridCellSize = 2.0*WT.kernelExtent*hmax)

        # Make a first pass looking for the maximum density (roughly).
        pos = nodes.positions()
        mass = nodes.mass()
        rhof = nodes.massDensity()
        H = nodes.Hfield()
        imin, imax = self.globalIDRange(n)
        nlocal = imax - imin
        nodes.numInternalNodes = nlocal

        # If necessary probe for a maximum density statistically.
        rangen = random.Random(randomseed + mpi.rank)
        if not rhomax:
            rhomax = 0.0
            nglobal = 0
            while nglobal < n:
                p = boundary.xmin + length*sph.Vector(rangen.random(), rangen.random(), rangen.random())
                use = boundary.contains(p, False)
                if use:
                    ihole = 0
                    while use and ihole < len(holes):
                        use = not holes[ihole].contains(p, True)
                        ihole += 1
                if use:
                    rhomax = max(rhomax, rhofunc(p))
                    i = 1
                else:
                    i = 0
                nglobal += mpi.allreduce(i, mpi.SUM)
            rhomax = mpi.allreduce(rhomax, mpi.MAX)
        print "MedialGenerator: selected a maximum density of ", rhomax

        # It's a bit tricky to properly use the Sobol sequence in parallel.  We handle this by searching for the lowest
        # seeds that give us the desired number of points.
        seeds = []
        seed = 0
        while mpi.allreduce(len(seeds)) < n:
            localseed = seed + mpi.rank
            [coords, newseed] = i4_sobol(ndim, localseed)
            p = boundary.xmin + length*sph.Vector(*tuple(coords))
            use = boundary.contains(p, False)
            if use:
                ihole = 0
                while use and ihole < len(holes):
                    use = not holes[ihole].contains(p, True)
                    ihole += 1
            if use:
                rhoi = rhofunc(p)
                if rangen.random() < rhoi/rhomax:
                    seeds.append(localseed)
            seed += mpi.procs

        # Drop the highest value seeds to ensure we have the correct number of total points.
        nglobal = mpi.allreduce(len(seeds), mpi.SUM)
        assert n + mpi.procs >= nglobal
        seeds.sort()
        seeds = [-1] + seeds
        while mpi.allreduce(len(seeds)) > n + mpi.procs:
            maxseed = mpi.allreduce(seeds[-1], mpi.MAX)
            assert maxseed > -1
            if seeds[-1] == maxseed:
                seeds = seeds[:-1]
        seeds = seeds[1:]

        # Load balance the number of seeds per domain.
        if len(seeds) > nlocal:
            extraseeds = seeds[nlocal:]
        else:
            extraseeds = []
        extraseeds = mpi.allreduce(extraseeds, mpi.SUM)
        seeds = seeds[:nlocal]
        for iproc in xrange(mpi.procs):
            ngrab = max(0, nlocal - len(seeds))
            ntaken = mpi.bcast(ngrab, root=iproc)
            if mpi.rank == iproc:
                seeds += extraseeds[:ngrab]
            extraseeds = extraseeds[ntaken:]
        assert len(extraseeds) == 0
        assert len(seeds) == nlocal
        assert mpi.allreduce(len(seeds)) == n

        # Initialize the desired number of generators in the boundary using the Sobol sequence.
        for i, seed in enumerate(seeds):
            [coords, newseed] = i4_sobol(ndim, seed)
            p = boundary.xmin + length*sph.Vector(*tuple(coords))
            rhoi = rhofunc(p)
            pos[i] = p
            rhof[i] = rhoi
            mass[i] = rhoi * boundvol/n  # Not actually correct, but mass will be updated in centroidalRelaxNodes
            hi = min(hmax, 2.0 * nNodePerh * (boundvol/n)**(1.0/ndim))
            assert hi > 0.0
            H[i] = sph.SymTensor.one / hi

        # Each domain has independently generated the correct number of points, but they are randomly distributed.
        # Before going further it's useful to try and spatially collect the points by domain.
        # We'll use the Spheral Peano-Hilbert space filling curve implementation to do this.
        if mpi.procs > 1:
            db = sph.DataBase()
            db.appendNodeList(nodes)
            maxNodes = max(maxNodesPerDomain, 2*n/mpi.procs)
            redistributor = sph.PeanoHilbertOrderRedistributeNodes(0.0,
                                                                   minNodesPerDomainFraction = 0.0,
                                                                   maxNodesPerDomainFraction = float(maxNodes)/n)
            redistributor.redistributeNodes(db)
            boundaries = [sph.NestedGridDistributedBoundary.instance()]
        else:
            boundaries = []

        # Iterate the points toward centroidal relaxation.
        vol, surfacePoint = centroidalRelaxNodes([(nodes, boundary, holes)],
                                                 W = WT,
                                                 rho = rhofunc,
                                                 gradrho = gradrho,
                                                 boundaries = boundaries,
                                                 maxIterations = maxIterations,
                                                 tessellationFileName = tessellationFileName)

        # Store the values the descendent generators will need.
        self.vol, self.surface, self.pos, self.m, self.H = [], [], [], [], []
        for i in xrange(nodes.numInternalNodes):
            self.vol.append(vol(0,i))
            self.surface.append(surfacePoint(0,i))
            self.pos.append(sph.Vector(pos[i]))
            self.m.append(vol(0,i) * rhofunc(pos[i]))
            self.H.append(sph.SymTensor(H[i]))
        assert mpi.allreduce(len(self.vol), mpi.SUM) == n
        assert mpi.allreduce(len(self.surface), mpi.SUM) == n
        assert mpi.allreduce(len(self.pos), mpi.SUM) == n
        assert mpi.allreduce(len(self.m), mpi.SUM) == n
        assert mpi.allreduce(len(self.H), mpi.SUM) == n

        return

#-------------------------------------------------------------------------------
# 2D Generator.  Seeds positions within the given boundary using the Sobol 
# sequence to create the seeds.
# Note as currently implemented this code does not necessarily generate 
# identical output for different numbers of processors!
#-------------------------------------------------------------------------------
class MedialGenerator2d(MedialGeneratorBase):

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
                 randomseed = 492739149274,
                 maxNodesPerDomain = 1000):

        # The base generator does most of the work.
        MedialGeneratorBase.__init__(self,
                                     ndim = 2,
                                     n = n,
                                     rho = rho,
                                     boundary = boundary,
                                     gradrho = gradrho,
                                     holes = holes,
                                     maxIterations = maxIterations,
                                     fracTol = fracTol,
                                     tessellationFileName = tessellationFileName,
                                     nNodePerh = nNodePerh,
                                     randomseed = randomseed,
                                     maxNodesPerDomain = maxNodesPerDomain)


        # Convert to our now regrettable standard coordinate storage for generators.
        self.x = [x.x + offset[0] for x in self.pos]
        self.y = [x.y + offset[1] for x in self.pos]
        del self.pos

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

#-------------------------------------------------------------------------------
# 3D Generator.  Seeds positions within the given boundary using the Sobol 
# sequence to create the seeds.
#-------------------------------------------------------------------------------
class MedialGenerator3d(NodeGeneratorBase):

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

        # Load our handy 3D aliases.
        import Spheral3d as sph

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
            [coords, seed] = i4_sobol(3, seed)
            p = boundary.xmin + length*Vector3d(coords[0], coords[1], coords[2])
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
        volbound = boundary.volume
        seed = 0
        i = 0
        while i < n:
            [coords, seed] = i4_sobol(3, seed)
            p = boundary.xmin + length*Vector3d(coords[0], coords[1], coords[2])
            use = boundary.contains(p, False)
            if use:
                ihole = 0
                while use and ihole < len(holes):
                    use = not holes[ihole].contains(p, True)
                    ihole += 1
            if use:
                rhoi = rhofunc(p)
                if rangen.uniform(0.0, 1.0) < rhoi/rhomax:
                    pos[i] = p
                    rhof[i] = rhoi
                    mass[i] = rhoi * volbound/n  # Not actually correct, but mass will be updated in centroidalRelaxNodes
                    hi = 2.0 * nNodePerh * (volbound/(n*4.0/3.0*pi))**(1.0/3.0)
                    assert hi > 0.0
                    H[i] = SymTensor3d(1.0/hi, 0.0, 0.0,
                                       0.0, 1.0/hi, 0.0,
                                       0.0, 0.0, 1.0/hi)
                    i += 1

        # Iterate the points toward centroidal relaxation.
        vol, surfacePoint = centroidalRelaxNodes([(nodes, boundary, holes)],
                                                 W = WT,
                                                 rho = rhofunc,
                                                 gradrho = gradrho,
                                                 maxIterations = maxIterations,
                                                 tessellationFileName = tessellationFileName)

        # Now we can fill out the usual Spheral generator info.
        self.x, self.y, self.z, self.m, self.H = [], [], [], [], []
        for i in xrange(n):
            self.x.append(pos[i].x)
            self.y.append(pos[i].y)
            self.z.append(pos[i].z)
            self.m.append(vol(0,i) * rhofunc(pos[i]))
            hi = nNodePerh * (vol(0,i)/(4.0/3.0*pi))**(1.0/3.0)
            assert hi > 0.0
            self.H.append(SymTensor3d(1.0/hi, 0.0, 0.0,
                                      0.0, 1.0/hi, 0.0,
                                      0.0, 0.0, 1.0/hi))
        assert len(self.x) == n
        assert len(self.y) == n
        assert len(self.z) == n
        assert len(self.m) == n
        assert len(self.H) == n

        # If the user provided a "rejecter", give it a pass
        # at the nodes.
        if rejecter:
            self.x, self.y, self.z, self.m, self.H = rejecter(self.x,
                                                              self.y,
                                                              self.z,
                                                              self.m,
                                                              self.H)

        # Have the base class break up the serial node distribution
        # for parallel cases.
        NodeGeneratorBase.__init__(self, True,
                                   self.x, self.y, self.z, self.m, self.H)
        return

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localPosition(self, i):
        assert i >= 0 and i < len(self.x)
        assert len(self.x) == len(self.y)
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
        return self.rhofunc(Vector3d(self.x[i], self.y[i], self.z[i]))
    
    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        assert i >= 0 and i < len(self.H)
        return self.H[i]
