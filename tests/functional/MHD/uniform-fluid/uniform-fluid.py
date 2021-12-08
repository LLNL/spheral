#-------------------------------------------------------------------------------
# The evolution of a uniform, magnetized conducting fluid.
#-------------------------------------------------------------------------------
from math import *
from Spheral import *
from SpheralTestUtilities import *
from SpheralVisitDump import dumpPhysicsState

# Load the mpi module if we"re parallel.
import loadmpi
mpi, procID, numProcs = loadmpi.loadmpi()

from GenerateNodeDistribution3d import *
from CubicNodeGenerator import GenerateCubicNodeDistribution

title("uniform magnetized conducting fluid test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(seed = "lattice",

            nx = 10,
            ny = 10,
            nz = 10,

            rho0 = 1.0,
            u0 = 0.0,
            B0 = Vector3d(1.0, 1.0, 1.0),
            nPerh = 1.25,
            mu0 = 1.0,

            gamma = 5.0/3.0,
            mu = 1.0,
            Qlimiter = True,
            balsaraCorrection = False,
            epsilon2 = 1e-2,
            negligibleSoundSpeed = 1e-5,
            csMultiplier = 1e-4,
            hmin = 1e-5,
            hmax = 1.0,
            hminratio = 0.05,
            HsmoothFraction = 0.0,
            cfl = 0.5,
            XSPH = True,
            epsilonTensile = 0.0,
            nTensile = 8,

            HEvolution = Hydro3d.HEvolutionType.IdealH,
            limitIdealH = False,

            neighborSearchType = Neighbor3d.NeighborSearchType.GatherScatter,
            numGridLevels = 20,
            topGridCellSize = 2.0,
            origin = Vector3d(0.0, 0.0, 0.0),

            goalTime = 0.1,
            dt = 0.0001,
            dtMin = 1.0e-5,
            dtMax = None,
            dtGrowth = 2.0,
            dtSample = 0.1,
            maxSteps = None,
            statsStep = 10,
            smoothIters = 0,
            sumForMassDensity = Hydro3d.MassDensityType.RigorousSumDensity,

            restoreCycle = None,

            L1v0 =   0.0889732,
            L1rho0 = 5.51975,
            L1eps0 = 0.04701,
            L1P0 =   1.66301,
            L1A0 =   0.00344783,

            graphics = False,
            )

#-------------------------------------------------------------------------------
# If we're using the cubic node generator, then scale things so we get a
# constant work per domain, and run to the same self-similar shock fraction
# of the node distribution.
#-------------------------------------------------------------------------------
if seed == "cubic":
    nxdomains = int(mpi.procs**(1.0/3.0) + 0.1)
    assert nxdomains**3 == mpi.procs
    nx *= nxdomains
    ny *= nxdomains
    nz *= nxdomains
    print nxdomains, nx, ny, nz

#-------------------------------------------------------------------------------
# A few derived variables.
#-------------------------------------------------------------------------------
xmin = (-0.5, -0.5, -0.5)
xmax = (0.5, 0.5, 0.5)

dataDir = "uniform-fluid-%ix%ix%i" % (nx, ny, nz)
visitDir = dataDir + "/visit"

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if not os.path.exists(visitDir):
        os.makedirs(visitDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS3d(gamma, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel3d(BSplineKernel3d(), 1000)
WTPi = TableKernel3d(BSplineKernel3d(), 1000)
output("WT")
output("WTPi")
kernelExtent = WT.kernelExtent()

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
# Perfectly conducting node list.
nodes = ConductingFluidNodeList("nodes", eos, WT, WTPi)
output("nodes")
nodes.HsmoothFraction = HsmoothFraction
nodes.XSPH = XSPH
nodes.nodesPerSmoothingScale = nPerh
nodes.epsilonTensile = epsilonTensile
nodes.nTensile = nTensile
nodes.hmin = hmin
nodes.hmax = hmax
nodes.hminratio = hminratio
output("nodes.HsmoothFraction")
output("nodes.nodesPerSmoothingScale")
output("nodes.epsilonTensile")
output("nodes.nTensile")
output("nodes.XSPH")
output("nodes.hmin")
output("nodes.hmax")
output("nodes.hminratio")

#-------------------------------------------------------------------------------
# Construct the neighbor object.
#-------------------------------------------------------------------------------
neighbor1 = NestedGridNeighbor3d(nodes,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes.registerNeighbor(neighbor1)

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    if seed == "cubic":
        from DistributeNodes import nullDistributeNodes3d
        generator = GenerateCubicNodeDistribution(nx, ny, nz, rho0,
                                                  xmin = xmin,
                                                  xmax = xmax,
                                                  nNodePerh = nPerh,
                                                  SPH = True)
        nullDistributeNodes3d((nodes, generator))
    else:
        from ParMETISDistributeNodes import distributeNodes3d
        generator = GenerateNodeDistribution3d(nx, ny, nz, rho0, seed,
                                               xmin = xmin,
                                               xmax = xmax,
                                               nNodePerh = nPerh,
                                               SPH = True)
        distributeNodes3d((nodes, generator))
    output("mpi.reduce(nodes.numInternalNodes, mpi.MIN)")
    output("mpi.reduce(nodes.numInternalNodes, mpi.MAX)")
    output("mpi.reduce(nodes.numInternalNodes, mpi.SUM)")

    # Set node specific thermal energies
    nodes.specificThermalEnergy(ScalarField3d("tmp", nodes, u0))

    # Set nodal magnetic inductions.
    for nodeID in xrange(nodes.numNodes):
        nodes.magneticInduction()[nodeID] = B0

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase3d()
output("db")
output("db.appendNodeList(nodes)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the artificial viscosities for the problem.
#-------------------------------------------------------------------------------
q = PriceMonaghanDissipation(1.0, 1.0, 1.0, 0.75, 1.0)
q.limiter = Qlimiter
q.balsaraShearCorrection = balsaraCorrection
q.epsilon2 = epsilon2
q.negligibleSoundSpeed = negligibleSoundSpeed
q.csMultiplier = csMultiplier
output("q")
output("q.Cl")
output("q.Cq")
output("q.limiter")
output("q.epsilon2")
output("q.negligibleSoundSpeed")
output("q.csMultiplier")
output("q.balsaraShearCorrection")

##-------------------------------------------------------------------------------
## Construct the hydro physics object.
##-------------------------------------------------------------------------------
hydro = Hydro3d(WT, WTPi, q)
hydro.cfl = cfl
hydro.HEvolution = HEvolution
hydro.sumForMassDensity = sumForMassDensity
hydro.HsmoothMin = hmin
hydro.HsmoothMax = hmax
#output("hydro")
#output("hydro.cfl")
#output("hydro.HEvolution")
#output("hydro.sumForMassDensity")
#output("hydro.HsmoothMin")
#output("hydro.HsmoothMax")
#output("hydro.kernel()")
#output("hydro.PiKernel()")
#output("hydro.valid()")

#-------------------------------------------------------------------------------
# Construct an MHD object.
#-------------------------------------------------------------------------------
MHD = MHD(WT, mu0)

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane1 = Plane3d(Vector3d(-0.5, 0.0, 0.0), Vector3d( 1.0, 0.0, 0.0))
xPlane2 = Plane3d(Vector3d( 0.5, 0.0, 0.0), Vector3d(-1.0, 0.0, 0.0))
yPlane1 = Plane3d(Vector3d( 0.0,-0.5, 0.0), Vector3d( 0.0, 1.0, 0.0))
yPlane2 = Plane3d(Vector3d( 0.0, 0.5, 0.0), Vector3d( 0.0,-1.0, 0.0))
zPlane1 = Plane3d(Vector3d( 0.0, 0.0,-0.5), Vector3d( 0.0, 0.0, 1.0))
zPlane2 = Plane3d(Vector3d( 0.0, 0.0, 0.5), Vector3d( 0.0, 0.0,-1.0))
xbc = PeriodicBoundary3d(xPlane1, xPlane2)
ybc = PeriodicBoundary3d(yPlane1, yPlane2)
zbc = PeriodicBoundary3d(zPlane1, zPlane2)
hydro.appendBoundary(xbc)
hydro.appendBoundary(ybc)
hydro.appendBoundary(zbc)
MHD.appendBoundary(xbc)
MHD.appendBoundary(ybc)
MHD.appendBoundary(zbc)

#-------------------------------------------------------------------------------
# Construct a time integrator.
#-------------------------------------------------------------------------------
integrator = SynchronousRK2Integrator3d(db)
integrator.appendPhysicsPackage(hydro)
integrator.appendPhysicsPackage(MHD)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.havePhysicsPackage(MHD)")
output("integrator.valid()")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
#raw_input()
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            initializeMassDensity = True)
output("control")

# Initial dump.
#dumpPhysicsState(integrator,
#                 "uniform-fluid-%ix%ix%i-visit" % (nx, ny, nz),
#                 visitDir)

#-------------------------------------------------------------------------------
# Advance one step.
#-------------------------------------------------------------------------------
control.step()
#control.advance(goalTime)
#dumpPhysicsState(integrator,
#                 "uniform-fluid-%ix%ix%i-visit" % (nx, ny, nz),
#                 visitDir)
v = nodes.velocity()
B = nodes.magneticInduction()
acc = nodes.DvelocityDt()

vMax = max([u.magnitude() for u in v.internalValues()])
vMin = min([u.magnitude() for u in v.internalValues()])
BMax = max([u.magnitude() for u in B.internalValues()])
BMin = min([u.magnitude() for u in B.internalValues()])
accMax = max([a.magnitude() for a in acc.internalValues()])
print 'max |v|: %g'%vMax
print 'min |v|: %g'%vMin
print 'max |B|: %g'%BMax
print 'min |B|: %g'%BMin
print 'Max |dv/dt|: %g'%accMax
