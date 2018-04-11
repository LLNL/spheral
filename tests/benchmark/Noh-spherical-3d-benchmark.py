#-------------------------------------------------------------------------------
# The spherical Noh test case.
#
# W.F. Noh 1987, JCP, 72, 78-120.
#-------------------------------------------------------------------------------
from math import *
from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from findLastRestart import *

# Load the mpi module if we"re parallel.
import loadmpi
mpi, procID, numProcs = loadmpi.loadmpi()

from GenerateNodeDistribution3d import *

title("3-D integrated hydro test -- spherical Noh problem")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(NodeListConstructor = SphNodeList3d,

            seed = "lattice",

            nx = 40,
            ny = 40,
            nz = 40,

            rho1 = 1.0,
            eps1 = 0.0,
            vr1 = -1.0,
            nPerh = 2.01,

            rmin = 0.0,
            rmax = 1.0,

            gamma = 5.0/3.0,
            mu = 1.0,
            #Qconstructor = MonaghanGingoldViscosity3d,
            Qconstructor = TensorMonaghanGingoldViscosity3d,
            IntegratorConstructor = PredictorCorrectorIntegrator3d,
            Cl = 1.0,
            Cq = 0.75,
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
            compatibleEnergy = True,
            gradhCorrection = False,

            HEvolution = Hydro3d.HEvolutionType.IdealH,
            limitIdealH = False,

            goalTime = 0.6,
            dt = 0.0001,
            dtMin = 1.0e-5,
            dtMax = None,
            dtGrowth = 2.0,
            dtSample = 0.1,
            maxSteps = 10,
            statsStep = 10,
            smoothIters = 0,
            sumForMassDensity = Hydro3d.MassDensityType.RigorousSumDensity,

            restoreCycle = None,
            restartStep = 20,

            graphics = "gnu",
            )

xmin = (0.0, 0.0, 0.0)
xmax = (1.0, 1.0, 1.0)

dataDir = "Noh-spherical-3d-%ix%ix%i" % (nx, ny, nz)
restartDir = dataDir + "/restarts"
visitDir = dataDir + "/visit"
restartBaseName = restartDir + "/Noh-spherical-3d-%ix%ix%i" % (nx, ny, nz)

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
    if not os.path.exists(visitDir):
        os.makedirs(visitDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# If we're restarting, find the set of most recent restart files.
#-------------------------------------------------------------------------------
restoreCycle = None # findLastRestart(restartBaseName)

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
nodes1 = NodeListConstructor("nodes", eos, WT, WTPi)
output("nodes1")
nodes1.HsmoothFraction = HsmoothFraction
nodes1.XSPH = XSPH
nodes1.nodesPerSmoothingScale = nPerh
nodes1.epsilonTensile = epsilonTensile
nodes1.nTensile = nTensile
nodes1.hmin = hmin
nodes1.hmax = hmax
nodes1.hminratio = hminratio
output("nodes1.HsmoothFraction")
output("nodes1.nodesPerSmoothingScale")
output("nodes1.epsilonTensile")
output("nodes1.nTensile")
output("nodes1.XSPH")
output("nodes1.hmin")
output("nodes1.hmax")
output("nodes1.hminratio")

#-------------------------------------------------------------------------------
# Construct the neighbor object.
#-------------------------------------------------------------------------------
neighbor1 = TreeNeighbor3d(nodes1,
                           kernelExtent = kernelExtent)
nodes1.registerNeighbor(neighbor1)

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    if numProcs > 1:
        from ParMETISDistributeNodes import distributeNodes3d
    else:
        from DistributeNodes import distributeNodes3d
    generator = GenerateNodeDistribution3d(nx, ny, nz, rho1, seed,
                                           xmin = xmin,
                                           xmax = xmax,
                                           rmin = rmin,
                                           rmax = rmax,
                                           nNodePerh = nPerh)
    distributeNodes3d((nodes1, generator))
    output("mpi.reduce(nodes1.numInternalNodes, mpi.MIN)")
    output("mpi.reduce(nodes1.numInternalNodes, mpi.MAX)")
    output("mpi.reduce(nodes1.numInternalNodes, mpi.SUM)")

    # Set node specific thermal energies
    nodes1.specificThermalEnergy(ScalarField3d("tmp", nodes1, eps1))

    # Set node velocities
    for nodeID in xrange(nodes1.numNodes):
        nodes1.velocity()[nodeID] = nodes1.positions()[nodeID].unitVector()*vr1

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase3d()
output("db")
output("db.appendNodeList(nodes1)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the artificial viscosities for the problem.
#-------------------------------------------------------------------------------
q = Qconstructor(Cl, Cq)
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

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
hydro = Hydro3d(WT, WTPi, q, compatibleEnergy, gradhCorrection)
hydro.cfl = cfl
hydro.HEvolution = HEvolution
hydro.sumForMassDensity = sumForMassDensity
hydro.HsmoothMin = hmin
hydro.HsmoothMax = hmax
output("hydro")
output("hydro.cfl")
output("hydro.HEvolution")
output("hydro.sumForMassDensity")
output("hydro.HsmoothMin")
output("hydro.HsmoothMax")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.valid()")

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane0 = Plane3d(Vector3d(0.0, 0.0, 0.0), Vector3d(1.0, 0.0, 0.0))
yPlane0 = Plane3d(Vector3d(0.0, 0.0, 0.0), Vector3d(0.0, 1.0, 0.0))
zPlane0 = Plane3d(Vector3d(0.0, 0.0, 0.0), Vector3d(0.0, 0.0, 1.0))
xbc0 = ReflectingBoundary3d(xPlane0)
ybc0 = ReflectingBoundary3d(yPlane0)
zbc0 = ReflectingBoundary3d(zPlane0)
hydro.appendBoundary(xbc0)
hydro.appendBoundary(ybc0)
hydro.appendBoundary(zbc0)
output("hydro.haveBoundary(xbc0)")
output("hydro.haveBoundary(ybc0)")
output("hydro.haveBoundary(zbc0)")

#-------------------------------------------------------------------------------
# Construct a predictor corrector integrator.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
output("integrator")
integrator.appendPhysicsPackage(hydro)
output("integrator.havePhysicsPackage(hydro)")
output("integrator.valid()")
integrator.lastDt = dt
output("integrator.lastDt")
if dtMin:
    integrator.dtMin = dtMin
    output("integrator.dtMin")
if dtMax:
    integrator.dtMax = dtMax
    output("integrator.dtMax")
integrator.dtGrowth = dtGrowth
output("integrator.dtGrowth")

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName)
output("control")

#-------------------------------------------------------------------------------
# Load a restart file if requested.
#-------------------------------------------------------------------------------
if restoreCycle is not None:
    control.loadRestartFile(restoreCycle)
else:
    control.iterateIdealH()
    control.smoothState(smoothIters)

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
control.step(maxSteps)
