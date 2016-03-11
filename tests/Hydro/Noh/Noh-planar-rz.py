from Spheral import *
from SpheralTestUtilities import *
from SpheralMatplotlibUtilities import *
from SpheralVisitDump import dumpPhysicsState
from GenerateNodeDistribution2d import *
from findLastRestart import *
from DistributeNodes import *
import mpi

title("3-D hydro test -- Planar Noh problem in RZ coordinates.")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
NodeListConstructor = SphNodeList3d

seed = "lattice"
nz, nr = 100, 20
xmin = (0.0, 0.0)
xmax = (1.0, 0.2)
rho1 = 1.0
eps1 = 0.0
v1 = -1.0
nPerh = 2.01

gamma = 5.0/3.0
mu = 1.0
#Qconstructor = MonaghanGingoldViscosity3d
Qconstructor = TensorMonaghanGingoldViscosity3d
Cl, Cq = 1.0, 0.75
Qlimiter = True
epsilon2 = 1e-4
negligibleSoundSpeed = 1e-5
csMultiplier = 1e-4
energyMultiplier = 0.1
HsmoothMin, HsmoothMax = 0.0001, 0.1
HsmoothFraction = 0.0
cfl = 0.5
XSPH = False

neighborSearchType = Neighbor3d.NeighborSearchType.GatherScatter
numGridLevels = 20
topGridCellSize = 2.0
origin = Vector3d(0.0, 0.0, 0.0)

goalTime = 0.02
dt = 0.0001
dtMin, dtMax = 1.0e-5, None
dtGrowth = 2.0
maxSteps = None
statsStep = 1
smoothIters = 0
HEvolution = Hydro3d.HEvolutionType.IdealH
sumForMassDensity = Hydro3d.MassDensityType.IntegrateDensity # RigorousSumDensity

restartStep = 10000
dataDir = "Noh-planar-rz-%ix%i" % (nz, nr)
restartDir = dataDir + "/restarts"
visitDir = dataDir + "/visit"
restartBaseName = restartDir + "/Noh-planar-rz"

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
restoreCycle = findLastRestart(restartBaseName)

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS3d(gamma, mu)

#-------------------------------------------------------------------------------
# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
#-------------------------------------------------------------------------------
WT = TableKernel3d(BSplineKernel3d(), 1000)
WTPi = TableKernel3d(BSplineKernel3d(), 1000)
kernelExtent = WT.kernelExtent()
output("WT")
output("WTPi")

#-------------------------------------------------------------------------------
# Create an empty NodeList
#-------------------------------------------------------------------------------
nodes1 = NodeListConstructor("nodes1", eos)
nodes1.HsmoothFraction = HsmoothFraction
nodes1.nodesPerSmoothingScale = nPerh
output("nodes1")
output("nodes1.HsmoothFraction")
output("nodes1.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Construct the neighbor object and associate it with the node list.
#-------------------------------------------------------------------------------
neighbor1 = NestedGridNeighbor3d(nodes1,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes1.registerNeighbor(neighbor1)

#-------------------------------------------------------------------------------
# Set node properties.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    generator = GenerateNodeDistributionRZ(nz, nr, rho1, seed,
                                           xmin = xmin,
                                           xmax = xmax,
                                           nNodePerh = nPerh,
                                           SPH = NodeListConstructor is SphNodeList3d)
    n1 = generator.globalNumNodes()

    print "Distributing nodes amongst processors."
    nodeInfo = distributeNodes3d([(nodes1, n1, generator)])
    output("mpi.reduce(nodes1.numInternalNodes, mpi.MIN)")
    output("mpi.reduce(nodes1.numInternalNodes, mpi.MAX)")
    output("mpi.reduce(nodes1.numInternalNodes, mpi.SUM)")
    assert len(nodeInfo[nodes1]["globalNodeListID"]) == nodes1.numInternalNodes

    # Specific thermal energies
    nodes1.setSpecificThermalEnergy(ScalarField3d("tmp", nodes1, eps1))

    # Node velocities
    nodes1.setVelocity(VectorField3d("tmp", nodes1, Vector3d(v1, 0.0, 0.0)))

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
q.epsilon2 = epsilon2
q.negligibleSoundSpeed = negligibleSoundSpeed
q.csMultiplier = csMultiplier
q.energyMultiplier = energyMultiplier
output("q")
output("q.Cl")
output("q.Cq")
output("q.limiter")
output("q.epsilon2")
output("q.negligibleSoundSpeed")
output("q.csMultiplier")
output("q.energyMultiplier")

#-------------------------------------------------------------------------------
# Set the XSPH and tensile corrections for the NodeList
#-------------------------------------------------------------------------------
nodes1.XSPH = XSPH
output("nodes1.XSPH")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
hydro = Hydro3d(WT, WTPi, q)
hydro.cfl = cfl
hydro.HEvolution = HEvolution
hydro.sumForMassDensity = sumForMassDensity
hydro.HsmoothMin = HsmoothMin
hydro.HsmoothMax = HsmoothMax
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
yPlane0 = Plane3d(Vector3d(0.0, xmax[1], 0.0), Vector3d(0.0, -1.0, 0.0))
xrefbc = ReflectingBoundary3d(xPlane0)
yrefbc = ReflectingBoundary3d(yPlane0)
cylbc = CylindricalBoundary(db)
hydro.appendBoundary(cylbc)
hydro.appendBoundary(xrefbc)
hydro.appendBoundary(yrefbc)
output("hydro.haveBoundary(cylbc)")
output("hydro.haveBoundary(xrefbc)")
output("hydro.haveBoundary(yrefbc)")

#-------------------------------------------------------------------------------
# Construct a predictor corrector integrator, and add the one physics package.
#-------------------------------------------------------------------------------
integrator = PredictorCorrectorIntegrator3d(db)
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

control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            initializeMassDensity = False)
output("control")

#-------------------------------------------------------------------------------
# Load a restart file if requested.
#-------------------------------------------------------------------------------
if restoreCycle:
    control.loadRestartFile(restoreCycle)
else:
    control.smoothState(smoothIters)
    dumpPhysicsState(integrator,
                     "Noh-planar-rz",
                     visitDir)
    dumpPhysicsState(integrator,
                     "Noh-planar-rz-ghost",
                     visitDir,
                     dumpGhosts = True)

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if control.time() < goalTime:
    control.advance(goalTime, maxSteps)
    control.dropRestartFile()
    dumpPhysicsState(integrator,
                     "Noh-planar-rz",
                     visitDir)
    dumpPhysicsState(integrator,
                     "Noh-planar-rz-ghost",
                     visitDir,
                     dumpGhosts = True)

#-------------------------------------------------------------------------------
# Plot the final state.
#-------------------------------------------------------------------------------
rplot = plotNodePositions2d(db)
rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db)

#-------------------------------------------------------------------------------
# Overplot the analytic solution.
#-------------------------------------------------------------------------------
import mpi
import NohAnalyticSolution
rlocal = [pos.x for pos in nodes1.positions().internalValues()]
r = mpi.reduce(rlocal, mpi.SUM)
h1 = 1.0/(nPerh*(xmax[0] - xmin[0])/nz)
answer = NohAnalyticSolution.NohSolution(1,
                                         r = r,
                                         v0 = -1.0,
                                         h0 = 1.0/h1)
plotAnswer(answer, control.time(), rhoPlot, velPlot, epsPlot, PPlot, HPlot)

#-------------------------------------------------------------------------------
# Compute the error.
#-------------------------------------------------------------------------------
def filterOnRadius(data, r, rmin, rmax):
    assert len(data) == len(r)
    return [data[i] for i in xrange(len(data)) if (r[i].x >= rmin and r[i].x <= rmax)]

rmin, rmax = 0.05, 0.35   # Throw away anything with r < rwall to avoid wall heating.
rhoprof = mpi.reduce(filterOnRadius(nodes1.massDensity().internalValues(), nodes1.positions().internalValues(), rmin, rmax), mpi.SUM)
Pprof = mpi.reduce(filterOnRadius(nodes1.pressure().internalValues(), nodes1.positions().internalValues(), rmin, rmax), mpi.SUM)
vprof = mpi.reduce(filterOnRadius([v.x for v in nodes1.velocity().internalValues()], nodes1.positions().internalValues(), rmin, rmax), mpi.SUM)
epsprof = mpi.reduce(filterOnRadius(nodes1.specificThermalEnergy().internalValues(), nodes1.positions().internalValues(), rmin, rmax), mpi.SUM)
hprof = mpi.reduce(filterOnRadius([1.0/H.xx for H in nodes1.Hfield().internalValues()], nodes1.positions().internalValues(), rmin, rmax), mpi.SUM)
if mpi.rank == 0:
    answer.r = [r for r in answer.r if r >= rmin and r <= rmax]
    rans, vans, epsans, rhoans, Pans, hans = answer.solution(control.time())
    import Pnorm
    print "\tQuantity \t\tL1 \t\t\tL2 \t\t\tLinf"
    for (name, data, ans) in [("Mass Density", rhoprof, rhoans),
                              ("Pressure", Pprof, Pans),
                              ("Velocity", vprof, vans),
                              ("Thermal E", epsprof, epsans),
                              ("h       ", hprof, hans)]:
        assert len(data) == len(ans)
        error = [data[i] - ans[i] for i in xrange(len(data))]
        Pn = Pnorm.Pnorm(error, hprof)
        print "\t%s \t\t%g \t\t%g \t\t%g" % (name,
                                             Pn.gridpnorm(1),
                                             Pn.gridpnorm(2),
                                             Pn.gridpnorm("inf"))
