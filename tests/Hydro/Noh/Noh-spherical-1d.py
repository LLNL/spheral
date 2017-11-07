from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from SpheralVisitDump import dumpPhysicsState

title("3-D spherical hydro test -- Spherical Noh problem")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
NodeListConstructor = SphNodeList3d

nx1 = 100
rho1 = 1.0
eps1 = 0.0
x0, x1 = 0.0, 1.0
nPerh = 2.01

vr0, vrSlope = -1.0, 0.0

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

goalTime = 0.6
dt = 0.0001
dtMin, dtMax = 1.0e-5, None
dtGrowth = 2.0
maxSteps = None
statsStep = 1
smoothIters = 0
HEvolution = Hydro3d.HEvolutionType.IdealH
sumForMassDensity = Hydro3d.MassDensityType.IntegrateDensity # RigorousSumDensity

restoreCycle = None
restartStep = 10000
restartBaseName = "Noh-spherical-1d-%i" % nx1

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
neighborTimer = SpheralTimer("Neighbor initialization.")
neighborTimer.start()
neighbor1 = NestedGridNeighbor3d(nodes1,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes1.registerNeighbor(neighbor1)
neighborTimer.stop()
neighborTimer.printStatus()

#-------------------------------------------------------------------------------
# Set node positions for this domain
#-------------------------------------------------------------------------------
from DistributeNodes import distributeNodesInSphericalRange3d
distributeNodesInSphericalRange3d([(nodes1, nx1, rho1, (x0, x1))])
output("nodes1.numNodes")

#-------------------------------------------------------------------------------
# Set node specific thermal energies
#-------------------------------------------------------------------------------
nodes1.setSpecificThermalEnergy(ScalarField3d(nodes1.specificThermalEnergy().name(), nodes1, eps1))

#-------------------------------------------------------------------------------
# Set node velocities
#-------------------------------------------------------------------------------
for ix in xrange(nodes1.numNodes):
    nodes1.velocity()[ix] = Vector3d(vr0 + vrSlope*nodes1.positions()[ix].x, 0.0, 0.0)

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
bc0 = SphericalBoundary(db)
hydro.appendBoundary(bc0)
output("hydro.haveBoundary(bc0)")

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

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
control.step(5)
control.advance(goalTime, maxSteps)
control.dropRestartFile()

#-------------------------------------------------------------------------------
# Plot the final state.
#-------------------------------------------------------------------------------
rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db, plotStyle="lines")

#-------------------------------------------------------------------------------
# Overplot the analytic solution.
#-------------------------------------------------------------------------------
import mpi
import NohAnalyticSolution
rlocal = [pos.x for pos in nodes1.positions().internalValues()]
r = mpi.reduce(rlocal, mpi.SUM)
h1 = 1.0/(nPerh*(x1 - x0)/nx1)
answer = NohAnalyticSolution.NohSolution(3,
                                         r = r,
                                         v0 = -1.0,
                                         h0 = 1.0/h1)
plotAnswer(answer, control.time(), rhoPlot, velPlot, epsPlot, PPlot, HPlot)

#-------------------------------------------------------------------------------
# Compute the error.
#-------------------------------------------------------------------------------
def filterOnRadius(data, r, rmin, rmax):
    assert len(data) == len(r)
    return [data[i] for i in xrange(len(data)) if (r[i].magnitude() >= rmin and r[i].magnitude() <= rmax)]

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
