from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from SpheralVisitDump import dumpPhysicsState
from GenerateNodeDistribution2d import *
from findLastRestart import *
from DistributeNodes import *
import mpi

title("3-D hydro test -- Spherical Noh problem in RZ coordinates.")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(NodeListConstructor = SphNodeListRZ,

            seed = "lattice", # "constantDTheta",
            theta = 0.5*pi,
            nRadial = 20,
            nTheta = 20,
            rmin = 0.0,
            rmax = 1.0,
            rho1 = 1.0,
            eps1 = 0.0,
            vr1 = -1.0,
            nPerh = 2.01,

            gamma = 5.0/3.0,
            mu = 1.0,
            Qconstructor = MonaghanGingoldViscosityRZ,
            Cl = 1.0,
            Cq = 0.75,
            Qlimiter = False,
            epsilon2 = 1e-2,
            negligibleSoundSpeed = 1e-5,
            csMultiplier = 1e-4,
            energyMultiplier = 0.1,
            hmin = 0.0001,
            hmax = 0.1,
            HsmoothFraction = 0.0,
            cfl = 0.5,
            XSPH = True,

            neighborSearchType = Neighbor2d.NeighborSearchType.GatherScatter,
            numGridLevels = 20,
            topGridCellSize = 2.0,
            origin = Vector2d(0.0, 0.0),

            goalTime = 0.6,
            dtSample = 0.1,
            dt = 0.0001,
            dtMin = 1.0e-5,
            dtMax = None,
            dtGrowth = 2.0,
            maxSteps = None,
            statsStep = 1,
            smoothIters = 0,
            HEvolution = Hydro2d.HEvolutionType.IdealH,
            sumForMassDensity = Hydro2d.MassDensityType.RigorousSumDensity, # IntegrateDensity, # 
            compatibleEnergy = False,

            restartStep = 10000,
            )

dataDir = "Noh-spherical-rz-%i" % nRadial
restartDir = dataDir + "/restarts"
visitDir = dataDir + "/visit"
restartBaseName = restartDir + "/Noh-spherical-rz"

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
eos = GammaLawGasMKS2d(gamma, mu)

#-------------------------------------------------------------------------------
# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
#-------------------------------------------------------------------------------
WT = TableKernel2d(BSplineKernel2d(), 1000)
WTPi = TableKernel2d(BSplineKernel2d(), 1000)
kernelExtent = WT.kernelExtent()
output("WT")
output("WTPi")

#-------------------------------------------------------------------------------
# Create an empty NodeList
#-------------------------------------------------------------------------------
nodes1 = NodeListConstructor("nodes1", eos, WT, WTPi)
nodes1.hmin = hmin
nodes1.hmax = hmax
nodes1.HsmoothFraction = HsmoothFraction
nodes1.nodesPerSmoothingScale = nPerh
nodes1.XSPH = XSPH
output("nodes1")
output("nodes1.HsmoothFraction")
output("nodes1.nodesPerSmoothingScale")
output("nodes1.XSPH")

#-------------------------------------------------------------------------------
# Construct the neighbor object and associate it with the node list.
#-------------------------------------------------------------------------------
neighbor1 = NestedGridNeighbor2d(nodes1,
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
    generator = GenerateNodeDistributionRZ(nRadial, nTheta, rho1, seed,
                                           xmin = (0.0, 0.0),
                                           xmax = (1.0, 1.0),
                                           theta = theta,
                                           nNodePerh = nPerh,
                                           SPH = NodeListConstructor is SphNodeListRZ)
##     generator = GenerateNodeDistributionRZ(nRadial, nTheta, rho1, seed,
##                                            rmin = rmin,
##                                            rmax = rmax,
##                                            theta = theta,
##                                            nNodePerh = nPerh,
##                                            SPH = NodeListConstructor is SphNodeListRZ)
    n1 = generator.globalNumNodes()

    print "Distributing nodes amongst processors."
    distributeNodes2d((nodes1, generator))
    output("mpi.reduce(nodes1.numInternalNodes, mpi.MIN)")
    output("mpi.reduce(nodes1.numInternalNodes, mpi.MAX)")
    output("mpi.reduce(nodes1.numInternalNodes, mpi.SUM)")

    # Specific thermal energies
    nodes1.specificThermalEnergy(ScalarField2d("tmp", nodes1, eps1))

    # Node velocities
    nodes1.velocity(VectorField2d("tmp", nodes1, Vector2d(-1.0, 0.0)))
##     for i in xrange(nodes1.numInternalNodes):
##         runit = nodes1.positions()[i].unitVector()
##         nodes1.velocity()[i] = runit*vr1

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase2d()
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
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
hydro = Hydro2d(WT, WTPi, q, compatibleEnergy)
hydro.cfl = cfl
hydro.HEvolution = HEvolution
hydro.sumForMassDensity = sumForMassDensity
hydro.HsmoothMin = hmin
hydro.HsmoothMax = hmax
output("hydro")
output("hydro.compatibleEnergyEvolution")
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
xPlane0 = Plane2d(Vector2d(0.0, 0.0), Vector2d(1.0, 0.0))
yPlane0 = Plane2d(Vector2d(0.0, 0.0), Vector2d(0.0, 1.0))
yPlane1 = Plane2d(Vector2d(0.0, 1.0), Vector2d(0.0, -1.0))
xbc = ReflectingBoundary2d(xPlane0)
ybc = ReflectingBoundary2d(yPlane0)
ybc1 = ReflectingBoundary2d(yPlane1)
hydro.appendBoundary(xbc)
hydro.appendBoundary(ybc)
hydro.appendBoundary(ybc1)
output("hydro.haveBoundary(xbc)")
output("hydro.haveBoundary(ybc)")

#-------------------------------------------------------------------------------
# Construct a predictor corrector integrator, and add the one physics package.
#-------------------------------------------------------------------------------
integrator = SynchronousRK2Integrator2d(db)
integrator.appendPhysicsPackage(hydro)
integrator.lastDt = dt
integrator.dtGrowth = dtGrowth
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.valid()")
output("integrator.lastDt")
output("integrator.dtGrowth")
if dtMin:
    integrator.dtMin = dtMin
    output("integrator.dtMin")
if dtMax:
    integrator.dtMax = dtMax
    output("integrator.dtMax")

control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            initializeMassDensity = (hydro.sumForMassDensity ==
                                                     Hydro2d.MassDensityType.RigorousSumDensity))
output("control")

#-------------------------------------------------------------------------------
# Load a restart file if requested.
#-------------------------------------------------------------------------------
if restoreCycle:
    control.loadRestartFile(restoreCycle)
else:
    control.iterateIdealH()
    control.smoothState(smoothIters)
    dumpPhysicsState(integrator,
                     "Noh-spherical-rz",
                     visitDir)

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
while control.time() < goalTime:
    nextGoalTime = min(int((control.time() + 1.001*dtSample)/dtSample)*dtSample,
                       goalTime)
    control.advance(nextGoalTime, maxSteps)
    control.dropRestartFile()
    dumpPhysicsState(integrator,
                     "Noh-spherical-rz",
                     visitDir)

#-------------------------------------------------------------------------------
# Plot the final state.
#-------------------------------------------------------------------------------
rplot = plotNodePositions2d(db)
rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db,
                                                    xFunction = "%s.x",
                                                    vecyFunction = "%s.x")
#plotRadialState(db)

for p in [rhoPlot, velPlot, epsPlot, PPlot, HPlot]:
    p("set xrange [0:]")
    p.refresh()

#-------------------------------------------------------------------------------
# Overplot the analytic solution.
#-------------------------------------------------------------------------------
import mpi
import NohAnalyticSolution
rlocal = [pos.x for pos in nodes1.positions().internalValues()]
r = mpi.reduce(rlocal, mpi.SUM)
h1 = 1.0/(nPerh*(rmax - rmin)/nRadial)
answer = NohAnalyticSolution.NohSolution(2,
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
