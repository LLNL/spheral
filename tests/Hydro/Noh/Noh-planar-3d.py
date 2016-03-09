5#-------------------------------------------------------------------------------
# The Planar Noh test case run in 3-D.
#
# W.F. Noh 1987, JCP, 72, 78-120.
#-------------------------------------------------------------------------------
from math import *
from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from findLastRestart import *
from SpheralVisitDump import dumpPhysicsState
import mpi

from GenerateNodeDistribution3d import *

title("2-D integrated hydro test -- Planar Noh problem")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(NodeListConstructor = AsphNodeList3d,

            seed = "lattice",

            nx = 100,
            ny = 20,
            nz = 20,

            rho1 = 1.0,
            eps1 = 0.0,
            vx = -1.0,
            vshear = 0.0,
            nPerh = 2.01,

            gamma = 5.0/3.0,
            mu = 1.0,
            Qconstructor = MonaghanGingoldViscosity3d,
            #Qconstructor = TensorMonaghanGingoldViscosity3d,
            Cl = 1.0,
            Cq = 1.5,
            Qlimiter = False,
            balsaraCorrection = False,
            epsilon2 = 1e-2,
            negligibleSoundSpeed = 1e-5,
            csMultiplier = 1e-4,
            energyMultiplier = 0.1,
            hmin = 1e-5,
            hmax = 1.0,
            hminratio = 0.05,
            HsmoothFraction = 0.0,
            cfl = 0.5,
            XSPH = True,

            HEvolution = Hydro3d.HEvolutionType.IdealH,
            limitIdealH = False,

            neighborSearchType = Neighbor3d.NeighborSearchType.GatherScatter,
            numGridLevels = 20,
            topGridCellSize = 2.0,
            origin = Vector3d(0.0, 0.0),

            goalTime = 0.6,
            dt = 0.0001,
            dtMin = 1.0e-5, 
            dtMax = None,
            dtGrowth = 2.0,
            dtSample = 0.1,
            rigorousBoundaries = False,
            maxSteps = None,
            statsStep = 1,
            smoothIters = 0,
            sumForMassDensity = Hydro3d.MassDensityType.RigorousSumDensity, # VolumeScaledDensity,
            compatibleEnergy = True,

            restoreCycle = None,
            restartStep = 1000,

            # Parameters for the test acceptance.,
            L1rho = 0.702502,
            L2rho = 1.01547,
            Linfrho = 3.07796,
            
            L1P = 0.194579,
            L2P = 0.341756,
            LinfP = 1.3412,
            
            L1v = 0.339706,
            L2v = 0.478579,
            Linfv = 1.10007,
            
            L1eps = 0.111857,
            L2eps = 0.170687,
            Linfeps = 0.493335,

            L1h = 0.00698501,
            L2h = 0.00828045,
            Linfh = 0.0149375,

            tol = 1.0e-5,

            graphics = "gnu",
            )

xmin = (0.0, 0.0, 0.0)
xmax = (1.0, 0.2, 0.2)

dataDir = "dumps-planar-%ix%ix%i" % (nx, ny, nz)
restartDir = dataDir + "/restarts"
visitDir = dataDir + "/visit"
restartBaseName = restartDir + "/Noh-planar-3d-%ix%i" % (nx, ny)

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
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel3d(BSplineKernel3d(), 100)
WTPi = TableKernel3d(BSplineKernel3d(), 100)
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
nodes1.hmin = hmin
nodes1.hmax = hmax
nodes1.hminratio = hminratio
output("nodes1.HsmoothFraction")
output("nodes1.nodesPerSmoothingScale")
output("nodes1.XSPH")
output("nodes1.hmin")
output("nodes1.hmax")
output("nodes1.hminratio")

#-------------------------------------------------------------------------------
# Construct the neighbor object.
#-------------------------------------------------------------------------------
neighbor1 = NestedGridNeighbor3d(nodes1,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes1.registerNeighbor(neighbor1)

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    from ParMETISDistributeNodes import distributeNodes3d
    generator1 = GenerateNodeDistribution3d(nx, ny, nz, rho1, seed,
                                            xmin = xmin,
                                            xmax = xmax,
                                            nNodePerh = nPerh)
    distributeNodes3d((nodes1, generator1))
    output("mpi.reduce(nodes1.numInternalNodes, mpi.MIN)")
    output("mpi.reduce(nodes1.numInternalNodes, mpi.MAX)")
    output("mpi.reduce(nodes1.numInternalNodes, mpi.SUM)")

    # Set node specific thermal energies
    nodes1.specificThermalEnergy(ScalarField3d("tmp", nodes1, eps1))

    # Set node velocities
    for i in xrange(nodes1.numInternalNodes):
        x = nodes1.positions()[i].x
        nodes1.velocity()[i] = Vector3d(vx,
                                        vshear*cos(2.0*pi*x),
                                        vshear*sin(2.0*pi*x))

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
hydro = Hydro3d(WT, WTPi, q, compatibleEnergy)
hydro.cfl = cfl
hydro.HEvolution = HEvolution
hydro.limitIdealH = limitIdealH
hydro.sumForMassDensity = sumForMassDensity
hydro.HsmoothMin = hmin
hydro.HsmoothMax = hmax
output("hydro")
output("hydro.cfl")
output("hydro.HEvolution")
output("hydro.sumForMassDensity")
output("hydro.HsmoothMin")
output("hydro.HsmoothMax")
output("hydro.compatibleEnergyEvolution")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.valid()")

packages = [hydro]

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane0 = Plane3d(Vector3d(*xmin), Vector3d(1.0, 0.0, 0.0))

yPlane0 = Plane3d(Vector3d(*xmin), Vector3d(0.0, 1.0, 0.0))
yPlane1 = Plane3d(Vector3d(*xmax), Vector3d(0.0, -1.0, 0.0))

zPlane0 = Plane3d(Vector3d(*xmin), Vector3d(0.0, 0.0, 1.0))
zPlane1 = Plane3d(Vector3d(*xmax), Vector3d(0.0, 0.0, -1.0))

xbc = ReflectingBoundary3d(xPlane0)
ybc = PeriodicBoundary3d(yPlane0, yPlane1)
zbc = PeriodicBoundary3d(zPlane0, zPlane1)

for p in packages:
    p.appendBoundary(xbc)
    p.appendBoundary(ybc)
    p.appendBoundary(zbc)

#-------------------------------------------------------------------------------
# Construct a time integrator.
#-------------------------------------------------------------------------------
integrator = SynchronousRK2Integrator3d(db)
for p in packages:
    integrator.appendPhysicsPackage(p)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.rigorousBoundaries = rigorousBoundaries
output("integrator")
output("integrator.valid()")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.rigorousBoundaries")

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName)
output("control")

# Smooth the initial conditions.
if restoreCycle is not None:
    control.loadRestartFile(restoreCycle)
else:
    control.iterateIdealH()
    control.smoothState(smoothIters)
    control.dropRestartFile()
    dumpPhysicsState(integrator,
                     "Noh-planar-3d-%ix%i" % (nx, ny),
                     visitDir)

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
hstats([nodes1])
while control.time() < goalTime:
    nextGoalTime = min(control.time() + dtSample, goalTime)
    control.advance(nextGoalTime, maxSteps)
    control.dropRestartFile()
    dumpPhysicsState(integrator,
                     "Noh-planar-3d-%ix%ix%i" % (nx, ny, nz),
                     visitDir)
hstats([nodes1])

#-------------------------------------------------------------------------------
# Plot the elongation (h1/h2) for the H tensors.
#-------------------------------------------------------------------------------
import Gnuplot
hratio0 = [H.eigenValues().minElement()/H.eigenValues().maxElement() for H in nodes1.Hfield().internalValues()]
x0 = [ri.x for ri in nodes1.positions().internalValues()]
hratio = mpi.allreduce(hratio0, mpi.SUM)
x = mpi.allreduce(x0, mpi.SUM)
hratioPlot = Gnuplot.Gnuplot()
cache = []
if mpi.rank == 0:
    hratioData = Gnuplot.Data(x, hratio,
                              title = "hmin/hmax ratio",
                              inline = True)
    hratioPlot.plot(hratioData)
    cache.extend([hratioData, hratioPlot])

# Plot the final state.
vxPlot = plotFieldList(db.fluidVelocity, xFunction="%s.x", yFunction="%s.x", plotStyle="points", winTitle="x velocity")
vyPlot = plotFieldList(db.fluidVelocity, xFunction="%s.x", yFunction="%s.y", plotStyle="points", winTitle="y velocity")
vzPlot = plotFieldList(db.fluidVelocity, xFunction="%s.x", yFunction="%s.z", plotStyle="points", winTitle="z velocity")
rhoPlot = plotFieldList(db.fluidMassDensity, xFunction="%s.y", plotStyle="points", winTitle="Mass Density")
P = db.fluidPressure
PPlot = plotFieldList(P, xFunction="%s.x", plotStyle="points", winTitle="Pressure")
Hinverse = db.fluidHinverse
hx = db.newFluidScalarFieldList()
hy = db.newFluidScalarFieldList()
hz = db.newFluidScalarFieldList()
xunit = Vector3d(1, 0, 0)
yunit = Vector3d(0, 1, 0)
zunit = Vector3d(0, 0, 1)
for Hfield, hxfield, hyfield, hzfield in zip(Hinverse.fields(),
                                             hx.fields(),
                                             hy.fields(),
                                             hz.fields()):
    n = Hfield.numElements()
    assert hxfield.numElements() == n and hyfield.numElements() == n and hzfield.numElements() == n
    for i in xrange(n):
        hxfield[i] = (Hfield[i]*xunit).magnitude()
        hyfield[i] = (Hfield[i]*yunit).magnitude()
        hzfield[i] = (Hfield[i]*zunit).magnitude()
hxPlot = plotFieldList(hx, xFunction="%s.x", plotStyle="points", winTitle="h_x")
hyPlot = plotFieldList(hy, xFunction="%s.x", plotStyle="points", winTitle="h_y")
hzPlot = plotFieldList(hy, xFunction="%s.z", plotStyle="points", winTitle="h_z")

# Overplot the analytic solution.
from NohAnalyticSolution import *
answer = NohSolution(1,
                     h0 = nPerh*1.0/nx)
plotAnswer(answer, control.time(),
           rhoPlot = rhoPlot,
           velPlot = vyPlot,
           PPlot = PPlot,
           HPlot = hxPlot)
