from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from findLastRestart import *
from SpheralVisitDump import dumpPhysicsState
from math import *
import mpi

from GenerateNodeDistribution3d import *

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("3-D integrated hydro test -- Shearing Planar Noh problem")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
seed = "lattice"
NodeListConstructor = AsphNodeList3d

nx, ny, nz = 20, 100, 20
nPerh = 1.51

xmin = (0.0, 0.0, 0.0)
xmax = (0.2, 1.0, 0.2)

rho1 = 1.0
eps1 = 0.0
vshear = 1.0
vy = -1.0

gamma = 5.0/3.0
mu = 1.0
Qconstructor = MonaghanGingoldViscosity3d
#Qconstructor = TensorMonaghanGingoldViscosity3d
Cl, Cq = 0.5, 1.0
Qlimiter = False
balsaraCorrection = True
epsilon2 = 1e-4
negligibleSoundSpeed = 1e-5
csMultiplier = 1e-4
HsmoothMin, HsmoothMax, HratioMin = 0.0004, 0.5, 0.1
HsmoothFraction = 0.0
cfl = 0.5
XSPH = True
epsilonTensile = 0.0
nTensile = 8

neighborSearchType = Neighbor3d.NeighborSearchType.GatherScatter
numGridLevels = 20
topGridCellSize = 0.5
origin = Vector3d(0.0, 0.0, 0.0)

goalTime = 0.6
dt = 0.0001
dtMin, dtMax = 1.0e-5, None
dtGrowth = 2.0
maxSteps = None
statsStep = 10
smoothIters = 0
HEvolution = Hydro3d.HEvolutionType.IdealH
sumForMassDensity = Hydro3d.MassDensityType.RigorousSumDensity

restartStep = 1000
dataDir = "shear-%ix%ix%i" % (nx, ny, nz)
restartDir = dataDir + "/restarts"
visitDir = dataDir + "/visit"
restartBaseName = restartDir + "/Noh-shear-3d-%ix%ix%i" % (nx, ny, nz)

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
# Gamma law gas EOS.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS3d(gamma, mu)

#-------------------------------------------------------------------------------
# Create an empty NodeList
#-------------------------------------------------------------------------------
nodes1 = NodeListConstructor("Regular Nodes", eos)
nodes1.HsmoothFraction = HsmoothFraction
nodes1.XSPH = XSPH
nodes1.nodesPerSmoothingScale = nPerh
nodes1.epsilonTensile = epsilonTensile
nodes1.nTensile = nTensile
output("nodes1")
output("nodes1.HsmoothFraction")
output("nodes1.nodesPerSmoothingScale")
output("nodes1.epsilonTensile")
output("nodes1.nTensile")
output("nodes1.XSPH")

#-------------------------------------------------------------------------------
# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
#-------------------------------------------------------------------------------
WT = TableKernel3d(BSplineKernel3d(), 1000)
WTPi = TableKernel3d(BSplineKernel3d(), 1000)
output("WT")
output("WTPi")
kernelExtent = WT.kernelExtent()

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
# Set node positions, masses, and H's for this domain.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    from DistributeNodes import distributeNodes3d
    print "Generating node distribution."
    generator1 = GenerateNodeDistribution3d(nx, ny, nz, rho1, seed,
                                            xmin = xmin,
                                            xmax = xmax,
                                            nNodePerh = nPerh,
                                            SPH = (NodeListConstructor is SphNodeList3d))
    n1 = generator1.globalNumNodes()

    print "Distributing nodes amongst processors."
    nodeInfo = distributeNodes3d([(nodes1, n1, generator1)])
    if mpi:
        output("mpi.reduce(nodes1.numInternalNodes, mpi.MIN)")
        output("mpi.reduce(nodes1.numInternalNodes, mpi.MAX)")
        output("mpi.reduce(nodes1.numInternalNodes, mpi.SUM)")
    else:
        output("nodes1.numInternalNodes")
    assert len(nodeInfo[nodes1]["globalNodeListID"]) == nodes1.numInternalNodes

    # Set initial node mass densities.
    nodes1.setMassDensity(ScalarField3d("tmp", nodes1, rho1))
    nodes1.updateWeight()

    # Set node specific thermal energies
    nodes1.setSpecificThermalEnergy(ScalarField3d("tmp", nodes1, eps1))

    # Set node velocities
    for i in xrange(nodes1.numInternalNodes):
        y = nodes1.positions()[i].y
        nodes1.velocity()[i] = Vector3d(vshear*cos(2.0*pi*y), vy, 0.0)

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
# Create boundary conditions.
xPlane0 = Plane3d(Vector3d(*xmin), Vector3d(1.0, 0.0, 0.0))
xPlane1 = Plane3d(Vector3d(*xmax), Vector3d(-1.0, 0.0, 0.0))
yPlane0 = Plane3d(Vector3d(*xmin), Vector3d(0.0, 1.0, 0.0))
zPlane0 = Plane3d(Vector3d(*xmin), Vector3d(0.0, 0.0, 1.0))
zPlane1 = Plane3d(Vector3d(*xmax), Vector3d(0.0, 0.0, -1.0))
xbc0 = PeriodicBoundary3d(xPlane0, xPlane1)
ybc0 = ReflectingBoundary3d(yPlane0)
zbc0 = PeriodicBoundary3d(zPlane0, zPlane1)
hydro.appendBoundary(xbc0)
hydro.appendBoundary(ybc0)
hydro.appendBoundary(zbc0)
output("hydro.haveBoundary(xbc0)")
output("hydro.haveBoundary(ybc0)")
output("hydro.haveBoundary(zbc0)")

#-------------------------------------------------------------------------------
# Construct a predictor corrector integrator, and add the physics packages.
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

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName)
output("control")

#-------------------------------------------------------------------------------
# Smooth the initial conditions/restore state.
#-------------------------------------------------------------------------------
if restoreCycle is not None:
    control.loadRestartFile(restoreCycle)
else:
    control.smoothState(smoothIters)

    # Viz the initial conditions.
    dumpPhysicsState(integrator,
                     "Noh-shear-2d-%ix%i-visit" % (nx, ny),
                     visitDir)

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
dtSample = 0.1
while control.time() < goalTime:
    nextGoalTime = min(int((control.time() + 1.0001*dtSample)/dtSample)*dtSample,
                       goalTime)
    control.advance(nextGoalTime, maxSteps)
    control.dropRestartFile()
    dumpPhysicsState(integrator,
                     "Noh-shear-3d-%ix%ix%i-visit" % (nx, ny, nz),
                     visitDir)

# Plot the elongation (h1/h2) for the H tensors.
import Gnuplot
hratio0 = [H.eigenValues().minElement()/H.eigenValues().maxElement() for H in nodes1.Hfield().internalValues()]
y0 = [ri.y for ri in nodes1.positions().internalValues()]
hratio = mpi.allreduce(hratio0, mpi.SUM)
y = mpi.allreduce(y0, mpi.SUM)
hratioPlot = Gnuplot.Gnuplot()
cache = []
if mpi.rank == 0:
    hratioData = Gnuplot.Data(y, hratio,
                              title = "hmin/hmax ratio",
                              inline = True)
    hratioPlot.plot(hratioData)
    cache.extend([hratioData, hratioPlot])

rPlot = plotNodePositions3d(db, colorNodeLists=0, colorDomains=1)

# Plot the final state.
vxPlot = plotFieldList(db.fluidVelocity, xFunction="%s.y", yFunction="%s.x", plotStyle="points", winTitle="x velocity")
vyPlot = plotFieldList(db.fluidVelocity, xFunction="%s.y", yFunction="%s.y", plotStyle="points", winTitle="y velocity")
vzPlot = plotFieldList(db.fluidVelocity, xFunction="%s.y", yFunction="%s.z", plotStyle="points", winTitle="z velocity")
rhoPlot = plotFieldList(db.fluidMassDensity, xFunction="%s.y", plotStyle="points", winTitle="Mass Density")
P = db.fluidPressure
PPlot = plotFieldList(P, xFunction="%s.y", plotStyle="points", winTitle="Pressure")
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
hxPlot = plotFieldList(hx, xFunction="%s.y", plotStyle="points", winTitle="h_x")
hyPlot = plotFieldList(hy, xFunction="%s.y", plotStyle="points", winTitle="h_y")
hzPlot = plotFieldList(hz, xFunction="%s.y", plotStyle="points", winTitle="h_z")

# Overplot the analytic solution.
from NohAnalyticSolution import *
answer = NohSolution(1,
                     h0 = nPerh*1.0/ny)
plotAnswer(answer, control.time(),
           rhoPlot = rhoPlot,
           velPlot = vyPlot,
           PPlot = PPlot,
           HPlot = hyPlot)
