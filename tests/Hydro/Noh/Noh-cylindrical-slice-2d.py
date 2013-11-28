from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from findLastRestart import *
from SpheralVisitDump import dumpPhysicsState
from math import *

# Load the mpi module if we're parallel.
try:
    import mpi
except:
    mpi = None

from GenerateNodeDistribution2d import *

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
title("2-D integrated hydro test -- cylindrical Noh problem")

seed = "constantDTheta" # "optimal" # "constantNTheta"
NodeListConstructor = AsphNodeList2d

theta = 0.25*pi
nRadial, nTheta = 40, 40
nPerh = 2.01

rmin, rmax = 0.0, 1.0
xmin, xmax = (0.0, 0.0), (rmax, rmax)

rho1 = 1.0
eps1 = 0.0
vr1 = -1.0
x0, y0 = 0.0, 0.0
x1, y1 = 1.0, 1.0

gamma = 5.0/3.0
mu = 1.0
#Qconstructor = MonaghanGingoldViscosity2d
Qconstructor = TensorMonaghanGingoldViscosity2d
Cl, Cq = 0.5, 1.0
Qlimiter = True
balsaraCorrection = 0
epsilon2 = 1e-4
negligibleSoundSpeed = 1e-5
csMultiplier = 1e-4
HsmoothMin, HsmoothMax, HratioMin = 0.0004, 0.5, 0.1
HsmoothFraction = 0.0
cfl = 0.5
XSPH = True
epsilonTensile = 0.0
nTensile = 4

neighborSearchType = Neighbor2d.NeighborSearchType.GatherScatter
numGridLevels = 16
topGridCellSize = 0.5
origin = Vector2d(0.0, 0.0)

goalTime = 0.6
dt = 0.0001
dtMin, dtMax = 1.0e-5, 0.1
dtGrowth = 2.0
useVelocityForDt = False
maxSteps = None
statsStep = 10
smoothIters = 0
HEvolution = Hydro2d.HEvolutionType.IdealH  #IntegrateH
sumForMassDensity = Hydro2d.MassDensityType.RigorousSumDensity

restartStep = 1000
dataDir = "cylindrical-slice-%ix%i" % (nRadial, nTheta)
restartDir = dataDir + "/restarts"
visitDir = dataDir + "/visit"
restartBaseName = restartDir + "Noh-cylindrical-2d-%ix%i" % (nRadial, nTheta)

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
# Material properties
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS2d(gamma, mu)

#-------------------------------------------------------------------------------
# Create an empty NodeList
#-------------------------------------------------------------------------------
nodes1 = NodeListConstructor("Regular Nodes", eos)
output("nodes1")
nodes1.HsmoothFraction = HsmoothFraction
nodes1.XSPH = XSPH
nodes1.nodesPerSmoothingScale = nPerh
nodes1.epsilonTensile = epsilonTensile
nodes1.nTensile = nTensile
output("nodes1.HsmoothFraction")
output("nodes1.nodesPerSmoothingScale")
output("nodes1.epsilonTensile")
output("nodes1.nTensile")
output("nodes1.XSPH")

#-------------------------------------------------------------------------------
# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
#-------------------------------------------------------------------------------
WT = TableKernel2d(BSplineKernel2d(), 1000)
WTPi = TableKernel2d(BSplineKernel2d(), 1000)
output("WT")
output("WTPi")
kernelExtent = WT.kernelExtent()

#-------------------------------------------------------------------------------
# Construct the neighbor object and associate it with the node list.
#-------------------------------------------------------------------------------
neighborTimer = SpheralTimer("Neighbor initialization.")
neighborTimer.start()
neighbor1 = NestedGridNeighbor2d(nodes1,
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
    import DistributeNodes
    print "Generating node distribution."
    generator = GenerateNodeDistribution2d(nRadial, nTheta, rho1, seed,
                                           rmin = rmin,
                                           rmax = rmax,
                                           theta = theta,
                                           nNodePerh = nPerh)
    n1 = generator.globalNumNodes()

    print "Distributing nodes amongst processors."
    nodeInfo = DistributeNodes.distributeNodes2d([(nodes1, n1, generator)])
    if mpi:
        output("mpi.reduce(nodes1.numInternalNodes, mpi.MIN)")
        output("mpi.reduce(nodes1.numInternalNodes, mpi.MAX)")
        output("mpi.reduce(nodes1.numInternalNodes, mpi.SUM)")
    else:
        output("nodes1.numInternalNodes")
    assert len(nodeInfo[nodes1]["globalNodeListID"]) == nodes1.numInternalNodes

    # Set node specific thermal energies
    nodes1.setSpecificThermalEnergy(ScalarField2d("tmp", nodes1, eps1))

    # Set node velocities
    for nodeID in xrange(nodes1.numNodes):
        nodes1.velocity()[nodeID] = nodes1.positions()[nodeID].unitVector()*vr1

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase2d()
output("db")
output("db.appendNodeList(nodes1)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct a standard Monaghan-Gingold artificial viscosity.
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
hydro = Hydro2d(WT, WTPi, q)
hydro.cfl = cfl
hydro.HEvolution = HEvolution
hydro.sumForMassDensity = sumForMassDensity
hydro.useVelocityMagnitudeForDt = useVelocityForDt
hydro.HsmoothMin = HsmoothMin
hydro.HsmoothMax = HsmoothMax
hydro.HratioMin = HratioMin
output("hydro")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.cfl")
output("hydro.HEvolution")
output("hydro.sumForMassDensity")
output("hydro.useVelocityMagnitudeForDt")
output("hydro.HsmoothMin")
output("hydro.HsmoothMax")
output("hydro.HratioMin")
output("hydro.valid()")

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
reflectPlane1 = Plane2d(Vector2d(0.0, 0.0),
                        Vector2d(0.0, 1.0))
reflectPlane2 = Plane2d(Vector2d(0.0, 0.0),
                        Vector2d(sin(theta), -cos(theta)))
reflectPlane3 = Plane2d(Vector2d(0.0, 0.0),
                        Vector2d(sin(pi - theta), -cos(pi - theta)))
reflectbc1 = ReflectingBoundary2d(reflectPlane1)
reflectbc2 = ReflectingBoundary2d(reflectPlane2)
reflectbc3 = ReflectingBoundary2d(reflectPlane3)
hydro.appendBoundary(reflectbc1)
hydro.appendBoundary(reflectbc2)
hydro.appendBoundary(reflectbc3)
output("hydro.haveBoundary(reflectbc1)")
output("hydro.haveBoundary(reflectbc2)")
output("hydro.haveBoundary(reflectbc3)")

#-------------------------------------------------------------------------------
# Construct a predictor corrector integrator, and add the physics packages.
#-------------------------------------------------------------------------------
integrator = PredictorCorrectorIntegrator2d(db)
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

# Smooth the initial conditions.
if restoreCycle is not None:
    control.loadRestartFile(restoreCycle)
else:
    control.smoothState(smoothIters)
    dumpPhysicsState(integrator,
                     "Noh-cylindrical-2d-%ix%i" % (nRadial, nTheta),
                     visitDir)

################################################################################
# Advance to the end time.
dtdump = 0.6
nextGoalTime = min(goalTime, control.time() + dtdump)
while control.time() < goalTime:
    while control.time() < nextGoalTime:
        control.advance(nextGoalTime, maxSteps)
        dumpPhysicsState(integrator,
                         "Noh-cylindrical-2d-%ix%i" % (nRadial, nTheta),
                         visitDir)
    nextGoalTime = min(goalTime, control.time() + dtdump)
    if restoreCycle < control.totalSteps:
        control.dropRestartFile()

# Plot the elongation (h1/h2) for the H tensors.
import Gnuplot
##hratio0 = [H.eigenValues().minElement()/H.eigenValues().maxElement() for H in nodes1.Hfield().internalValues()]
##r0 = [ri.magnitude() for ri in nodes1.positions().internalValues()]
##hratio = mpi.allreduce(hratio0, mpi.SUM)
##r = mpi.allreduce(r0, mpi.SUM)
##hratioPlot = Gnuplot.Gnuplot()
##cache = []
##if mpi.rank == 0:
##    hratioData = Gnuplot.Data(r, hratio,
##                              title = "hmin/hmax ratio",
##                              inline = True)
##    hratioPlot.plot(hratioData)
##    cache.extend([hratioData, hratioPlot])

rPlot = plotNodePositions2d(db, colorNodeLists=0, colorDomains=1)

# Plot the final state.
rhoPlot, vrPlot, epsPlot, PPlot, HPlot = plotRadialState(db)
del HPlot
Hinverse = db.fluidHinverse
hr = db.newFluidScalarFieldList()
ht = db.newFluidScalarFieldList()
for Hfield, hrfield, htfield in zip(Hinverse.fields(),
                                    hr.fields(),
                                    ht.fields()):
    n = Hfield.numElements()
    assert hrfield.numElements() == n
    assert htfield.numElements() == n
    positions = Hfield.nodeList().positions()
    for i in xrange(n):
        runit = positions[i].unitVector()
        tunit = Vector2d(-(positions[i].y), positions[i].x).unitVector()
        hrfield[i] = (Hfield[i]*runit).magnitude()
        htfield[i] = (Hfield[i]*tunit).magnitude()
hrPlot = plotFieldList(hr, xFunction="%s.magnitude()", plotStyle="points", winTitle="h_r")
htPlot = plotFieldList(ht, xFunction="%s.magnitude()", plotStyle="points", winTitle="h_t")

# Overplot the analytic solution.
import NohAnalyticSolution
answer = NohAnalyticSolution.NohSolution(2,
                                         h0 = nPerh*rmax/nRadial)
plotAnswer(answer, control.time(),
           rhoPlot = rhoPlot,
           velPlot = vrPlot,
           epsPlot = epsPlot,
           PPlot = PPlot,
           HPlot = hrPlot)

if mpi.rank == 0:
    r, hrans, htans = answer.hrtsolution(control.time())
    htData = Gnuplot.Data(r, htans,
                          title = "Solution",
                          with = "lines",
                          inline = "true")
    htPlot.replot(htData)
