from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
NodeListConstructor = SphNodeList1d

nxtot = 100
numNodeLists = 4
nx1 = nxtot/numNodeLists
assert nx1*numNodeLists == nxtot

rho1 = 1.0
eps1 = 0.0
x0, x1 = 0.0, 1.0
dxnodes = (x1 - x0)/numNodeLists
nPerh = 2.01

vx1 = -1.0

gamma = 5.0/3.0
mu = 1.0
Qconstructor = MonaghanGingoldViscosity1d
#Qconstructor = TensorMonaghanGingoldViscosity1d
Cl, Cq = 1.0, 0.75
Qlimiter = 0
epsilon2 = 1e-4
negligibleSoundSpeed = 1e-5
csMultiplier = 1e-4
energyMultiplier = 0.1
HsmoothMin, HsmoothMax = 0.0001, 0.1
HsmoothFraction = 0.0
cfl = 0.5
XSPH = True
epsilonTensile = 0.0
nTensile = 8
epsilonTensileGradient = 0.01

neighborSearchType = Neighbor1d.NeighborSearchType.GatherScatter
numGridLevels = 20
topGridCellSize = 2.0
origin = Vector1d(0.0)

goalTime = 0.6
dt = 0.0001
dtMin, dtMax = 1.0e-5, None #0.09
dtGrowth = 2.0
maxSteps = None
statsStep = 1
smoothIters = 0
HEvolution = Hydro1d.HEvolutionType.IdealH
sumForMassDensity = Hydro1d.MassDensityType.RigorousSumDensity

restoreCycle = None
restartStep = 10000
restartBaseName = "Noh-planar-1d-%i" % nx1

#-------------------------------------------------------------------------------
title("1-D integrated hydro test -- planar Noh problem broken up into multiple NodeLists")

eos = GammaLawGasMKS1d(gamma, mu)

# Create the NodeLists.
nodeLists = []
nodeInfo = []
for i in xrange(numNodeLists):
    nodes = NodeListConstructor("Nodes %i" % i, eos)
    nodes.HsmoothFraction = HsmoothFraction
    nodes.nodesPerSmoothingScale = nPerh
    nodes.XSPH = XSPH
    nodes.epsilonTensile = epsilonTensile
    nodes.nTensile = nTensile
    nodes.epsilonTensileGradient = epsilonTensileGradient
    nodeLists.append(nodes)
    nodeInfo.append((nodes, nx1, rho1, (x0 + i*dxnodes, x0 + (i + 1)*dxnodes)))
    output("nodes.name()")
    output("nodes.HsmoothFraction")
    output("nodes.nodesPerSmoothingScale")
    output("nodes.XSPH")
    output("nodes.epsilonTensile")
    output("nodes.nTensile")
    output("nodes.epsilonTensileGradient")

# Set node positions for this domain
from DistributeNodes import distributeNodesInRange1d
distributeNodesInRange1d(nodeInfo)

# Set the node properties.
for nodes in nodeLists:

    # Set node specific thermal energies
    nodes.setSpecificThermalEnergy(ScalarField1d(nodes.specificThermalEnergy().name(), nodes, eps1))

    # Set node velocities
    nodes.setVelocity(VectorField1d(nodes.velocity().name(), nodes, Vector1d(vx1)))

# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
WT = TableKernel1d(BSplineKernel1d(), 100)
WTPi = TableKernel1d(BSplineKernel1d(), 100)
output("WT")
output("WTPi")
kernelExtent = WT.kernelExtent()

# Construct the neighbor objects.
neighborTimer = SpheralTimer("Neighbor initialization.")
neighborTimer.start()
_cache = []
for nodes in nodeLists:
    neighbor = NestedGridNeighbor1d(nodes,
                                    neighborSearchType,
                                    numGridLevels,
                                    topGridCellSize,
                                    origin,
                                    kernelExtent)
    nodes.registerNeighbor(neighbor)
    _cache.append(neighbor)
neighborTimer.stop()
neighborTimer.printStatus()

# Construct a DataBase to hold our node list
db = DataBase1d()
output("db")
for nodes in nodeLists:
    output("db.appendNodeList(nodes)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

# Construct the artificial viscosities for the problem.
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

# Construct the hydro physics object.
hydro = Hydro1d(WT, WTPi, q)
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

# Create boundary conditions.
xPlane0 = Plane1d(Vector1d(x0), Vector1d(1.0))
xbc0 = ReflectingBoundary1d(xPlane0)
hydro.appendBoundary(xbc0)
output("hydro.haveBoundary(xbc0)")

# Construct a predictor corrector integrator, and add the one physics package.
integrator = PredictorCorrectorIntegrator1d(db)
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
                            restartBaseName = restartBaseName)
output("control")

# Smooth the initial conditions.
if restoreCycle is not None:
    control.loadRestartFile(restoreCycle)
else:
    control.smoothState(smoothIters)

##################################################################################
# Advance to the end time.
control.step(5)
control.advance(goalTime, maxSteps)

# Plot the final state.
rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db, plotStyle="lines")

