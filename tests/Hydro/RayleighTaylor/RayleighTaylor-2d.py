#-------------------------------------------------------------------------------
# RayleighTaylor-2d.py
#
# Performs a 2-D simulation of the the growth a a Rayleigh-Taylor instability
# in an initially hydrostatically balanced atmosphere in a gravitational field.
# An initial spatial perturbation is applied at the interface by displacing the
# nodes according to:
#
#  dy = a0 exp(-((y - y0)/b0)^2) cos(2 pi x lambda0/xbox + pi)
#     a0: max amplitude of perturbation
#     b0: scale for depth of pertubation away from interface
#     y0: y position of interface
#     lambda0: number of perturbations in box
#     xbox: x width of the box
#
# Kull, H.J., 1991, Physics Reports (Review section of Physics Letters),
#   206, No. 5, 197-325.
#-------------------------------------------------------------------------------
from math import *
from Spheral import *
from SpheralTestUtilities import *
from findLastRestart import *
from SpheralVisitDump import dumpPhysicsState

# Load the mpi module.
import loadmpi
mpi, procID, numProcs = loadmpi.loadmpi()

from GenerateNodeDistribution2d import *

#-------------------------------------------------------------------------------
# Identify ourselves.
#-------------------------------------------------------------------------------
title("2-D Rayleigh-Taylor gravitational instability.")

#-------------------------------------------------------------------------------
# Functor to provide the exponential initial density field.
#-------------------------------------------------------------------------------
class ExponentialDensity:
    def __init__(self,
                 y0,
                 rho0,
                 alpha):
        self.y0 = y0
        self.rho0 = rho0
        self.alpha = alpha
        return
    def __call__(self, r):
        return self.rho0*exp(self.alpha*(r.y - self.y0))

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
seed = "lattice"
NodeListConstructor = SphNodeList2d

nx1, ny1 = 50, 50
nx2, ny2 = 50, 50
nPerh = 2.01
numybound = 2.5*nPerh

y0 = 0.5            # Position of the interface
xmin1 = (0.0, 0.0)
xmax1 = (0.5, y0)
xmin2 = (0.0, y0)
xmax2 = (0.5, 1.0)

gamma = 5.0/3.0
mu = 1.0
S = 3.0             # density jump at interface
A = (S - 1)/(S + 1) # Atwood number
eps0 = 1.0          # initial specific thermal energy (hot/low density side)
rho0 = 1.0          # initial mass density (cold/high density side)
g0 = -2.0           # constant gravitational acceleration

lambda0 = 1.0       # number of wavelengths in box
a0 = 0.001          # initial amplitude of perturbation
b0 = 0.1            # exponential scale for decay of perturbation away from
                    # interface

#Qconstructor = MonaghanGingoldViscosity2d
Qconstructor = TensorMonaghanGingoldViscosity2d
Cl, Cq = 0.5, 1.0
Qlimiter = True
balsaraCorrection = False
epsilon2 = 1e-4
negligibleSoundSpeed = 1e-5
csMultiplier = 1e-4
HsmoothMin, HsmoothMax, HratioMin = 0.000001, 0.5, 0.1
cfl = 0.5
XSPH = True
epsilonTensile = 0.0
nTensile = 4

neighborSearchType = Neighbor2d.NeighborSearchType.GatherScatter
numGridLevels = 20
topGridCellSize = 1.0
origin = Vector2d(0.0, 0.0)

goalTime = 3.7
dtSample = 0.05
dt = 0.001
dtMin, dtMax = 1.0e-5, None
dtGrowth = 2.0
maxSteps = None
statsStep = 10
smoothIters = 0
HEvolution = Hydro2d.HEvolutionType.IdealH
sumForMassDensity = Hydro2d.MassDensityType.VolumeScaledDensity # IntegrateDensity # 

restartStep = 1000
dataDir = "RayleighTaylor-%3.1fwave-2d/%ix%i-XSPH=%i" % (lambda0, nx1, ny1 + ny2, XSPH)
restartDir = dataDir + "/restarts"
visitDir = dataDir + "/visit"
restartBaseName = restartDir + "/RayleighTaylor-2d-%ix%i" % (nx1, ny1 + ny2)

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
# Material parameters.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS2d(gamma, mu)

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
# Create the NodeLists.
#-------------------------------------------------------------------------------
nodes1 = NodeListConstructor("Low density nodes", eos, WT, WTPi)
nodes2 = NodeListConstructor("High density nodes", eos, WT, WTPi)
nodeSet = [nodes1, nodes2]
for nodes in nodeSet:
    nodes.XSPH = XSPH
    nodes.nodesPerSmoothingScale = nPerh
    nodes.epsilonTensile = epsilonTensile
    nodes.nTensile = nTensile
    output("nodes.name()")
    output("  nodes.nodesPerSmoothingScale")
    output("  nodes.epsilonTensile")
    output("  nodes.nTensile")
    output("  nodes.XSPH")

#-------------------------------------------------------------------------------
# Construct the neighbor object and associate it with the node list.
#-------------------------------------------------------------------------------
_cache = []
for nodes in nodeSet:
    neighbor = NestedGridNeighbor2d(nodes,
                                    neighborSearchType,
                                    numGridLevels,
                                    topGridCellSize,
                                    origin,
                                    kernelExtent)
    nodes.registerNeighbor(neighbor)
    _cache.append(neighbor)

#-------------------------------------------------------------------------------
# Set node positions, masses, and H's for this domain.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    from DistributeNodes import distributeNodes2d
    print "Generating node distribution."
    generator1 = GenerateNodeDistribution2d(nx1, ny1,
                                            ExponentialDensity(y0,
                                                               rho0/S,
                                                               g0/((gamma - 1.0)*eps0)),
                                            seed,
                                            xmin = xmin1,
                                            xmax = xmax1,
                                            nNodePerh = nPerh,
                                            SPH = NodeListConstructor is SphNodeList2d)
    n1 = generator1.globalNumNodes()
    generator2 = GenerateNodeDistribution2d(nx2, ny2,
                                            ExponentialDensity(y0,
                                                               rho0,
                                                               g0*S/((gamma - 1.0)*eps0)),
                                            seed,
                                            xmin = xmin2,
                                            xmax = xmax2,
                                            nNodePerh = nPerh,
                                            SPH = NodeListConstructor is SphNodeList2d)
    n2 = generator2.globalNumNodes()

    print "Distributing nodes amongst processors."
    distributeNodes2d([(nodes1, generator1),
                       (nodes2, generator2)])
    for nodes in nodeSet:
        output("mpi.reduce(nodes.numInternalNodes, mpi.MIN)")
        output("mpi.reduce(nodes.numInternalNodes, mpi.MAX)")
        output("mpi.reduce(nodes.numInternalNodes, mpi.SUM)")

    # Set node specific thermal energies
    nodes1.specificThermalEnergy(ScalarField2d("tmp", nodes1, eps0))
    nodes2.specificThermalEnergy(ScalarField2d("tmp", nodes2, eps0/S))

    # Perturb the initial node positions.
    for nodes in nodeSet:
        for i in xrange(nodes.numInternalNodes):
            ri = nodes.positions()[i]
            dy = a0*exp(-((ri.y - y0)/b0)**2)*cos(2.0*pi*lambda0*ri.x/xmax1[0] + pi)
            nodes.positions()[i].y += dy

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase2d()
output("db")
output("db.appendNodeList(nodes1)")
output("db.appendNodeList(nodes2)")
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
# Construct the gravitational acceleration object.
#-------------------------------------------------------------------------------
nodeIndicies1 = vector_of_int()
nodeIndicies2 = vector_of_int()
nodeIndicies1.extend(range(nodes1.numInternalNodes))
nodeIndicies2.extend(range(nodes2.numInternalNodes))
gravity1 = ConstantAcceleration2d(Vector2d(0.0, g0),
                                  nodes1,
                                  nodeIndicies1)
gravity2 = ConstantAcceleration2d(Vector2d(0.0, g0),
                                  nodes2,
                                  nodeIndicies2)

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
hydro = Hydro2d(WT, WTPi, q)
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
# Create boundary conditions and huck them onto the hydro object.
#-------------------------------------------------------------------------------
xPlane0 = Plane2d(Vector2d(*xmin1), Vector2d(1.0, 0.0))
xPlane1 = Plane2d(Vector2d(*xmax2), Vector2d(-1.0, 0.0))
xbc0 = ReflectingBoundary2d(xPlane0)
xbc1 = ReflectingBoundary2d(xPlane1)

dy1 = (xmax1[1] - xmin1[1])/ny1
dy2 = (xmax2[1] - xmin2[1])/ny2
ynodes1 = vector_of_int()
ynodes2 = vector_of_int()
for i in xrange(nodes1.numInternalNodes):
    if nodes1.positions()[i].y < xmin1[1] + numybound*dy1:
        ynodes1.append(i)
for i in xrange(nodes2.numInternalNodes):
    if nodes2.positions()[i].y > xmax2[1] - numybound*dy2:
        ynodes2.append(i)
ybc1 = ConstantVelocityBoundary2d(nodes1, ynodes1)
ybc2 = ConstantVelocityBoundary2d(nodes2, ynodes2)
print ("Selected (%i, %i) constant velocity nodes." %
       (mpi.allreduce(len(ynodes1), mpi.SUM),
        mpi.allreduce(len(ynodes2), mpi.SUM)))

for bc in [xbc0, xbc1, ybc1, ybc2]:
    hydro.appendBoundary(bc)

#-------------------------------------------------------------------------------
# Construct a predictor corrector integrator, and add the physics packages.
#-------------------------------------------------------------------------------
integrator = PredictorCorrectorIntegrator2d(db)
integrator.appendPhysicsPackage(gravity1)
integrator.appendPhysicsPackage(gravity2)
integrator.appendPhysicsPackage(hydro)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
output("integrator")
output("integrator.havePhysicsPackage(gravity1)")
output("integrator.havePhysicsPackage(gravity2)")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.lastDt")
output("integrator.dtGrowth")
output("integrator.valid()")

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            initializeMassDensity = False)
output("control")

#-------------------------------------------------------------------------------
# Smooth the initial conditions.
#-------------------------------------------------------------------------------
if restoreCycle is not None:
    control.loadRestartFile(restoreCycle)
else:
    control.iterateIdealH(hydro,
                          idealHTolerance = 2.0e-10)
    control.smoothState(smoothIters)
    dumpPhysicsState(integrator,
                     "RayleighTaylor-2d",
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
                     "RayleighTaylor-2d",
                     visitDir)

#-------------------------------------------------------------------------------
# Plot the final positions.
#-------------------------------------------------------------------------------
##try:
##    import pylab
##    from SpheralMatplotlibUtilities import *
##    rplot = plotNodePositions2d(db)
##    rhoplot = plotFieldList(db.fluidMassDensity,
##                            xFunction = "%s.y",
##                            plotStyle = "o",
##                            winTitle = "Mass density")
##    Pplot = plotFieldList(db.fluidPressure,
##                          xFunction = "%s.y",
##                          plotStyle = "o",
##                          winTitle = "Pressure")
##    pylab.show()
##except:
##    print "Unable to display inline graphics."
##    pass
