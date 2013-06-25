#-------------------------------------------------------------------------------
# A goofy hacked version of the 1-D planar Sod problem set to run with 3-D
# Spheral++ objects, using our hacked support for running a single line of
# 3-D nodes in the x direction.
#-------------------------------------------------------------------------------
from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from SodAnalyticSolution import *

################################################################################
# Generic problem parameters
nx1, nx2 = 400, 100
rho1, rho2 = 1.0, 0.25
P1, P2 = 1.0, 0.1795
x0, x1, x2 = -0.5, 0.0, 0.5

nPerh1, nPerh2 = 2.01, 2.01

m1 = (x1 - x0)*rho1/nx1
m2 = (x2 - x1)*rho2/nx2

gamma = 5.0/3.0
mu = 1.0
Qconstructor = MonaghanGingoldViscosity3d
#Qconstructor = TensorMonaghanGingoldViscosity3d
Cl, Cq = 1.0, 0.75
Qlimiter = 0
epsilon2 = 1e-4
negligibleSoundSpeed = 1e-5
csMultiplier = 1e-4
energyMultiplier = 1.0
HsmoothMin, HsmoothMax = 0.0001, 1.0
cfl = 0.5

neighborSearchType = Neighbor3d.NeighborSearchType.GatherScatter
numGridLevels = 15
topGridCellSize = 2.0
origin = Vector3d(-1.0, 0.0, 0.0)

goalTime = 0.15
dt = 1e-4
dtMin, dtMax = 1.0e-5, 0.1
dtGrowth = 2.0
maxSteps = None
statsStep = 10
smoothIters = 0
HEvolution = Integrator3d.HEvolutionType.IdealH
sumForMassDensity = Integrator3d.MassDensityType.RigorousSum

restoreCycle = None
restartStep = 200
restartBaseName = "Sod-planar-3d-restart"

################################################################################
title('1-D integrated hydro test -- planar Sod problem')

eos = GammaLawGasMKS3d(gamma, mu)

# Create two empty NodeLists
nodes1 = AsphNodeListThreeOneDimension3d(eos)
nodes1.nodesPerSmoothingScale = nPerh1
output('nodes1')
output('nodes1.nodesPerSmoothingScale')

nodes2 = AsphNodeListThreeOneDimension3d(eos)
nodes2.nodesPerSmoothingScale = nPerh2
output('nodes2')
output('nodes2.nodesPerSmoothingScale')

# Set node positions
from DistributeNodes import distributeNodesInRange1d
distributeNodesInRange1d([(nodes1, nx1, (x0, x1)),
                          (nodes2, nx2, (x1, x2))])
nNodesThisDomain1 = nodes1.numInternalNodes
nNodesThisDomain2 = nodes2.numInternalNodes
output('nodes1.numNodes')
output('nodes2.numNodes')

# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
# Note that we use the hacked TableKernel constructor here to indicate we're
# going to be running a 1-D problem with 3-D objects.
WT = TableKernel3d(BSplineKernel3d(), BSplineKernel1d(), 100)
WTPi = WT
output('WT')
output('WTPi')
kernelExtent = WT.kernelExtent()

# Set node masses
nodes1.setMass(ScalarField3d(nodes1, m1))
nodes2.setMass(ScalarField3d(nodes2, m2))

# Set node specific thermal energies
eps1 = P1/((gamma - 1.0)*rho1)
eps2 = P2/((gamma - 1.0)*rho2)
nodes1.setSpecificThermalEnergy(ScalarField3d(nodes1, eps1))
nodes2.setSpecificThermalEnergy(ScalarField3d(nodes2, eps2))

# Set the smoothing scales.
dx1 = (x1 - x0)/nx1
dx2 = (x2 - x1)/nx2
h1 = 1.0/(nPerh1*dx1)
h2 = 1.0/(nPerh2*dx2)
H1 = SymTensor3d(h1, 0.0, 0.0,
                 0.0, 1.0, 0.0,
                 0.0, 0.0, 1.0)
H2 = SymTensor3d(h2, 0.0, 0.0,
                 0.0, 1.0, 0.0,
                 0.0, 0.0, 1.0)
nodes1.setHfield(SymTensorField3d(nodes1, H1))
nodes2.setHfield(SymTensorField3d(nodes2, H2))

# Set the mass densities if required.
nodes1.setMassDensity(ScalarField3d(nodes1, rho1))
nodes2.setMassDensity(ScalarField3d(nodes2, rho2))

# Construct the neighbor object and associate it with the node list.
neighborTimer = SpheralTimer('Neighbor initialization.')
neighborTimer.start()
neighbor1 = NestedGridNeighbor3d(nodes1,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
neighbor2 = NestedGridNeighbor3d(nodes2,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes1.registerNeighbor(neighbor1)
nodes2.registerNeighbor(neighbor2)
neighborTimer.stop()
neighborTimer.printStatus()

# Create boundary conditions.  We need at least this much to create the initial
# mass density field.
xPlane0 = Plane3d(Vector3d(x0, 0.0, 0.0), Vector3d( 1.0, 0.0, 0.0))
xPlane1 = Plane3d(Vector3d(x2, 0.0, 0.0), Vector3d(-1.0, 0.0, 0.0))
xbc0 = ReflectingBoundary3d(xPlane0)
xbc1 = ReflectingBoundary3d(xPlane1)
xbc2 = ThreeOneDimensionBoundary3d()

# Construct a DataBase to hold our node list
db = DataBase3d()
output('db')
output('db.appendNodeList(nodes1)')
output('db.appendNodeList(nodes2)')
output('db.numNodeLists')
output('db.numFluidNodeLists')

# Construct the artificial viscosities for the problem.
q = Qconstructor(Cl, Cq)
q.limiter = Qlimiter
q.epsilon2 = epsilon2
q.negligibleSoundSpeed = negligibleSoundSpeed
q.csMultiplier = csMultiplier
q.energyMultiplier = energyMultiplier
output('q')
output('q.Cl')
output('q.Cq')
output('q.limiter')
output('q.epsilon2')
output('q.negligibleSoundSpeed')
output('q.csMultiplier')
output('q.energyMultiplier')

# Construct the hydro physics object.
hydro = Hydro3d(WT, WTPi, q)
hydro.cfl = cfl
output('hydro')
output('hydro.kernel')
output('hydro.PiKernel')
output('hydro.cfl')
output('hydro.valid()')

# Construct a predictor corrector integrator, and add the one physics package.
integrator = PredictorCorrectorIntegrator3d(db)
output('integrator')
integrator.appendPhysicsPackage(hydro)
output('integrator.havePhysicsPackage(hydro)')
output('integrator.valid()')
integrator.HsmoothMin = HsmoothMin
integrator.HsmoothMax = HsmoothMax
output('integrator.HsmoothMin')
output('integrator.HsmoothMax')
integrator.lastDt = dt
output('integrator.lastDt')
integrator.dtGrowth = dtGrowth
if dtMin:
    integrator.dtMin = dtMin
    output('integrator.dtMin')
if dtMax:
    integrator.dtMax = dtMax
    output('integrator.dtMax')
integrator.HEvolution = HEvolution
integrator.sumForMassDensity = sumForMassDensity
if (sumForMassDensity == Integrator3d.MassDensityType.Sum or
    sumForMassDensity == Integrator3d.MassDensityType.RigorousSum):
    integrator.setKernel(WT)
output('integrator.HEvolution')
output('integrator.sumForMassDensity')

control = SpheralController(integrator, WT,
                            boundaryConditions = [xbc0, xbc1, xbc2],
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName)
output('control')

# Smooth the initial conditions.
if restoreCycle:
    control.loadRestartFile(restoreCycle)
else:
    control.smoothState(smoothIters)

##################################################################################
# Advance to the end time.
control.advance(goalTime, maxSteps)
control.dropRestartFile()

# Plot the final state.
rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db)

# Now overplot the analytic solution.
answer = SodSolution(nPoints=nx1 + nx2,
                     gamma = gamma,
                     rho1 = rho1,
                     P1 = P1,
                     rho2 = rho2,
                     P2 = P2,
                     x0 = x0,
                     x1 = x1,
                     x2 = x2,
                     h1 = 1.0/h1,
                     h2 = 1.0/h2)

plotAnswer(answer, control.time(),
           rhoPlot, velPlot, epsPlot, PPlot, HPlot)
