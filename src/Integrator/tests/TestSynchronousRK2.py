from Spheral import *
from SpheralTestUtilities import *

################################################################################
def plotState(nodes, color='black', plotGhosts=0):
    from SpheralGistUtilities import *
    from gist import *

    if plotGhosts:
        nx = nodes.numNodes
    else:
        nx = nodes.numInternalNodes

    window(0)
    xNodes = nodePositions1d(nodes)[:nx]
    rhoNodes = array(nodes.massDensity[:nx])
    plg(rhoNodes, xNodes, color=color)
    pltitle('Mass density')

    window(1)
    vNodes = array([0.0]*nx)
    for i in xrange(nx):
        vNodes[i] = nodes.velocity[i].x
    plg(vNodes, xNodes, color=color)
    pltitle('Velocity')

    window(2)
    pressure = nodes.pressure
    PNodes = array(pressure[:nx])
    plg(PNodes, xNodes, color=color)
    pltitle('Pressure')

    window(3)
    HNodes = array([0.0]*nx)
    for i in xrange(nx):
        HNodes[i] = 1.0/nodes.Hfield[i].xx
    plg(HNodes, xNodes, color=color)
    pltitle('Smoothing scale')

################################################################################
# Generic problem parameters
nx1, nx2 = 50, 50
rho1, rho2 = 1.0, 1.0
m1, m2 = 0.5*rho1/nx1, 0.5*rho2/nx2
P1, P2 = 1.0, 1.0
x0, x1, x2 = -0.5, 0.0, 0.5

gamma = 1.4
mu = 1.0
Cl = 0.75
Cq = 1.5
epsilon2 = 1e-2
HsmoothMin, HsmoothMax = 0.0001, 0.1
cfl = 0.1

neighborSearchType = 3 # GatherScatter
numGridLevels = 10
topGridCellSize = 0.25
origin = Vector1d(0.0)

goalTime = 0.5
dtMin, dtMax = 1e-5, 0.1
dt = 0.0001
maxSteps = 500
smoothIters = 0
sumForMassDensity = 0

################################################################################
title('1-D integrated hydro test -- planar Sod problem')

nx = nx1 + nx2

eos = GammaLawGasMKS1d(gamma, mu)

nodes1 = SphNodeList1d(nx1 + nx2, eos)
output('nodes1.numNodes')

W = BSplineKernel1d()
#W = W4SplineKernel1d()
#W = GaussianKernel1d()
#W = SuperGaussianKernel1d()
#W = PiGaussianKernel1d(1.0)
output('W')
kernelExtent = W.kernelExtent

# Set node positions
dx1 = (x1 - x0)/nx1
for ix in xrange(nx1 + nx2):
    nodeID = ix
    nodes1.positions[nodeID] = x0 + (ix + 0.5)*dx1

dx2 = (x2 - x1)/nx2
for ix in xrange(nx2):
    nodeID = ix + nx1
    nodes1.positions[nodeID] = x1 + (ix + 0.5)*dx2

# Set node masses
nodes1.mass[:nx1] = [m1]*nx1
nodes1.mass[nx1:] = [m2]*nx2

# Set node specific thermal energies
eps1 = P1/((gamma - 1.0)*rho1)
for nodeID in xrange(nx1):
    nodes1.specificThermalEnergy[nodeID] = eps1

eps2 = P2/((gamma - 1.0)*rho2)
for nodeID in xrange(nx1, nx1 + nx2):
    nodes1.specificThermalEnergy[nodeID] = eps2

# Set node velocities
for nodeID in xrange(nodes1.numNodes):
    nodes1.velocity[nodeID] = (0.0)

# Set the smoothing scales.
h1 = 1.0/(2.01*dx1)
h2 = 1.0/(2.01*dx2)
for i in xrange(nx1):
    nodes1.Hfield[i].xx = h1
for i in xrange(nx2, nx1 + nx2):
    nodes1.Hfield[i].xx = h2

# Set the mass densities if required.
nodes1.massDensity[:nx1] = [rho1]*nx1
nodes1.massDensity[nx1:] = [rho2]*nx2

# Construct the neighbor object and associate it with the node list.
neighborTimer = SpheralTimer('Neighbor initialization.')
neighborTimer.start()
neighbor1 = NestedGridNeighbor1d(nodes1,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes1.neighbor = neighbor1
neighborTimer.stop()
neighborTimer.printStatus()

# Create boundary conditions.  We need at least this much to create the initial
# mass density field.
xPlane0 = Plane1d((x0), ( 1.0))
xPlane1 = Plane1d((x2), (-1.0))
xbc0 = ReflectingBoundary1d(xPlane0)
xbc1 = ReflectingBoundary1d(xPlane1)

# Construct a DataBase to hold our node list
db = DataBase1d()
output('db')
output('db.appendNodeList(nodes1)')
output('db.numNodeLists')
output('db.numFluidNodeLists')

# Construct a standard Monaghan-Gingold artificial viscosity.
q = MonaghanGingoldViscosity1d(Cl, Cq)
q.epsilon2 = epsilon2
output('q')
output('q.epsilon2')

# Construct the hydro physics object.
hydro = FakeHydro1d(q)
hydro.cfl = cfl
output('hydro')
output('hydro.valid')
output('hydro.cfl')

# Construct a synchronous RK2 integrator, and add the one physics package and
# boundary condtion.
integrator = CheapSynchronousRK2Integrator1d(db)
output('integrator')
integrator.appendPhysicsPackage(hydro)
integrator.appendBoundary(xbc0)
integrator.appendBoundary(xbc1)
output('integrator.havePhysicsPackage(hydro)')
output('integrator.haveBoundary(xbc0)')
output('integrator.haveBoundary(xbc1)')
output('integrator.valid')
#output('integrator.initialize()')
integrator.HsmoothMin = HsmoothMin
integrator.HsmoothMax = HsmoothMax
output('integrator.HsmoothMin')
output('integrator.HsmoothMax')
integrator.lastDt = dt
output('integrator.lastDt')
integrator.dtMin = dtMin
output('integrator.dtMin')
integrator.dtMax = dtMax
output('integrator.dtMax')
integrator.sumForMassDensity = sumForMassDensity
output('integrator.sumForMassDensity')

control = SpheralController(integrator, W, boundaryConditions=[xbc0, xbc1])
output('control')

# Smooth the initial conditions.
control.smoothState(smoothIters)

##################################################################################
# Plot the initial conditions
plotState(nodes1)

# Advance to the end time.
control.advance(goalTime, maxSteps)

# Plot the final state.
plotState(nodes1, 'blue')
