from Spheral import *
from SpheralTestUtilities import *
from SpheralController import *

################################################################################
# Generic problem parameters

# Initial conditions for the particle.
m0 = 1.0
r0 = 0.025
v0 = 1.0
rho0 = 1.0
eps0 = 0.0
h0 = 1.0/0.02

# Planes to define the periodic boundary condition.
xPlane0 = Plane1d((0.0), (1.0))
xPlane1 = Plane1d((1.0), (-1.0))

# Equation of state.
gamma = 5.0/3.0
mu = 1.0

# Q
Cl, Cq = 0.0, 0.0
epsilon2 = 1e-8

HsmoothMin, HsmoothMax = 0.02, 0.02
cfl = 0.5

# Neighbor search parameters.
neighborSearchType = 3 # GatherScatter
numGridLevels = 10
topGridCellSize = 0.25
origin = Vector1d(0.0)

# Miscellaneous parameters.
goalTime = 5.0
dt = 0.05
dtMin, dtMax = dt, dt
maxSteps = None
statsStep = 1
smoothIters = 0
sumForMassDensity = 2

################################################################################
title('1-D periodic boundary test.')

eos = GammaLawGasMKS1d(gamma, mu)

# Create an SphNodeList with 1 node.
nodes = SphNodeList1d(eos, 1)

# Set the initial conditions for this node.
nodes.mass[0] = m0
nodes.positions[0].x = r0
nodes.velocity[0].x = v0
nodes.massDensity[0] = rho0
nodes.specificThermalEnergy[0] = eps0
nodes.Hfield[0].xx = h0

# Create our interpolation kernels.
W = BSplineKernel1d()
WT = TableKernel1d()
WT.setTableData(W, 1000)
kernelExtent = W.kernelExtent

# Construct the neighbor object and associate it with the node list.
neighbor = NestedGridNeighbor1d(nodes,
                                neighborSearchType,
                                numGridLevels,
                                topGridCellSize,
                                origin,
                                kernelExtent)
nodes.neighbor = neighbor

# Create the periodic boundary condition.
xbc0 = PeriodicBoundary1d(xPlane0, xPlane1)

# Construct a DataBase to hold our node list
db = DataBase1d()
db.appendNodeList(nodes)

# Construct a standard Monaghan-Gingold artificial viscosity.
qMG = MonaghanGingoldViscosity1d(Cl, Cq)
qMG.epsilon2 = epsilon2

# Construct the hydro physics object.
hydro = Hydro1d(qMG)
hydro.kernel = WT
hydro.PiKernel = WT
hydro.cfl = cfl

# Construct a synchronous RK2 integrator, and add the one physics package and
# boundary condtion.
integrator = CheapSynchronousRK2Integrator1d(db)
integrator.appendPhysicsPackage(hydro)
integrator.HsmoothMin = HsmoothMin
integrator.HsmoothMax = HsmoothMax
integrator.lastDt = dt
integrator.kernel = WT
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.sumForMassDensity = sumForMassDensity

control = SpheralController(integrator, WT,
                            boundaryConditions = [xbc0],
                            statsStep = statsStep)

##################################################################################
# Advance to the end time one step at a time, taking snapshots of the particle
# position as you go.
xInternal = []
xGhost = []
while control.time() < goalTime:
    control.advance(goalTime, 1)
    output("control.time(), nodes.numInternalNodes, nodes.numGhostNodes")
    xInternal.append((control.time(), nodes.positions[0].x))
    if nodes.numGhostNodes > 0:
        xGhost.append((control.time(), nodes.positions[1].x))

# Make Numeric arrays out of the history of our positions.
import numpy
tInternalArray = numpy.array([x[0] for x in xInternal])
xInternalArray = numpy.array([x[1] for x in xInternal])
tGhostArray = numpy.array([x[0] for x in xGhost])
xGhostArray = numpy.array([x[1] for x in xGhost])

# Now we can plot them with Gnuplot.
import Gnuplot
dInternal = Gnuplot.Data(tInternalArray, xInternalArray,
                         with = 'points',
                         title = 'Internal node')
dGhost = Gnuplot.Data(tGhostArray, xGhostArray,
                      with = 'points',
                      title = 'Ghost node')
plot = Gnuplot.Gnuplot()
plot.replot(dInternal)
plot.replot(dGhost)

