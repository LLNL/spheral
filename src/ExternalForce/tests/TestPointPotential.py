from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from math import *

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
rmin, rmax = (0.0,0.0,0.0), (1.0,1.0,1.0)

# Point mass exerting external force.
G  = 1.0
M0 = 1.0
Rc = 1e-2
R0 = (0.0, 0.0)
deltaPotentialFraction = 0.01

# Particle positions and masses.
m1 = 1.0e-10
r1 = Vector2d(-0.005, 0.0)

# Particle velocities.  We want a circular orbit.
#v = sqrt(2.0*G*M0/sqrt(r1.magnitude2() + Rc*Rc))
#v = sqrt(G*M0/sqrt(r1.magnitude2() + Rc*Rc))
v = sqrt(G*M0*r1.magnitude2()* (r1.magnitude2() + Rc*Rc)**(-1.5))
v1 = (0.0, -v)

# Particle H tensors.
H1 = SymTensor2d(100.0, 0.0,
                 0.0, 100.0)

gamma = 5.0/3.0
mu = 1.0
Cl, Cq = 1.0, 0.75
epsilon2 = 1e-8
HsmoothMin, HsmoothMax = 0.004, 0.5 #0.0001, 0.1
HsmoothFraction = 0.0 #0.1
cfl = 0.5

neighborSearchType = 3 # GatherScatter
numGridLevels = 16
topGridCellSize = 0.5
origin = Vector2d(0.0, 0.0)

orbitTime = 2*pi*r.magnitude()/v
dt = orbitTime / 100
goalTime = orbitTime * 10
dtMin, dtMax = 1.0e-5*orbitTime, orbitTime
dtGrowth = 2.0
maxSteps = None
statsStep = 10
smoothIters = 0

restoreCycle = 0
restartStep = 1000000
restartBaseName = "PointPotentialOrbit"
restartFileConstructor = FlatFileIO2d

#-------------------------------------------------------------------------------
title('2-D Point Potential Gravity test -- single particle orbit.')

eos = GammaLawGasCGS2d(gamma, mu)

# Create an empty NodeList
nodes1 = SphNodeList2d(eos, 1)
output('nodes1')
nodes1.HsmoothFraction = HsmoothFraction

# Set particle positions, velocities, and H tensors.
nodes1.positions[0] = r1
nodes1.velocity[0] = v1
nodes1.Hfield[0] = H1

# Set node specific thermal energies
nodes1.specificThermalEnergy[:] = [0]*nodes1.numNodes

# Manually initialize the mass density for the nodes.
nodes1.massDensity[0] = 1.0
nodes1.mass[0] = 1.0
nodes1.massDensity[1] = 1.0
nodes1.mass[1] = 1.0

# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
W = BSplineKernel2d()
output('W')

# Set the table kernel.
WT = TableKernel2d()
WT.setTableData(W, 100)

# Construct the neighbor object and associate it with the node list.
neighborTimer = SpheralTimer('Neighbor initialization.')
neighborTimer.start()
neighbor1 = NestedGridNeighbor2d(nodes1,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 2) # Kernel extent (sans Kernel... :-P)
nodes1.neighbor = neighbor1
neighborTimer.stop()
neighborTimer.printStatus()

# Note: No boundary conditions.

# Construct a DataBase to hold our node list
db = DataBase2d()
output('db')
output('db.appendNodeList(nodes1)')
output('db.numNodeLists')
output('db.numFluidNodeLists')

# Gimme gravity.
gravity = PointPotential2d(G, M0, Rc, R0)
gravity.deltaPotentialFraction = deltaPotentialFraction
output("gravity.G")
output("gravity.mass")
output("gravity.coreRadius")
output("gravity.origin")
output("gravity.deltaPotentialFraction")

# Construct an integrator, and add the physics package.
#integrator = PredictorCorrectorIntegrator2d(db)
integrator = CheapSynchronousRK2Integrator2d(db)
output('integrator')
integrator.setKernel(WT)
integrator.appendPhysicsPackage(gravity)
output('integrator.havePhysicsPackage(gravity)')
output('integrator.valid()')
integrator.HsmoothMin = HsmoothMin
integrator.HsmoothMax = HsmoothMax
output('integrator.HsmoothMin')
output('integrator.HsmoothMax')
integrator.lastDt = dt
output('integrator.lastDt')
if dtMin:
    integrator.dtMin = dtMin
    output('integrator.dtMin')
if dtMax:
    integrator.dtMax = dtMax
    output('integrator.dtMax')
integrator.dtGrowth = dtGrowth
output('integrator.dtGrowth')

# Tell the integrator to integrate to update the mass density instead of 
# solving.
# FIXME: Uhhhhh... we don't need no steenking kernel.
integrator.sumForMassDensity = 2

restartObjects = ['nodes1', 'integrator', 'control']
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            printStep = 100,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restartFileConstructor = restartFileConstructor,
                            restartObjects = restartObjects)
output('control')

# Smooth the initial conditions.
if restoreCycle:
    control.loadRestartFile(restoreCycle)
else:
    control.smoothState(smoothIters)

##################################################################################
# Advance to the end time.
position1 = []
while control.time() < goalTime:
    control.advance(goalTime, 1000)
    position1.append(nodes1.positions[0][:])
xPosition1 = array([x[0] for x in position1])
yPosition1 = array([x[1] for x in position1])

# Plot the final state.
import Gnuplot
gdata = Gnuplot.Data(xPosition1, yPosition1, with='points')
plot = Gnuplot.Gnuplot()
plot.title = 'particle 1 position'
plot.xlabel = 'x'
plot.ylabel = 'y'
plot.plot(gdata)
plot('set size square')
plot.refresh()

