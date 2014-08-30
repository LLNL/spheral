from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from math import *

from Gadget import *

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
rmin, rmax = (0.0,0.0,0.0), (1.0,1.0,1.0)

# Particle positions.
r1 = (-1.0, 0.0, 0.0)
r2 = (1.0, 0.0, 0.0)

# Particle velocities.  We want a circular orbit.
v = sqrt(6.672e-8/4.0)
v1 = (0.0, -v, 0.0)
v2 = (0.0, v, 0.0)

# Particle H tensors.
H1 = SymTensor3d(100.0, 0.0, 0.0,
                   0.0, 100.0, 0.0,
                   0.0,   0.0, 100.0)
H2 = SymTensor3d(100.0, 0.0, 0.0,
                   0.0, 100.0, 0.0,
                   0.0,   0.0, 100.0)

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
origin = Vector3d(0.0, 0.0, 0.0)

orbitTime = 2*pi/v
dt = orbitTime / 90
goalTime = orbitTime * 2
dtMin, dtMax = dt, dt
#dtMin, dtMax = 1.0e-5, 0.1
dtGrowth = 2.0
maxSteps = None
statsStep = 10
smoothIters = 0

restoreCycle = 0
restartStep = 40
restartBaseName = "2-particle-Gravity"
restartFileConstructor = FlatFileIO3d

#-------------------------------------------------------------------------------
title('3-D Gadget Gravity test -- two particle problem')

eos = GammaLawGasCGS3d(gamma, mu)

# Create an empty NodeList
nodes1 = SphNodeList3d(2, eos)
output('nodes1')
nodes1.HsmoothFraction = HsmoothFraction

# Set particle positions, velocities, and H tensors.
nodes1.positions[0] = r1
nodes1.velocity[0] = v1
nodes1.Hfield[0] = H1
nodes1.positions[1] = r2
nodes1.velocity[1] = v2
nodes1.Hfield[1] = H2

# Set node specific thermal energies
nodes1.specificThermalEnergy[:] = [0]*nodes1.numNodes

# Manually initialize the mass density for the nodes.
nodes1.massDensity[0] = 1.0
nodes1.mass[0] = 1.0
nodes1.massDensity[1] = 1.0
nodes1.mass[1] = 1.0

# Set the table kernel.
WT = TableKernel3d(BSplineKernel3d(), 100)
output('WT')

# Construct the neighbor object and associate it with the node list.
neighborTimer = SpheralTimer('Neighbor initialization.')
neighborTimer.start()
neighbor1 = NestedGridNeighbor3d(nodes1,
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
db = DataBase3d()
output('db')
output('db.appendNodeList(nodes1)')
output('db.numNodeLists')
output('db.numFluidNodeLists')

# Gimme gravity.
gravity = GadgetGravityForce()

# Construct a synchronous RK2 integrator, and add the physics package.
integrator = CheapSynchronousRK2Integrator3d(db)
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
integrator.sumForMassDensity = 2;

control = SpheralController(integrator, WT,
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
position1 = []
position2 = []
for t in xrange(0, 100):
    control.advance(goalTime, 1)
    position1.append(nodes1.positions[0][:])
    position2.append(nodes1.positions[1][:])
xPosition1 = array([x[0] for x in position1])
yPosition1 = array([x[1] for x in position1])

# Plot the final state.
import Gnuplot
gdata = Gnuplot.Data(xPosition1, yPosition1, with='points')
plot = Gnuplot.Gnuplot()
plot('set size square')
plot.title = 'particle 1 position'
plot.xlabel = 'x'
plot.ylabel = 'y'
plot.plot(gdata)

