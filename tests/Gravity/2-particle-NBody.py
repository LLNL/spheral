from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from math import *

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(
    # Particle positions.
    r1 = Vector3d(-1.0, 0.0, 0.0),
    r2 = Vector3d(1.0, 0.0, 0.0),

    # Particle velocities.  We want a circular orbit.
    v = sqrt(6.672e-8/4.0),

    # Particle H tensors.
    H1 = SymTensor3d(100.0, 0.0, 0.0,
                     0.0, 100.0, 0.0,
                     0.0,   0.0, 100.0),
    H2 = SymTensor3d(100.0, 0.0, 0.0,
                     0.0, 100.0, 0.0,
                     0.0,   0.0, 100.0),

    gamma = 5.0/3.0,
    mu = 1.0,
    Cl = 1.0,
    Cq = 0.75,
    epsilon2 = 1e-8,
    HsmoothMin = 0.004,
    HsmoothMax = 0.5,
    HsmoothFraction = 0.0,
    cfl = 0.5,

    neighborSearchType = Neighbor3d.NeighborSearchType.GatherScatter,
    numGridLevels = 16,
    topGridCellSize = 0.5,
    origin = Vector3d(0.0, 0.0, 0.0),
    )

v1 = Vector3d(0.0, -v, 0.0)
v2 = Vector3d(0.0, v, 0.0)
orbitTime = 2*pi/v
dt = orbitTime / 90
goalTime = orbitTime * 2
dtMin, dtMax = dt, dt
dtGrowth = 2.0
maxSteps = None
statsStep = 10
smoothIters = 0

restoreCycle = 0
restartStep = 40
restartBaseName = "2-particle-Gravity"

#-------------------------------------------------------------------------------
title('3-D N-Body Gravity test -- two particle problem')

# Set the table kernel.
WT = TableKernel3d(BSplineKernel3d(), 1000)
output('WT')

# Make an equation of state.
eos = GammaLawGasCGS3d(gamma, mu)

# Create an empty NodeList
nodes1 = SphNodeList3d("nodes1", eos, WT, WT, 2)
output('nodes1')
nodes1.HsmoothFraction = HsmoothFraction

# Set particle positions, velocities, and H tensors.
nodes1.positions()[0] = r1
nodes1.velocity()[0] = v1
nodes1.Hfield()[0] = H1
nodes1.positions()[1] = r2
nodes1.velocity()[1] = v2
nodes1.Hfield()[1] = H2

# Manually initialize the mass density for the nodes.
nodes1.massDensity()[0] = 1.0
nodes1.mass()[0] = 1.0
nodes1.massDensity()[1] = 1.0
nodes1.mass()[1] = 1.0

# Construct the neighbor object and associate it with the node list.
neighbor1 = NestedGridNeighbor3d(nodes1,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 2) # Kernel extent (sans Kernel... :-P)
nodes1.registerNeighbor(neighbor1)

# Note: No boundary conditions.

# Construct a DataBase to hold our node list
db = DataBase3d()
output('db')
output('db.appendNodeList(nodes1)')
output('db.numNodeLists')
output('db.numFluidNodeLists')

# Gimme gravity.
plummerLength = 0.001
dvMax = 2.0
gravity = NBodyGravity3d(plummerLength, dvMax)

# Construct a synchronous RK2 integrator, and add the physics package.
integrator = PredictorCorrectorIntegrator3d(db)
output('integrator')
integrator.appendPhysicsPackage(gravity)
output('integrator.havePhysicsPackage(gravity)')
output('integrator.valid()')
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
    position1.append(Vector3d(nodes1.positions()[0]))
    position2.append(Vector3d(nodes1.positions()[1]))

# Plot the final state.
import Gnuplot
gdata1 = Gnuplot.Data([r.x for r in position1],
                      [r.y for r in position1],
                      with = 'linespoints',
                      title = 'Particle 1')
gdata2 = Gnuplot.Data([r.x for r in position2],
                      [r.y for r in position2],
                      with = 'linespoints',
                      title = 'Particle 2')
plot = Gnuplot.Gnuplot()
plot.plot(gdata1)
plot.replot(gdata2)
plot('set size square')
plot.xlabel = 'x'
plot.ylabel = 'y'
plot.refresh()

