#-------------------------------------------------------------------------------
# Set up a pair of equal mass N-body points in a simple circular orbit of each
# other.
#-------------------------------------------------------------------------------
from Spheral3d import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from NodeHistory import *
from SpheralVisitDump import dumpPhysicsState
from math import *

print "3-D N-Body Gravity test -- two particle problem"

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(

    # Initial particle stuff
    r0 = 1.0,                      # (AU) Start stuff out at 1 AU from barycenter
    m0 = 1.0,                      # (earth masses) particle mass
    plummerLength = 1.0e-3,        # (AU) Plummer softening scale
    opening = 0.5,                 # (dimensionless, OctTreeGravity) opening parameter for tree walk
    fdt = 0.1,                     # (dimensionless, OctTreeGravity) timestep multiplier

    # Problem control
    steps = None,
    numOrbits = 25,                 # How many orbits do we want to follow?

    # Which N-body method should we use?
    nbody = OctTreeGravity,
    timeStepChoice = AccelerationRatio,
    integratorConstructor = CheapSynchronousRK2Integrator,

    # Output
    dataDir = "Two-Earth-Nbody",
    baseName = "2_particle_nbody",
    restoreCycle = None,
    restartStep = 100,
    numViz = 100,
    dtverbose = False,
    )

# Convert to MKS units.
AU = 149597870700.0  # m
Mearth = 5.9722e24   # kg
r0 *= AU
m0 *= Mearth
plummerLength *= AU

# Compute the velocity necessary for an elliptical orbit.
G = MKS().G
a = 2.0*r0
M = 2*m0
orbitTime = 2.0*pi*sqrt(a**3/(G*M))
v0 = 0.5*sqrt(0.5*M*G*(2.0/r0-1.0/a)) #vis-viva equation

# Miscellaneous problem control parameters.
dt = orbitTime / 100
goalTime = orbitTime * numOrbits
dtMin, dtMax = 0.01*dt, 100.0*dt
dtGrowth = 2.0
maxSteps = None
statsStep = 10
smoothIters = 0
vizTime = goalTime / numViz

restartDir = os.path.join(dataDir, "restarts")
visitDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, baseName + "_restart")

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
# For now we have set up a fluid node list, even though this is collisionless
# problem.  Fix at some point!
# In the meantime, set up the hydro objects this script isn't really going to
# need.
#-------------------------------------------------------------------------------
WT = TableKernel(BSplineKernel(), 1000)
eos = GammaLawGasMKS3d(gamma = 5.0/3.0, mu = 1.0)

#-------------------------------------------------------------------------------
# Make the NodeList, and set our initial properties.
#-------------------------------------------------------------------------------
nodes = makeFluidNodeList("nodes", eos,
                          numInternal = 2,
                          xmin = Vector(-100*r0, -100*r0, -100*r0),
                          xmax = Vector( 100*r0,  100*r0,  100*r0))
mass = nodes.mass()
pos = nodes.positions()
vel = nodes.velocity()

mass[0] = m0
mass[1] = m0

pos[0] = Vector(-r0, 0.0, 0.0)
pos[1] = Vector( r0, 0.0, 0.0)

vel[0] = Vector(0.0, -v0, 0.0)
vel[1] = Vector(0.0,  v0, 0.0)

# These are fluid variables we shouldn't need.  Just set them to valid values.
H = nodes.Hfield()
rho = nodes.massDensity()
H[0] = r0*SymTensor.one
H[1] = r0*SymTensor.one
rho[0] = 1.0
rho[1] = 1.0

#-------------------------------------------------------------------------------
# DataBase
#-------------------------------------------------------------------------------
db = DataBase()
db.appendNodeList(nodes)

#-------------------------------------------------------------------------------
# Gimme gravity.
#-------------------------------------------------------------------------------
if nbody is NBodyGravity:
    gravity = NBodyGravity(plummerSofteningLength = plummerLength,
                           maxDeltaVelocity = 1e-2*v0,
                           G = G)
elif nbody is OctTreeGravity:
    gravity = OctTreeGravity(G = G,
                             softeningLength = plummerLength,
                             opening = opening,
                             ftimestep = fdt,
                             timeStepChoice = timeStepChoice)
elif nbody is FractalGravity:
    gravity = FractalGravity(G = G,
                             xmin = Vector(-1.5*r0, -1.5*r0, -1.5*r0),
                             xmax = Vector( 1.5*r0,  1.5*r0,  1.5*r0),
                             periodic = False,
                             ngrid = 64,
                             nlevelmax = 1,
                             minHighParticles = 10,
                             padding = 0,
                             maxDeltaVelocity = 1e-2*v0)

#-------------------------------------------------------------------------------
# Construct a time integrator.
#-------------------------------------------------------------------------------
integrator = integratorConstructor(db)
integrator.appendPhysicsPackage(gravity)
integrator.lastDt = 1e3    # seconds
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.verbose = dtverbose

#-------------------------------------------------------------------------------
# Build the problem controller to follow the problem evolution.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            vizBaseName = os.path.join(visitDir, baseName),
                            vizTime = vizTime,
                            vizMethod = dumpPhysicsState)

#-------------------------------------------------------------------------------
# Build a diagnostic to maintain the history of our points.
#-------------------------------------------------------------------------------
def sampleMethod(nodes, indices):
    m = nodes.mass()
    pos = nodes.positions()
    vel = nodes.velocity()
    assert nodes.numInternalNodes == 2
    r = (pos[1] - pos[0]).magnitude()
    return (m[0], pos[0].x, pos[0].y, pos[0].z, vel[0].x, vel[0].y, vel[0].z, 
            m[1], pos[1].x, pos[1].y, pos[1].z, vel[1].x, vel[1].y, vel[1].z,
            0.5*(m[0]*vel[0].magnitude2() + m[1]*vel[1].magnitude2()),
            -2.0*G*m[0]*m[1]/r)

sampleNodes = [0, 1]  # We're going to sample both of our nodes!
history = NodeHistory(nodes, sampleNodes, sampleMethod,
                      os.path.join(dataDir, "node_history.txt"),
                      header = "# Orbit history of a 2 earth (no sun) system.",
                      labels = ("m1", "x1", "y1", "z1", "vx1", "vy1", "vz1",
                                "m1", "x1", "y1", "z1", "vx1", "vy1", "vz1",
                                "KE", "PE"))
control.appendPeriodicTimeWork(history.sample, vizTime)

#-------------------------------------------------------------------------------
# If we're restarting, read in the restart file.
#-------------------------------------------------------------------------------
if restoreCycle:
    control.loadRestartFile(restoreCycle)

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)
else:
    control.advance(goalTime)

# Plot the final state.
x1 = [stuff[1]/AU for stuff in history.sampleHistory]
y1 = [stuff[2]/AU for stuff in history.sampleHistory]
x2 = [stuff[8]/AU for stuff in history.sampleHistory]
y2 = [stuff[9]/AU for stuff in history.sampleHistory]

import SpheralGnuPlotUtilities
gdata1 = Gnuplot.Data(x1, y1,
                      with_ = 'linespoints',
                      title = 'Particle 1')
gdata2 = Gnuplot.Data(x2, y2,
                      with_ = 'linespoints',
                      title = 'Particle 2')
plot = generateNewGnuPlot()
plot.plot(gdata1)
plot.replot(gdata2)
plot('set size square; set xrange [-4.1:4.1]; set yrange [-2.1:2.1]')
plot.xlabel = 'x'
plot.ylabel = 'y'
plot.refresh()
