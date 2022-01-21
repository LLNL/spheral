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

print "Klemperer Rosette problem with 4 bodies"


'''
            m2
          /    \
        m1 ---- m1
          \    /
            m2
            
'''

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(

    # Initial particle stuff
    r1  = 1.0,                      # (AU) Start stuff out at 1 AU from barycenter
    m1  = 1.0,                      # (earth masses) particle mass
    ya  = 1.5,                      # ratio of r1/r2 - smaller mass (m1) is always farther out
    plummerLength = 1.0e-4,        # (AU) Plummer softening scale
    opening = 0.5,                 # (dimensionless, OctTreeGravity) opening parameter for tree walk
    fdt = 0.1,                     # (dimensionless, OctTreeGravity) timestep multiplier

    # Problem control
    steps = None,
    numOrbits = 2,                 # How many orbits do we want to follow?

    # Which N-body method should we use?
    nbody = OctTreeGravity,
    timeStepChoice = AccelerationRatio,
    integratorConstructor = CheapSynchronousRK2Integrator,

    # Output
    dataDir = "Klemperer-Nbody",
    baseName = "klemperer_nbody",
    restoreCycle = None,
    restartStep = 100,
    numViz = 100,
    dtverbose = False,
    )

# Convert to MKS units.
AU = 149597870700.0  # m
Mearth = 5.9722e24   # kg
r1 *= AU
r2 = r1/ya
m1 *= Mearth
mu = (8.0-pow(1.0+ya*ya,3.0/2.0))/(8.0-pow(1.0+1.0/(ya*ya),3.0/2.0))
m2 = m1/mu
plummerLength *= AU

# Compute the velocity necessary for a circular orbit.
G = MKS().G
omega = sqrt(G/r1*(2.0*m2*r1/pow(r1*r1+r2*r2,3.0/2.0)+m1/(4.0*r1*r1)))
orbitTime = 2.0*pi/omega
v1 = omega*r1
v2 = omega*r2

# Miscellaneous problem control parameters.
dt = orbitTime / 90
goalTime = orbitTime * numOrbits
dtMin, dtMax = 1e3, 1000.0*dt
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
                          numInternal = 4,
                          xmin = Vector(-100*r1, -100*r1, -100*r1),
                          xmax = Vector( 100*r1,  100*r1,  100*r1))
mass = nodes.mass()
pos = nodes.positions()
vel = nodes.velocity()

mass[0] = m1
mass[1] = m1
mass[2] = m2
mass[3] = m2

pos[0] = Vector(-r1, 0.0, 0.0)
pos[1] = Vector( r1, 0.0, 0.0)
pos[2] = Vector(0.0, -r2, 0.0)
pos[3] = Vector(0.0,  r2, 0.0)

vel[0] = Vector(0.0,  v1, 0.0)
vel[1] = Vector(0.0, -v1, 0.0)
vel[2] = Vector(-v2, 0.0, 0.0)
vel[3] = Vector( v2, 0.0, 0.0)

# These are fluid variables we shouldn't need.  Just set them to valid values.
H = nodes.Hfield()
rho = nodes.massDensity()
H[0] = r1*SymTensor.one
H[1] = r1*SymTensor.one
H[2] = r1*SymTensor.one
H[3] = r1*SymTensor.one
rho[0] = 1.0
rho[1] = 1.0
rho[2] = 1.0
rho[3] = 1.0

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
                           maxDeltaVelocity = 0.5e-1*v1,
                           G = G)
elif nbody is OctTreeGravity:
    gravity = OctTreeGravity(G = G,
                             softeningLength = plummerLength,
                             opening = opening,
                             ftimestep = fdt,
                             timeStepChoice = timeStepChoice)
elif nbody is FractalGravity:
    gravity = FractalGravity(G = G,
                             xmin = Vector(-1.5*r1, -1.5*r1, -1.5*r1),
                             xmax = Vector( 1.5*r1,  1.5*r1,  1.5*r1),
                             periodic = False,
                             ngrid = 64,
                             nlevelmax = 1,
                             minHighParticles = 10,
                             padding = 0,
                             maxDeltaVelocity = 1e-2*v1)

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
    H = nodes.Hfield()
    assert nodes.numInternalNodes == 4
    r = []
    ke = []
    for j in xrange(4):
        r.append(pos[j].magnitude())
        ke.append(0.5*m[j]*vel[j].magnitude()*vel[j].magnitude())
    return (m[0], pos[0].x, pos[0].y, pos[0].z, vel[0].x, vel[0].y, vel[0].z,r[0],
            m[1], pos[1].x, pos[1].y, pos[1].z, vel[1].x, vel[1].y, vel[1].z,r[1],
            m[2], pos[2].x, pos[2].y, pos[2].z, vel[2].x, vel[2].y, vel[2].z,r[2],
            m[3], pos[3].x, pos[3].y, pos[3].z, vel[3].x, vel[3].y, vel[3].z,r[3],ke[0]+ke[1]+ke[2]+ke[3])

sampleNodes = [0, 1]  # We're going to sample both of our nodes!
history = NodeHistory(nodes, sampleNodes, sampleMethod,
                      os.path.join(dataDir, "node_history.txt"),
                      header = "# Orbit history of a 4 body system.",
                      labels = ("m1", "x1", "y1", "z1", "vx1", "vy1", "vz1", "r1"
                                "m2", "x2", "y2", "z2", "vx2", "vy2", "vz2", "r2"
                                "m3", "x3", "y3", "z3", "vx3", "vy3", "vz3", "r3"
                                "m4", "x4", "y4", "z4", "vx4", "vy4", "vz4", "r4", "KE"))
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
x2 = [stuff[9]/AU for stuff in history.sampleHistory]
y2 = [stuff[10]/AU for stuff in history.sampleHistory]
x3 = [stuff[17]/AU for stuff in history.sampleHistory]
y3 = [stuff[18]/AU for stuff in history.sampleHistory]
x4 = [stuff[25]/AU for stuff in history.sampleHistory]
y4 = [stuff[26]/AU for stuff in history.sampleHistory]

import SpheralGnuPlotUtilities
gdata1 = Gnuplot.Data(x1, y1,
                      with_ = 'linespoints',
                      title = 'Particle 1')
gdata2 = Gnuplot.Data(x2, y2,
                      with_ = 'linespoints',
                      title = 'Particle 2')
gdata3 = Gnuplot.Data(x3, y3,
                      with_ = 'linespoints',
                      title = 'Particle 3')
gdata4 = Gnuplot.Data(x4, y4,
                      with_ = 'linespoints',
                      title = 'Particle 4')
plot = generateNewGnuPlot()
plot.plot(gdata1)
plot.replot(gdata2)
plot.replot(gdata3)
plot.replot(gdata4)
plot('set size square; set xrange [-1.1:1.1]; set yrange [-1.1:1.1]')
plot.xlabel = 'x'
plot.ylabel = 'y'
plot.refresh()
