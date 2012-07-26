#-------------------------------------------------------------------------------
# Create a spherical distribution of collisionless points, which will of course 
# promptly collapse under their own self-gravity.
#-------------------------------------------------------------------------------
from Spheral3d import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from NodeHistory import *
from SpheralVisitDump import dumpPhysicsState
from GenerateNodeDistribution3d import *
from VoronoiDistributeNodes import distributeNodes3d
from math import *

print "3-D N-Body Gravity test -- collisionless sphere."

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(

    # Initial particle stuff
    nr = 50,                       # radial number of particles to seed
    r0 = 1.0,                      # (AU) initial radius of the sphere
    M0 = 1.0,                      # (earth masses) total mass of the sphere
    plummerLength = 1.0e-2,        # (AU) Plummer softening scale
    opening = 1.0,                 # (dimensionless, OctTreeGravity) opening parameter for tree walk
    fdt = 0.1,                     # (dimensionless, OctTreeGravity) timestep multiplier

    # Problem control
    steps = None,                  # Optionally advance a fixed number of steps
    numCollapseTimes = 1.0,        # What time to advance to in units of the collapse time for the sphere

    # Which N-body method should we use?
    nbody = OctTreeGravity,

    # Output
    dataDir = "Collisionless_Sphere_Collapse",
    baseNameRoot = "sphere_collapse_%i",
    restoreCycle = None,
    restartStep = 100,
    numViz = 200,
    )

# Convert to MKS units.
AU = 149597870700.0  # m
Mearth = 5.9722e24   # kg
r0 *= AU
M0 *= Mearth
plummerLength *= AU

# Compute the analytically expected collapse time.
G = MKS().G
rho0 = M0/(4.0/3.0*pi*r0*r0*r0)
tdyn = sqrt(3.0*pi/(16*G*rho0))
collapseTime = tdyn/sqrt(2.0)

# Miscellaneous problem control parameters.
dt = collapseTime / 100
goalTime = collapseTime * numCollapseTimes
dtMin, dtMax = 0.1*dt, 100.0*dt
dtGrowth = 2.0
maxSteps = None
statsStep = 10
smoothIters = 0
vizTime = goalTime / numViz

baseName = baseNameRoot % nr
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
                          topGridCellSize = 100*r0)

if restoreCycle is None:
    generator = GenerateNodeDistribution3d(2*nr, 2*nr, 2*nr, rho0,
                                           distributionType = "lattice",
                                           xmin = (-r0, -r0, -r0),
                                           xmax = ( r0,  r0,  r0),
                                           rmin = 0.0,
                                           rmax = r0)
    distributeNodes3d((nodes, generator))
    output("mpi.reduce(nodes.numInternalNodes, mpi.MIN)")
    output("mpi.reduce(nodes.numInternalNodes, mpi.MAX)")
    output("mpi.reduce(nodes.numInternalNodes, mpi.SUM)")

    # Renormalize the node masses to get our desired total mass.
    mass = nodes.mass()
    msum = mpi.allreduce(sum(nodes.mass().internalValues()), mpi.SUM)
    fmass = M0/msum
    print "Renormalizing masses by ", fmass
    for i in xrange(nodes.numInternalNodes):
        mass[i] *= fmass

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
                           maxDeltaVelocity = 1e-2*r0/tdyn,
                           G = G)
# elif nbody is FractalGravity:
#     gravity = FractalGravity(G = G,
#                              xmin = Vector(-1.5*r0, -1.5*r0, -1.5*r0),
#                              xmax = Vector( 1.5*r0,  1.5*r0,  1.5*r0),
#                              periodic = False,
#                              ngrid = 64,
#                              nlevelmax = 1,
#                              minHighParticles = 10,
#                              padding = 0,
#                              maxDeltaVelocity = 1e-2*v0)
elif nbody is OctTreeGravity:
    gravity = OctTreeGravity(G = G,
                             softeningLength = plummerLength,
                             opening = opening,
                             ftimestep = fdt)

#-------------------------------------------------------------------------------
# Construct a time integrator.
#-------------------------------------------------------------------------------
integrator = SynchronousRK2Integrator(db)
integrator.appendPhysicsPackage(gravity)
integrator.lastDt = dtMin
integrator.dtMin = dtMin
integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth

#-------------------------------------------------------------------------------
# Build the problem controller to follow the problem evolution.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            vizBaseName = baseName,
                            vizDir = visitDir,
                            vizTime = vizTime,
                            vizMethod = dumpPhysicsState)

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
    print "Advancing to %g sec = %g years" % (goalTime, goalTime/(365.24*24*3600))
    control.advance(goalTime)
