#ATS:DEM3dSPBC = test(        SELF, "--clearDirectories True  --checkError True  --checkConservation True --normalRestitutionCoefficient 1.0 --steps 100", label="DEM perfectly elastic 2 particle collision -- 3-D (serial)")

import os, sys, shutil, mpi
from math import *
from Spheral3d import *
from SpheralTestUtilities import *
from findLastRestart import *
from GenerateNodeDistribution3d import *
from GenerateDEMfromSPHGenerator import GenerateDEMfromSPHGenerator3d

sys.path.insert(0, '..')
from DEMConservationTracker import TrackConservation3d as TrackConservation

if mpi.procs > 1:
    from PeanoHilbertDistributeNodes import distributeNodes3d
else:
    from DistributeNodes import distributeNodes3d

title("DEM Restitution Coefficient Test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(vImpact = 1.0,                            # impact velocity
            omega0 = 0.1,                             # initial angular velocity it we're doing that
            g0 = 0.0,                                 # grav acceleration

            h0 = 1.00,                                # initial height above the solid bc plane
            radius = 0.95,                            # particle radius
            normalSpringConstant=10000.0,             # spring constant for LDS model
            normalRestitutionCoefficient=1.00,        # restitution coefficient to get damping const
            tangentialSpringConstant=2857.0,          # spring constant for LDS model
            tangentialRestitutionCoefficient=0.55,    # restitution coefficient to get damping const
            dynamicFriction = 1.0,                    # static friction coefficient sliding
            staticFriction = 1.0,                     # dynamic friction coefficient sliding
            rollingFriction = 1.05,                   # static friction coefficient for rolling
            torsionalFriction = 1.3,                  # static friction coefficient for torsion
            cohesiveTensileStrength = 0.0,            # units of pressure
            shapeFactor = 0.5,                        # in [0,1] shape factor from Zhang 2018, 0 - no torsion or rolling

            nPerh = 1.01,                             # we need this as an input for thing but doesn't affect DEM
            
            # integration
            IntegratorConstructor = VerletIntegrator, # Verlet one integrator to garenteee conservation
            stepsPerCollision = 100,                   # replaces CFL for DEM
            goalTime = 3.0,
            dt = 1e-8,
            dtMin = 1.0e-8, 
            dtMax = 0.1,
            dtGrowth = 2.0,
            steps = None,
            maxSteps = None,
            statsStep = 10,
            domainIndependent = False,
            rigorousBoundaries = False,
            dtverbose = False,
            
            # output control
            vizCycle = None,
            vizTime = 0.1,
            clearDirectories = False,
            restoreCycle = None,
            restartStep = 1000,
            redistributeStep = 500,
            dataDir = "dumps-DEM-particle-boundary-3d", 

             # ats parameters
            checkError = False,                    # turn on error checking for restitution coefficient
            boolCheckSlidingFriction=False,        # checks sliding friction reduces relative rotation
            boolCheckRollingFriction=False,        # checks rolling friction reduces relative rotation 
            boolCheckTorsionalFriction=False,      # checks torsional friction reduces relative rotation
            boolCheckTorsionalObjectivity=False,   # checks to make sure torsion is objective
            checkRestart = False,                  # turn on error checking for restartability
            checkConservation = False,             # turn on error checking for momentum conservation
            restitutionErrorThreshold = 0.02,      # relative error actual restitution vs nominal
            conservationErrorThreshold = 1e-15,    # relative error for momentum conservation
            torsionalObjectivityThreshold = 1e-10  # relative error bounds on torsion objectivity test
            )

#-------------------------------------------------------------------------------
# check for bad inputs
#-------------------------------------------------------------------------------
assert mpi.procs == 1 
assert nPerh >= 1
assert g0 <= 0.0
assert h0 > radius
assert shapeFactor <= 1.0 and shapeFactor >= 0.0
assert dynamicFriction >= 0.0
assert staticFriction >= 0.0
assert torsionalFriction >= 0.0
assert rollingFriction >= 0.0
assert cohesiveTensileStrength >= 0.0
assert sum([boolCheckSlidingFriction,
            boolCheckRollingFriction,
            boolCheckTorsionalFriction,
            boolCheckTorsionalObjectivity]) <= 1

if boolCheckSlidingFriction:
    shapeFactor = 0.0
    
#-------------------------------------------------------------------------------
# file things
#-------------------------------------------------------------------------------
testName = "DEM-SingleParticleBoundaryCollision-3d"

dataDir = os.path.join(dataDir,
                  "restitutionCoefficient=%s" % normalRestitutionCoefficient,
                  "boolCheckSlidingFriction=%s" % boolCheckSlidingFriction,
                  "boolCheckRollingFriction=%s" % boolCheckRollingFriction,
                  "boolCheckTorsionalFriction=%s" % boolCheckTorsionalFriction,
                  "boolCheckTorsionalObjectivity=%s" % boolCheckTorsionalObjectivity)

restartDir = os.path.join(dataDir, "restarts")
vizDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, testName)
vizBaseName = testName

if vizCycle is None and vizTime is None:
    vizBaseName=None

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
if mpi.rank == 0:
    if clearDirectories and os.path.exists(dataDir):
        shutil.rmtree(dataDir)
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
    if not os.path.exists(vizDir):
        os.makedirs(vizDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# If we're restarting, find the set of most recent restart files.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    restoreCycle = findLastRestart(restartBaseName)

#-------------------------------------------------------------------------------
# This doesn't really matter kernel filler for neighbor algo
#-------------------------------------------------------------------------------
WT = TableKernel(WendlandC2Kernel(), 1000)

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
units = CGuS()
nodes1 = makeDEMNodeList("nodeList1",
                          hmin = 1.0e-30,
                          hmax = 1.0e30,
                          hminratio = 100.0,
                          nPerh = nPerh,
                          kernelExtent = WT.kernelExtent)
nodeSet = [nodes1]
for nodes in nodeSet:
    output("nodes.name")
    output("nodes.hmin")
    output("nodes.hmax")
    output("nodes.hminratio")
    output("nodes.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Set the node properties. (gen 2 particles visit doesn't like just one)
#-------------------------------------------------------------------------------
generator0 = GenerateNodeDistribution3d(2, 1, 1,
                                        rho = 1.0,
                                        distributionType = "lattice",
                                        xmin = (-2.0,  -1.0, -1+h0),
                                        xmax = (2.0,  1.0, 1+h0),
                                        nNodePerh = nPerh)

generator1 = GenerateDEMfromSPHGenerator3d(WT,
                                           generator0,
                                           nPerh=nPerh)
distributeNodes3d((nodes1, generator1))
 
#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
for nodes in nodeSet:
    db.appendNodeList(nodes)
output("db.numNodeLists")
output("db.numDEMNodeLists")
output("db.numFluidNodeLists")


#-------------------------------------------------------------------------------
# DEM
#-------------------------------------------------------------------------------
dem = DEM(db,
          normalSpringConstant = normalSpringConstant,
          normalRestitutionCoefficient = normalRestitutionCoefficient,
          tangentialSpringConstant = tangentialSpringConstant,
          tangentialRestitutionCoefficient = tangentialRestitutionCoefficient,
          dynamicFrictionCoefficient = dynamicFriction,
          staticFrictionCoefficient = staticFriction,
          rollingFrictionCoefficient = rollingFriction,
          torsionalFrictionCoefficient = torsionalFriction,
          cohesiveTensileStrength =cohesiveTensileStrength,
          shapeFactor = shapeFactor,
          stepsPerCollision = stepsPerCollision)

packages = [dem]



solidWall = PlanarWall(Vector(0.0, 0.0, 0.0), Vector(  0.0, 0.0, 1.0))
dem.appendSolidBoundary(solidWall)

#-------------------------------------------------------------------------------
# PhysicsPackage : gravity
#-------------------------------------------------------------------------------
gravity = ConstantAcceleration(a0 = Vector(0.0,0.0,g0),
                               nodeList = nodes1)
packages += [gravity]


#-------------------------------------------------------------------------------
# initial conditions
#-------------------------------------------------------------------------------

velocity = nodes1.velocity()
velocity[0] = Vector(0.0,0.0,-vImpact)
velocity[1] = Vector(0.0,0.0,-vImpact)

particleRadius = nodes1.particleRadius()
particleRadius[0] = radius
particleRadius[1] = radius

omega = dem.omega
omega[0][0] = Vector(0.0,0.0, omega0)
omega[0][1] = Vector(0.0,0.0,-omega0)

if boolCheckSlidingFriction:
    omega[0][0] = Vector(0.0,0.0,omega0)
    omega[0][1] = Vector(0.0,0.0,omega0)
elif boolCheckRollingFriction:  
    omega[0][0] = Vector(0.0,0.0, omega0)
    omega[0][1] = Vector(0.0,0.0,-omega0)
elif boolCheckTorsionalFriction:
    omega[0][0] = Vector( omega0,0.0,0.0)
    omega[0][1] = Vector(-omega0,0.0,0.0)
elif boolCheckTorsionalObjectivity:
    omega[0][0] = Vector( omega0,0.0,0.0)
    omega[0][1] = Vector( omega0,0.0,0.0)
    

#-------------------------------------------------------------------------------
# Construct a time integrator, and add the physics packages.
#-------------------------------------------------------------------------------

integrator = IntegratorConstructor(db)
for p in packages:
    integrator.appendPhysicsPackage(p)
integrator.lastDt = dt
integrator.dtMin = dtMin
integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.domainDecompositionIndependent = domainIndependent
integrator.verbose = dtverbose
integrator.rigorousBoundaries = rigorousBoundaries

integrator.cullGhostNodes = False

output("integrator")
output("integrator.havePhysicsPackage(dem)")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.domainDecompositionIndependent")
output("integrator.rigorousBoundaries")
output("integrator.verbose")


#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
from SpheralPointmeshSiloDump import dumpPhysicsState
control = SpheralController(integrator, WT,
                            iterateInitialH = False,
                            initializeDerivatives = True,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            vizBaseName = vizBaseName,
                            vizMethod=dumpPhysicsState,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime)
output("control")

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------

if not steps is None:
    if checkRestart:
        control.setRestartBaseName(restartBaseName + "_CHECK")
    control.step(steps)
else:
    control.advance(goalTime, maxSteps)

#-------------------------------------------------------------------------------
# Great success?
#-------------------------------------------------------------------------------
if checkRestart:
# check reproducibility when restarted
#-------------------------------------------------------------
    control.setRestartBaseName(restartBaseName)
    state0 = State(db, integrator.physicsPackages())
    state0.copyState()
    control.loadRestartFile(control.totalSteps)
    state1 = State(db, integrator.physicsPackages())
    if not state1 == state0:
        raise ValueError, "The restarted state does not match!"
    else:
        print "Restart check PASSED."

if checkError:
# check our restitution coefficient is correct
#-------------------------------------------------------------
    vijPostImpact = velocity[0].x - velocity[1].x
    vijPreImpact = 2.0*vImpact
    restitutionEff = vijPostImpact/vijPreImpact
    restitutionError = abs(restitutionEff + normalRestitutionCoefficient)/normalRestitutionCoefficient
    if  restitutionError > restitutionErrorThreshold:
        raise ValueError, ("relative restitution coefficient error, %g, exceeds bounds" % restitutionError)


if boolCheckSlidingFriction or boolCheckRollingFriction or boolCheckTorsionalFriction:
# check for non-physical behavior
#-------------------------------------------------------------
    if omega[0][0].magnitude()+omega[0][1].magnitude() > 2*omega0:
        raise ValueError, "particles are rotating faster post-collision"

if boolCheckTorsionalObjectivity:
# to satify objectivity omega (along axis) should not change when equal
#-------------------------------------------------------------
    omegaError = (2*omega0 - omega[0][0][0] - omega[0][1][0]) / (2*omega0)
    if omegaError > torsionalObjectivityThreshold:
        raise ValueError, ("torsional objectivity failure with relative angular velocity error, %g, exceeds bounds" % omegaError)
