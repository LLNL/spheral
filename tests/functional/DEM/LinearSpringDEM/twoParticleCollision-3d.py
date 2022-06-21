#ATS:DEM3d0 = test(        SELF, "--clearDirectories True  --checkError True  --checkConservation True --restitutionCoefficient=1.0 --steps 100", label="DEM perfectly elastic 2 particle collision -- 3-D (serial)")
#ATS:DEM3d1 = test(        SELF, "--clearDirectories True  --checkError True  --checkConservation True --restitutionCoefficient=0.8 --steps 100", label="DEM perfectly inelastic 2 particle collision -- 3-D (serial)")
#ATS:DEM3d2 = test(        SELF, "--clearDirectories True  --checkError True --boolCheckSlidingFriction True --checkConservation True --restitutionCoefficient=0.8 --steps 100", label="DEM inelastic 2 particle collision - sliding friction -- 3-D (serial)")
#ATS:DEM3d3 = test(        SELF, "--clearDirectories True  --checkError True --boolCheckRollingFriction True --checkConservation True --restitutionCoefficient=0.8 --steps 100", label="DEM inelastic 2  particle collision - rolling friction -- 3-D (serial)")
#ATS:DEM3d4 = test(        SELF, "--clearDirectories True  --checkError True --boolCheckTorsionalFriction True --checkConservation True --restitutionCoefficient=0.8 --steps 100", label="DEM inelastic 2  particle collision - torsional friction -- 3-D (serial)")
#ATS:DEM3d5 = test(        SELF, "--clearDirectories True  --checkError True --boolCheckTorsionalObjectivity True --checkConservation True --restitutionCoefficient=0.8 --steps 100", label="DEM inelastic 2  particle collision - torsional objectivity -- 3-D (serial)")

import os, sys, shutil, mpi
from math import *
from Spheral3d import *
from SpheralTestUtilities import *
from findLastRestart import *
from GenerateNodeDistribution3d import *
from DEMConservationTracker import TrackConservation3d as TrackConservation

if mpi.procs > 1:
    from PeanoHilbertDistributeNodes import distributeNodes3d
else:
    from DistributeNodes import distributeNodes3d

title("DEM Restitution Coefficient Test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(vImpact = 1.0,                 # impact velocity
            omega0 = 0.1,                  # initial angular velocity it we're doing that

            radius = 0.25,                 # particle radius
            normalSpringConstant=10000.0,           # spring constant for LDS model
            normalRestitutionCoefficient=0.55,      # restitution coefficient to get damping const
            tangentialSpringConstant=2857.0,        # spring constant for LDS model
            tangentialRestitutionCoefficient=0.55,  # restitution coefficient to get damping const
            dynamicFriction = 1.0,
            staticFriction = 1.0,
            rollingFriction = 1.05,
            torsionalFriction = 1.3,
            shapeFactor = 0.5,

            nPerh = 1.01,                  # this should basically always be 1 for DEM
            
            # integration
            IntegratorConstructor = VerletIntegrator,
            stepsPerCollision = 50,  # replaces CFL for DEM
            goalTime = None,
            dt = 1e-8,
            dtMin = 1.0e-8, 
            dtMax = 0.1,
            dtGrowth = 2.0,
            steps = 100,
            maxSteps = None,
            statsStep = 10,
            domainIndependent = False,
            rigorousBoundaries = False,
            dtverbose = False,
            
            # output control
            vizCycle = 10,#None,
            vizTime = None,
            clearDirectories = False,
            restoreCycle = None,
            restartStep = 1000,
            redistributeStep = 500,
            dataDir = "dumps-DEM-3d",

            # test rotation on top of restitution coefficient
            boolCheckSlidingFriction=False,
            boolCheckRollingFriction=False,
            boolCheckTorsionalFriction=False,
            boolCheckTorsionalObjectivity=False,

             # ats parameters
            checkError = False,                # turn on error checking for restitution coefficient
            checkRestart = False,              # turn on error checking for restartability
            checkConservation = False,         # turn on error checking for momentum conservation
            restitutionErrorThreshold = 0.01,  # relative error actual restitution vs nominal
            conservationErrorThreshold = 1e-15 # relative error for momentum conservation
            )

# assert sum([checkError,
#             boolCheckSlidingFriction,
#             boolCheckRollingFriction,
#             boolCheckTorsionalFriction,
#             boolCheckTorsionalObjectivity]) <= 1

#-------------------------------------------------------------------------------
# file things
#-------------------------------------------------------------------------------
testName = "DEM-twoParticleCollision-3d"
restartDir = os.path.join(dataDir, "restarts")
vizDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, testName)
vizBaseName = testName

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
# Set the node properties.
#-------------------------------------------------------------------------------
generator1 = GenerateNodeDistribution3d(2, 1, 1,
                                        rho = 1.0,
                                        distributionType = "lattice",
                                        xmin = (0.0,  0.0, 0.0),
                                        xmax = (1.0,  0.5, 0.5),
                                        nNodePerh = nPerh)


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
hydro = DEM(db,
            normalSpringConstant = normalSpringConstant,
            normalRestitutionCoefficient = normalRestitutionCoefficient,
            tangentialSpringConstant = tangentialSpringConstant,
            tangentialRestitutionCoefficient = tangentialRestitutionCoefficient,
            dynamicFrictionCoefficient = dynamicFriction,
            staticFrictionCoefficient = staticFriction,
            rollingFrictionCoefficient = rollingFriction,
            torsionalFrictionCoefficient = torsionalFriction,
            shapeFactor = shapeFactor,
            stepsPerCollision = stepsPerCollision)

packages = [hydro]


#-------------------------------------------------------------------------------
# initial conditions
#-------------------------------------------------------------------------------

velocity = nodes1.velocity()
velocity[0] = Vector(vImpact,0.0,0.0)
velocity[1] = Vector(-vImpact,0.0,0.0)

particleRadius = nodes1.particleRadius()
particleRadius[0] = radius
particleRadius[1] = radius
omega = hydro.omega

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
output("integrator.havePhysicsPackage(hydro)")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.domainDecompositionIndependent")
output("integrator.rigorousBoundaries")
output("integrator.verbose")

#-------------------------------------------------------------------------------
# Periodic Work Function : track conservation
#-------------------------------------------------------------------------------
conservation = TrackConservation(db,
                                  hydro,
                                  verbose=False)
                                  
periodicWork = [(conservation.periodicWorkFunction,1)]


#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            iterateInitialH = False,
                            initializeDerivatives = True,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            vizBaseName = vizBaseName,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime,
                            periodicWork=periodicWork)
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
    restitutionError = abs(restitutionEff + restitutionCoefficient)/restitutionCoefficient
    if  restitutionError > restitutionErrorThreshold:
        raise ValueError, "relative restitution coefficient error, %g, exceeds bounds" % restitutionError

if checkConservation:
    if  conservation.deltaLinearMomentumX() > conservationErrorThreshold:
        raise ValueError, "linear momentum - x conservation error, %g, exceeds bounds" % conservation.deltaLinearMomentumX()
    if  conservation.deltaLinearMomentumY() > conservationErrorThreshold:
        raise ValueError, "linear momentum - y conservation error, %g, exceeds bounds" % conservation.deltaLinearMomentumY()
    if  conservation.deltaLinearMomentumZ() > conservationErrorThreshold:
        raise ValueError, "linear momentum - z conservation error, %g, exceeds bounds" % conservation.deltaLinearMomentumZ()
    if  conservation.deltaRotationalMomentumX() > conservationErrorThreshold:
        raise ValueError, "rotational momentum - x conservation error, %g, exceeds bounds" % conservation.deltaRotationalMomentumX()
    if  conservation.deltaRotationalMomentumY() > conservationErrorThreshold:
        raise ValueError, "rotational momentum - y conservation error, %g, exceeds bounds" % conservation.deltaRotationalMomentumY()
    if  conservation.deltaRotationalMomentumZ() > conservationErrorThreshold:
        raise ValueError, "rotational momentum -z conservation error, %g, exceeds bounds" % conservation.deltaRotationalMomentumZ()