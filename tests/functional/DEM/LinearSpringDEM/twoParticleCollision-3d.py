#ATS:DEM3d = test(        SELF, "--clearDirectories True  --checkError True  --restitutionCoefficient=1.0 --steps 100", label="DEM idividual particle collision -- 1-D (serial)")

import os, sys, shutil, mpi
from math import *
from Spheral3d import *
from SpheralTestUtilities import *
from findLastRestart import *
from GenerateNodeDistribution3d import *

if mpi.procs > 1:
    from PeanoHilbertDistributeNodes import distributeNodes3d
else:
    from DistributeNodes import distributeNodes3d

title("DEM Restitution Coefficient Test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(vImpact = 1.0,                 # impact velocity
            normalSpringConstant=10000.0,  # spring constant for LDS model
            restitutionCoefficient=0.8,    # restitution coefficient to get damping const
            radius = 0.25,                 # particle radius
            nPerh = 1.01,                  # this should basically always be 1 for DEM

            # integration
            IntegratorConstructor = CheapSynchronousRK2Integrator,
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

            # ats parameters
            checkError = False,
            checkRestart = False,
            restitutionErrorThreshold = 0.01, # relative error
            )

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
if restoreCycle is None:
    generator1 = GenerateNodeDistribution3d(2, 1, 1,
                                            rho = 1.0,
                                            distributionType = "lattice",
                                            xmin = (0.0,  0.0, 0.0),
                                            xmax = (1.0,  0.5, 0.5),
                                            nNodePerh = nPerh)


    distributeNodes3d((nodes1, generator1))

    # initial conditions
    velocity = nodes1.velocity()
    velocity[0] = Vector(vImpact,0.0,0.0)
    velocity[1] = Vector(-vImpact,0.0,0.0)

    particleRadius = nodes1.particleRadius()

    particleRadius[0] = radius
    particleRadius[1] = radius

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
            normalSpringConstant,
            restitutionCoefficient,
            stepsPerCollision = stepsPerCollision)

packages = [hydro]

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
                            SPH = SPH)
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