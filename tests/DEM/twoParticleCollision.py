#ATS:DEM0 = test(        SELF, "--clearDirectories True  --checkError True  --restartStpe 10 --step 50", label="DEM idividual particle collision -- 1-D (serial)")
#ATS:DEM1 = testif(DEM0, SELF, "--clearDirectories False --checkError False  --restartStep 10 --restoreCycle 10 --steps 10 --checkRestart True", label="DEM idividual particle collision -- 1-D (serial) RESTART CHECK")

import shutil
from math import *
from Spheral2d import *
from SpheralTestUtilities import *
#from SpheralGnuPlotUtilities import *
from findLastRestart import *
from GenerateNodeDistribution2d import *

import mpi
import DistributeNodes

title("Testing Out DEM Machinery")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(vImpact = 1.0,
            youngsModulus=10000.0,
            restitutionCoefficient=0.8,
            radius = 1.00,
            nPerh = 1.00,

            IntegratorConstructor = CheapSynchronousRK2Integrator,
            goalTime = None,
            steps = 50,
            vizCycle = 1,
            vizTime = None,
            dt = 1e-8,
            dtMin = 1.0e-8, 
            dtMax = 0.1,
            dtGrowth = 2.0,
            maxSteps = None,
            statsStep = 10,
            smoothIters = 0,
            domainIndependent = False,
            rigorousBoundaries = False,
            dtverbose = False,
            cfl = 0.05,

            useVoronoiOutput = False,
            clearDirectories = False,
            restoreCycle = None,
            restartStep = 1000,
            redistributeStep = 500,
            dataDir = "dumps-DEM-2d",
            
            serialDump = False, #whether to dump a serial ascii file at the end for viz

            checkError = False,
            checkRestart = False,
            restitutionErrorThreshold = 0.01, # relative error
            )

testName = "DEM-twoParticleCollision-2d"
restartDir = os.path.join(dataDir, "restarts")
vizDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, testName)
vizBaseName = testName

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
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
# Make the NodeList.
#-------------------------------------------------------------------------------
units = CGuS()
nodes1 = makeDEMNodeList("nodeList1",
                          hmin = 1.0e-30,
                          hmax = 1.0e30,
                          hminratio = 100.0,
                          nPerh = nPerh)
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
    generator1 = GenerateNodeDistribution2d(2, 1,
                                            rho = 1.0,
                                            distributionType = "lattice",
                                            xmin = (0.0,  0.0),
                                            xmax = (1.0,  0.5),
                                            nNodePerh = nPerh)

    if mpi.procs > 1:
        from VoronoiDistributeNodes import distributeNodes2d
    else:
        from DistributeNodes import distributeNodes2d

    distributeNodes2d((nodes1, generator1))

    # initial conditions
    velocity = nodes1.velocity()
    velocity[0] = Vector(vImpact,0.0)
    velocity[1] = Vector(-vImpact,0.0)
    
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
# This doesn't really matter kernel filler for neighbor algo
#-------------------------------------------------------------------------------
WT = TableKernel(WendlandC2Kernel(), 1000)

#-------------------------------------------------------------------------------
# DEM
#-------------------------------------------------------------------------------
DLS = DampedLinearSpring(youngsModulus,
                         restitutionCoefficient)
hydro = DEM(db, 
            WT, 
            cfl = cfl)

hydro.appendContactModel(DLS)

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