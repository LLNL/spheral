#ATS:DEM2d0COMP = test(              SELF, "--clearDirectories True  --checkError True  --restartStep 10 --steps 100", label="DEM composite particle collision -- 2-D (serial)")
#ATS:DEM2d1COMP = testif(DEM2d0COMP, SELF, "--clearDirectories False --checkError False  --restartStep 10 --restoreCycle 10 --steps 90 --checkRestart True", label="DEM composite particle collision -- 2-D (serial) RESTART CHECK")
#ATS:DEM2d2COMP = test(              SELF, "--clearDirectories True  --checkError True  --restitutionCoefficient=1.0 --steps 100", label="DEM composite particle collision -- 2-D (serial)")

import os, sys, shutil, mpi
from math import *
from Spheral2d import *
from SpheralTestUtilities import *
from findLastRestart import *
from GenerateNodeDistribution2d import *

if mpi.procs > 1:
    from PeanoHilbertDistributeNodes import distributeNodes2d
else:
    from DistributeNodes import distributeNodes2d

title("DEM Composite-Particle Restitution Coefficient Test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(vImpact = 1.0,                 # impact velocity
            normalSpringConstant=10000.0,  # spring constant for LDS model
            restitutionCoefficient=0.8,    # restitution coefficient to get damping const
            radius = 0.25,                 # particle radius
            
            neighborSearchBuffer = 0.1,    # multiplicative buffer to radius for neighbor search algo

            # integration
            IntegratorConstructor = Verlet,
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
            dataDir = "dumps-DEM-2d-CompositeParticleCollision",

            # ats parameters
            checkError = False,
            checkRestart = False,
            restitutionErrorThreshold = 0.01, # relative error
            )

#-------------------------------------------------------------------------------
# file things
#-------------------------------------------------------------------------------
testName = "DEM-CompositeParticleCollision-2d"
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
    generator1 = GenerateNodeDistribution2d(3, 1,
                                            rho = 1.0,
                                            distributionType = "lattice",
                                            xmin = (-0.5,  0.0),
                                            xmax = ( 1.0,  0.5))


    distributeNodes2d((nodes1, generator1))

    # initial conditions
    velocity = nodes1.velocity()
    particleRadius = nodes1.particleRadius()
    position = nodes1.positions()

    velocity[0] = Vector(vImpact,0.0)
    velocity[1] = Vector(vImpact,0.0)
    velocity[2] = Vector(-vImpact,0.0)

    position[0] = Vector(0.0,position[0][1])

    particleRadius[0] = radius
    particleRadius[1] = radius
    particleRadius[2] = radius

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

eqOverlap=hydro.equilibriumOverlap
neighborIndices = hydro.neighborIndices
uniqueIndices = hydro.uniqueIndices
shearDisplacement = hydro.shearDisplacement
DDtShearDisplacement = hydro.DDtShearDisplacement
neighborIndices[0][0]=vector_of_int([uniqueIndices(0,1)])
eqOverlap[0][0]=vector_of_double([0.25])

print(neighborIndices(0,0))
print eqOverlap(0,0)
print eqOverlap(0,0)
print uniqueIndices(0,0)
print uniqueIndices(0,1)
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

def printFunc(cycle,time,dt):
    print [eqOverlap(0,0),eqOverlap(0,1),eqOverlap(0,2)]
    print [uniqueIndices(0,0),uniqueIndices(0,1),uniqueIndices(0,2)]
    print [neighborIndices(0,0),neighborIndices(0,1),neighborIndices(0,2)]
    print [shearDisplacement(0,0),shearDisplacement(0,1),shearDisplacement(0,2)]
    print [DDtShearDisplacement(0,0),DDtShearDisplacement(0,1),DDtShearDisplacement(0,2)]
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
                            SPH = SPH,
                            periodicWork=[(printFunc,1)])
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