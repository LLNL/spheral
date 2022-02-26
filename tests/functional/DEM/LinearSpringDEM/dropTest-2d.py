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

title("DEM 2d Drop Test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(numParticlePerLength = 50,     # number of particles on a side of the box
            normalSpringConstant=1000.0,  # spring constant for LDS model
            restitutionCoefficient=0.8,    # restitution coefficient to get damping const
            radius = 0.25,                 # particle radius
            nPerh = 1.01,                  # this should basically always be 1 for DEM

            # integration
            IntegratorConstructor = CheapSynchronousRK2Integrator,
            stepsPerCollision = 50,  # replaces CFL for DEM
            goalTime = 25.0,
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
            dataDir = "dumps-DEM-2d",
            )

#-------------------------------------------------------------------------------
# file things
#-------------------------------------------------------------------------------
testName = "DEM-dropTest-2d"
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
    generator1 = GenerateNodeDistribution2d(numParticlePerLength, numParticlePerLength,
                                            rho = 1.0,
                                            distributionType = "lattice",
                                            xmin = (-0.5,  0.525),
                                            xmax = ( 0.5,  1.525),
                                            nNodePerh = nPerh)


    distributeNodes2d((nodes1, generator1))

    # set our particle radius
    radius =1.0/numParticlePerLength/2.5
    particleRadius = nodes1.particleRadius()
    nodes1.particleRadius(ScalarField("radii", nodes1, radius))

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
# PhysicsPackage : DEM
#-------------------------------------------------------------------------------
dem = DEM(db,
          normalSpringConstant,
          restitutionCoefficient,
          stepsPerCollision = stepsPerCollision)

packages = [dem]

#-------------------------------------------------------------------------------
# PhysicsPackage : gravity
#-------------------------------------------------------------------------------
allindices = vector_of_int(range(nodes1.numInternalNodes))

gravity = ConstantAcceleration(a0 = Vector(0.0,-1.0),
                               nodeList = nodes1,
                               indices = allindices)
packages += [gravity]

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
plane1 = Plane(Vector(0.0, 0.0), Vector(  1.0, 1.0))
plane2 = Plane(Vector(0.0, 0.0), Vector( -1.0, 1.0))
bc1 = ReflectingBoundary(plane1)
bc2 = ReflectingBoundary(plane2)
bcSet = [bc1, bc2]

for p in packages:
    for bc in bcSet:
        p.appendBoundary(bc)

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
                            redistributeStep=redistributeStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            vizBaseName = vizBaseName,
                            vizMethod = dumpPhysicsState,
                            vizGhosts=True,
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

