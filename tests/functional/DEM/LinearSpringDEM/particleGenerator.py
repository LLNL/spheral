#ATS:DEM3dImpact = test(          SELF, "--clearDirectories True --checkConservation True   --goalTime 1.0", label="DEM impacting squares -- 3-D (parallel)", np=8)

import os, sys, shutil, mpi, random
from math import *
from Spheral3d import *
from SpheralTestUtilities import *
from findLastRestart import *
from GenerateNodeDistribution3d import *
from GenerateDEMfromSPHGenerator import GenerateDEMfromSPHGenerator3d

from DEMConservationTracker import TrackConservation3d as TrackConservation

if mpi.procs > 1:
    from PeanoHilbertDistributeNodes import distributeNodes3d
else:
    from DistributeNodes import distributeNodes3d

title("DEM 3d Impacting Squares")
# This tests the conservation properties of the DEM package when
# distribution across multiple processors

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(numParticlePerLength = 3,                 # number of particles on a side of the box
            radius = 0.25,                            # particle radius
            normalSpringConstant=1000.0,              # spring constant for LDS model
            normalRestitutionCoefficient=0.55,        # restitution coefficient to get damping const
            tangentialSpringConstant=285.70,          # spring constant for LDS model
            tangentialRestitutionCoefficient=0.55,    # restitution coefficient to get damping const
            dynamicFriction = 1.0,                    # static friction coefficient sliding
            staticFriction = 1.0,                     # dynamic friction coefficient sliding
            rollingFriction = 1.05,                   # static friction coefficient for rolling
            torsionalFriction = 1.3,                  # static friction coefficient for torsion
            cohesiveTensileStrength = 0.0,            # units of pressure
            shapeFactor = 0.1,                        # in [0,1] shape factor from Zhang 2018, 0 - no torsion or rolling
            
            # integration
            IntegratorConstructor = VerletIntegrator,
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
            vizTime = None, 
            clearDirectories = False,
            restoreCycle = None,
            restartStep = 10000,
            redistributeStep = 1000000000,
            dataDir = "dumps-DEM-impactingSquares-3d",

            # ats
            checkRestart = False,
            checkConservation = False,             # turn on error checking for momentum conservation
            conservationErrorThreshold = 2e-14,    # relative error for momentum conservation  
            )

#-------------------------------------------------------------------------------
# file things
#-------------------------------------------------------------------------------
testName = "DEM-ImpactingSquares-3d"
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
# Make the NodeList.
#-------------------------------------------------------------------------------
units = CGuS()
nodes1 = makeDEMNodeList("nodeList1",)
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
# if restoreCycle is None:
#     generator0 = GenerateNodeDistribution3d(numParticlePerLength, numParticlePerLength, numParticlePerLength,
#                                             rho = 1.0,
#                                             distributionType = "lattice",
#                                             xmin = (0, 0.0, 1.0),
#                                             xmax = ( 1.0,1.0, 2.0),
#                                             nNodePerh = nPerh)
    
#     # # really simple bar shaped particle
#     # def DEMParticleGenerator(xi,yi,zi,Hi,mi,Ri):
#     #     xout = [xi]
#     #     yout = [yi]
#     #     zout = [zi]
#     #     mout = [mi/1.1]
#     #     Rout = [Ri/1.1]
#     #     return xout,yout,zout,mout,Rout

#     generator1 = GenerateDEMfromSPHGenerator3d(WT,
#                                                generator0,
#                                                #DEMParticleGenerator=DEMParticleGenerator,
#                                                nPerh=nPerh)

#     distributeNodes3d((nodes1, generator1))
   
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
          normalSpringConstant = normalSpringConstant,
          normalRestitutionCoefficient = normalRestitutionCoefficient,
          tangentialSpringConstant = tangentialSpringConstant,
          tangentialRestitutionCoefficient = tangentialRestitutionCoefficient,
          dynamicFrictionCoefficient = dynamicFriction,
          staticFrictionCoefficient = staticFriction,
          rollingFrictionCoefficient = rollingFriction,
          torsionalFrictionCoefficient = torsionalFriction,
          cohesiveTensileStrength = cohesiveTensileStrength,
          shapeFactor = shapeFactor,
          stepsPerCollision = stepsPerCollision)

packages = [dem]


# #-------------------------------------------------------------------------------
# # PhysicsPackage : gravity
# #-------------------------------------------------------------------------------
gravity = ConstantAcceleration(a0 = Vector(0.0, 0.0, -1.0),
                               nodeList = nodes1)
packages += [gravity]

# #-------------------------------------------------------------------------------
# # Create boundary conditions.
# #-------------------------------------------------------------------------------
plane1 = Plane(Vector(0.0, 0.0, 0.0), Vector(  0.0, 0.0, 1.0))
bc1 = ReflectingBoundary(plane1)
bcSet = [bc1]

for p in packages:
    for bc in bcSet:
        p.appendBoundary(bc)

#-------------------------------------------------------------------------------
# Fields and Variables
#-------------------------------------------------------------------------------
numNodeLists = db.numNodeLists
nodeLists = db.nodeLists()

position = db.DEMPosition
mass = db.DEMMass
velocity = db.DEMVelocity
H = db.DEMHfield
radius = db.DEMParticleRadius
compositeParticleIndex = db.DEMCompositeParticleIndex

uniqueIndex = db.DEMUniqueIndex
omega = dem.omega

#-------------------------------------------------------------------------------
# Initial Conditions
#-------------------------------------------------------------------------------
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
# Periodic Work Function: Track conseravation
#-------------------------------------------------------------------------------

conservation = TrackConservation(db,
                                  dem,
                                  verbose=True)
                                  
periodicWork = [(conservation.periodicWorkFunction,1)]


#-------------------------------------------------------------------------------
# Periodic Work Function: DEM portal that controls the inflow
#-------------------------------------------------------------------------------
class DEMInflow:
    def __init__(self,
                 nodeList,            # node list we'll add new nodes to
                 db,                  # we'll need a couple of the database methods 
                 dem,                 # dem object has some fields we need
                 radius,              # radius of the inflow region
                 position,            # position (center) of inflow
                 normal,              # inflow direction
                 inflowVelocity,      # velocity of particles inflowing
                 maxParticleRadius):  # max radius of particles
            
        self.nodeList = nodeList
        self.db = db
        self.dem = dem
        self.radius = radius
        self.position = position
        self.normal = normal
        self.inflowVelocity = inflowVelocity
        self.maxParticleRadius = maxParticleRadius

    	self.dtGen = maxParticleRadius/inflowVelocity.magnitude()*1.5
    	self.tGen = -self.dtGen
  	
    def addParticles(self,particles):
        
        numParticles = len(particles)
    	cId0 = max(self.db.DEMCompositeParticleIndex.max(),-1)
    	uId0 = max(self.db.DEMUniqueIndex.max(),-1)

        if mpi.rank == 0:

            # set the particle unique and composite indices
            numNewParticles = 0
    	    for particle in particles:
    	        numNewParticles += particle.numSubParticles

            # fields
    	    mas = self.nodeList.mass()
            rad = self.nodeList.particleRadius()
            pos = self.nodeList.positions()
            vel = self.nodeList.velocity()
            cId = self.nodeList.compositeParticleIndex()
            uId = self.nodeList.uniqueIndex()
            H = self.nodeList.Hfield()
            # initial new nodes and fill in the data
            k = self.nodeList.numInternalNodes
            cIdi = cId0+1
            uIdi = uId0+1

            self.nodeList.numInternalNodes += numNewParticles

            for i in range(numParticles):
                for j in range(particles[i].numSubParticles):
                    print k, uId[k], uIdi
                    mas[k] = particles[i].mass[j]
                    rad[k] = particles[i].radius[j]
                    pos[k] = particles[i].position[j]
                    vel[k] = self.inflowVelocity+Vector((random.random()-0.5),(random.random()-0.5),(random.random()-0.5))
                    #H[k] = 1.0/(2.0 * rad[k] *1.1) * SymTensor.one
                    cId[k] = 1*cIdi
                    uId[k] = 1*uIdi
                    k += 1
                    uIdi += 1
                cIdi += 1   

        self.db.setDEMHfieldFromParticleRadius(uId0+1)
        self.db.reinitializeNeighbors()
    	self.db.updateConnectivityMap()
    	self.dem.updateContactMap(db)
        self.dem.resizePairFieldLists()
    	self.dem.initializeOverlap(db,cId0+1)

    def __call__(self,cycle,time,dt):

        print self.nodeList.numInternalNodes
        if (time - self.tGen > self.dtGen):

            particleMass  = 1.0
            particleScale = self.maxParticleRadius*2
            spawnPoint    = self.position

            particle = CubeParticle(particleMass,
                                    spawnPoint,
                                    particleScale)

            self.addParticles([particle])

            self.tGen += self.dtGen

class CubeParticle:
    def __init__(self,
                 cubeMass,
                 cubeCenterOfMass,
                 cubeSideLength,
                 overLapPercentage = 40,
                 numParticlesPerLength = 4):

        self.numSubParticles = numParticlesPerLength ** 3
        self.mass=[]
        self.radius=[]
        self.position=[]

        f = 1+overLapPercentage/200.0
        
        massi = cubeMass/float(self.numSubParticles)
        radiusi = f*cubeSideLength/float(2*numParticlePerLength+overLapPercentage/100)

        s0 = -cubeSideLength/2.0 + radiusi
        s1 = -s0
        ds  = (s1-s0)/float(numParticlesPerLength-1.0)

        for i in range(self.numSubParticles):

            ix = int(i%numParticlesPerLength)
            iy = int(i/numParticlesPerLength%numParticlesPerLength)
            iz = int(i/numParticlesPerLength**2%numParticlesPerLength)

            xi = s0 + ix*ds + cubeCenterOfMass[0]
            yi = s0 + iy*ds + cubeCenterOfMass[1]
            zi = s0 + iz*ds + cubeCenterOfMass[2]

            self.mass.append(massi)
            self.radius.append(radiusi)
            self.position.append(Vector(xi,yi,zi))
            

inflow = DEMInflow(nodes, db, dem,
                 1.0,
                 Vector(0.0,0.0, 10.1),
                 Vector(0.0,0.0,-1.0),
                 Vector(0.0,0.0, -1.0),
                 1.0)
periodicWork += [(inflow,1)]   
# MPI = 1 


    #  	
# in the finalize step do the DEM solid boundary
    
# hash map for octree with spheres 
# want to be able to set the mass flow rate
# octree after first particle insertion
# random point in a disk




#-------------------------------------------------------------------------------
# quadtree to initialize particle and set a fixed inflow mass
#
#-------------------------------------------------------------------------------
# use an quadtree to efficiently pack
# make a certain number of available slots and randomly assign to those slots
#-------------------------------------------------------------------------------
# maybe something like particleGenerator.generateParticle() when some criterion
# occurs. Have a position and radius for an inflow condition with a flow rate 
# max flow rate limited by the size of the radius and the particle size.    
#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
from SpheralPointmeshSiloDump import dumpPhysicsState
control = SpheralController(integrator,
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
                            periodicWork=periodicWork)
output("control")

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------

if not steps is None:
    control.step(steps)
else:
    control.advance(goalTime, maxSteps)

