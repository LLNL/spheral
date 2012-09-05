#-------------------------------------------------------------------------------
# The Taylor anvil impact problem -- impact of a solid cylinder on an unyielding
# surface.
#-------------------------------------------------------------------------------
from Numeric import *
from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from SpheralController import *
from findLastRestart import *
from SpheralVisitDump import SpheralVisitDump
from math import *

from Strength import *

# Load the mpi module if we're parallel.
mpi, procID, numProcs = loadmpi()

#-------------------------------------------------------------------------------
# Generic problem parameters
# All CGS units.
#-------------------------------------------------------------------------------
seed = "cylindrical"

zlength, rlength = 4.694, 0.381
nz, nr = 129, 11 # 258, 21
nPerh = 2.01

rmin = 0.0
rmax = rlength
thetamin = 0.0
thetamax = 0.5*pi
zmin = 0.0
zmax = zlength

rho0 = 16.69
v0 = Vector3d(0.0, 0.0, -2.5e4)

Qconstructor = MonaghanGingoldViscosity3d
Cl, Cq = 1.0, 1.0
Qlimiter = False
balsaraCorrection = False
epsilon2 = 1e-4
negligibleSoundSpeed = 1e-5
csMultiplier = 1e-4
HsmoothMin, HsmoothMax = 1e-5, 0.5
cfl = 0.5
useVelocityMagnitudeForDt = False
XSPH = False
epsilonTensile = 0.0
nTensile = 4
hybridMassDensityThreshold = 0.01

neighborSearchType = Neighbor3d.NeighborSearchType.GatherScatter
numGridLevels = 20
topGridCellSize = 1.0
origin = Vector3d(0.0, 0.0, 0.0)

goalTime = 150.0e-6
dtSample = 1e-6
dt = 1e-10
dtMin, dtMax = 1e-12, 1e-5
dtGrowth = 2.0
maxSteps = None
statsStep = 10
smoothIters = 0
HEvolution = Hydro3d.HEvolutionType.IdealH
sumForMassDensity = Hydro3d.MassDensityType.IntegrateDensity # HybridDensity # CorrectedSumDensity

restartStep = 1000
restartBaseName = "dumps-3d/TaylorImpact-%i-%i" % (nr, nz)
restoreCycle = findLastRestart(restartBaseName)

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("3-D Taylor anvil impact strength test")

#-------------------------------------------------------------------------------
# Tantalum material properties.
#-------------------------------------------------------------------------------
eosTantalum = GruneisenEquationOfStateCGS3d(16.69,   # reference density  
                                            -0.2,     # etamin             
                                            4.0,     # etamax             
                                            0.341e6, # C0                 
                                            1.2,     # S1                 
                                            0.0,     # S2                 
                                            0.0,     # S3                 
                                            1.67,    # gamma0             
                                            0.42,    # b                  
                                            180.948) # atomic weight
coldFitTantalum = NinthOrderPolynomialFit(-6.86446714e9,
                                          -1.17070812e10,
                                           9.70252276e11,
                                          -4.12557402e11,
                                           1.13401823e11,
                                          -1.86584799e10,
                                           0.0,
                                           0.0,
                                           0.0,
                                           0.0)
meltFitTantalum = NinthOrderPolynomialFit(9.24414908e10,
                                          2.53949977e11,
                                          1.06113848e12,
                                         -5.28947636e11,
                                          1.67906438e11,
                                         -2.92459765e10,
                                          0.0,
                                          0.0,
                                          0.0,
                                          0.0)
strengthTantalum = SteinbergGuinanStrengthCGS3d(eosTantalum,
                                                6.900000e11,        # G0
                                                1.4500e-12,         # A
                                                1.3000e-04,         # B
                                                7.7000e9,           # Y0
                                                1.1e10,             # Ymax
                                                1.0e-3,             # Yp
                                                10.0000,            # beta
                                                0.0,                # gamma0
                                                0.1,                # nhard
                                                coldFitTantalum,
                                                meltFitTantalum)

#-------------------------------------------------------------------------------
# Create the NodeLists.
#-------------------------------------------------------------------------------
nodes = SphSolidNodeList3d("Tantalum", eosTantalum, strengthTantalum)
nodes.nodesPerSmoothingScale = nPerh
nodes.epsilonTensile = epsilonTensile
nodes.nTensile = nTensile
nodes.XSPH = XSPH
output('nodes.name()')
output('nodes.nodesPerSmoothingScale')
output('nodes.epsilonTensile')
output('nodes.nTensile')
output('nodes.XSPH')

#-------------------------------------------------------------------------------
# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
#-------------------------------------------------------------------------------
WT = TableKernel3d(BSplineKernel3d(), 1000)
WTPi = TableKernel3d(BSplineKernel3d(), 1000)
output('WT')
output('WTPi')
kernelExtent = WT.kernelExtent()

#-------------------------------------------------------------------------------
# Construct the neighbor objects and associate them with the node lists.
#-------------------------------------------------------------------------------
neighborTimer = SpheralTimer('Neighbor initialization.')
neighborTimer.start()
neighbor = NestedGridNeighbor3d(nodes,
                                neighborSearchType,
                                numGridLevels,
                                topGridCellSize,
                                origin,
                                kernelExtent)
nodes.registerNeighbor(neighbor)
neighborTimer.stop()
neighborTimer.printStatus()

#-------------------------------------------------------------------------------
# Set node properties (positions, masses, H's, etc.)
#-------------------------------------------------------------------------------
if restoreCycle is None:
    from GenerateNodeDistribution3d import *
    from DistributeNodes import distributeNodes3d
    print "Generating node distribution."
    generator = GenerateNodeDistribution3d(nr,
                                           nz,
                                           0,
                                           rho = rho0,
                                           distributionType = seed,
                                           rmin = rmin,
                                           rmax = rmax,
                                           thetamin = thetamin,
                                           thetamax = thetamax,
                                           zmin = zmin,
                                           zmax = zmax,
                                           nNodePerh = nPerh)
    nTantalum = generator.globalNumNodes()
    nodeInfo = distributeNodes3d([(nodes, nTantalum, generator)])
    if mpi:
        output('mpi.reduce(nodes.numInternalNodes, mpi.MIN)')
        output('mpi.reduce(nodes.numInternalNodes, mpi.MAX)')
        output('mpi.reduce(nodes.numInternalNodes, mpi.SUM)')
    else:
        output('nodes.numInternalNodes')
    assert len(nodeInfo[nodes]['globalNodeListID']) == nodes.numInternalNodes

    # Set the node mass densities.
    nodes.setMassDensity(ScalarField3d("tmp", nodes, rho0))

    # Set node specific thermal energies
    nodes.setSpecificThermalEnergy(ScalarField3d("tmp", nodes,
                                                 eosTantalum.specificThermalEnergy(rho0, 300.0)))

    # Set the node velocities.
    nodes.setVelocity(VectorField3d("tmp", nodes, v0))

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
plane0 = Plane3d(Vector3d(0.0, 0.0, 0.0), Vector3d(1.0, 0.0, 0.0))
plane1 = Plane3d(Vector3d(0.0, 0.0, 0.0), Vector3d(0.0, 1.0, 0.0))
plane2 = Plane3d(Vector3d(0.0, 0.0, 0.0), Vector3d(0.0, 0.0, 1.0))
bc0 = ReflectingBoundary3d(plane0)
bc1 = ReflectingBoundary3d(plane1)
bc2 = ReflectingBoundary3d(plane2)

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase3d()
db.appendNodeList(nodes)
output('db')
output('db.numNodeLists')
output('db.numFluidNodeLists')

#-------------------------------------------------------------------------------
# Construct the artificial viscosities for the problem.
#-------------------------------------------------------------------------------
q = Qconstructor(Cl, Cq)
q.limiter = Qlimiter
q.balsaraShearCorrection = balsaraCorrection
q.epsilon2 = epsilon2
q.negligibleSoundSpeed = negligibleSoundSpeed
q.csMultiplier = csMultiplier
output('q')
output('q.Cl')
output('q.Cq')
output('q.limiter')
output('q.epsilon2')
output('q.negligibleSoundSpeed')
output('q.csMultiplier')
output('q.balsaraShearCorrection')

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
hydro = Hydro3d(WT, WTPi, q)
hydro.cfl = cfl
hydro.useVelocityMagnitudeForDt = True
hydro.HEvolution = HEvolution
hydro.sumForMassDensity = sumForMassDensity
hydro.HsmoothMin = HsmoothMin
hydro.HsmoothMax = HsmoothMax
hydro.hybridMassDensityThreshold = hybridMassDensityThreshold
output('hydro')
output('hydro.cfl')
output('hydro.useVelocityMagnitudeForDt')
output('hydro.HEvolution')
output('hydro.sumForMassDensity')
output('hydro.HsmoothMin')
output('hydro.HsmoothMax')
output('hydro.kernel()')
output('hydro.PiKernel()')
output('hydro.valid()')
output('hydro.hybridMassDensityThreshold')

#-------------------------------------------------------------------------------
# Construct a strength physics object.
#-------------------------------------------------------------------------------
strength = Strength3d()
output("strength")

#-------------------------------------------------------------------------------
# Construct a predictor corrector integrator, and add the physics packages.
#-------------------------------------------------------------------------------
integrator = PredictorCorrectorIntegrator3d(db)
integrator.appendPhysicsPackage(hydro)
integrator.appendPhysicsPackage(strength)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
output('integrator')
output('integrator.havePhysicsPackage(hydro)')
output('integrator.havePhysicsPackage(strength)')
output('integrator.valid()')
output('integrator.lastDt')
output('integrator.dtMin')
output('integrator.dtMax')
output('integrator.dtGrowth')

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            boundaryConditions = [bc0, bc1, bc2],
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            initializeMassDensity = False)
output('control')

#-------------------------------------------------------------------------------
# Smooth the initial conditions/restore state.
#-------------------------------------------------------------------------------
if restoreCycle is not None:
    control.loadRestartFile(restoreCycle)
else:
    control.smoothState(smoothIters)

    # Viz the initial conditions.
    P = db.fluidPressure
    cs = db.fluidSoundSpeed
    Hi = db.fluidHinverse
    vx = ScalarField3d("x velocity", nodes)
    vy = ScalarField3d("y velocity", nodes)
    vz = ScalarField3d("z velocity", nodes)
    for i in xrange(nodes.numInternalNodes):
        vx[i] = nodes.velocity()[i].x
        vy[i] = nodes.velocity()[i].y
        vz[i] = nodes.velocity()[i].z
    dumper = SpheralVisitDump(db,
                              "TaylorImpact-3d-visit",
                              "dumps-3d",
                              listOfFieldLists = [db.fluidMassDensity,
                                                  db.fluidVelocity,
                                                  db.fluidWeight,
                                                  db.fluidSpecificThermalEnergy,
                                                  P,
                                                  cs,
                                                  Hi],
                              listOfFields = [nodes.deviatoricStress(),
                                              nodes.plasticStrain(),
                                              vx,
                                              vy,
                                              vz]
                              )
    dumper.dump(control.time(), control.totalSteps)

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
while control.time() < goalTime:
    nextGoalTime = min(control.time() + dtSample, goalTime)
    control.advance(nextGoalTime, maxSteps)
    control.dropRestartFile()

    # Viz the current state.
    P = db.fluidPressure
    cs = db.fluidSoundSpeed
    Hi = db.fluidHinverse
    vx = ScalarField3d("x velocity", nodes)
    vy = ScalarField3d("y velocity", nodes)
    vz = ScalarField3d("z velocity", nodes)
    for i in xrange(nodes.numInternalNodes):
        vx[i] = nodes.velocity()[i].x
        vy[i] = nodes.velocity()[i].y
        vz[i] = nodes.velocity()[i].z
    dumper = SpheralVisitDump(db,
                              "TaylorImpact-3d-visit",
                              "dumps-3d",
                              listOfFieldLists = [db.fluidMassDensity,
                                                  db.fluidVelocity,
                                                  db.fluidWeight,
                                                  db.fluidSpecificThermalEnergy,
                                                  P,
                                                  cs,
                                                  Hi],
                              listOfFields = [nodes.deviatoricStress(),
                                              nodes.plasticStrain(),
                                              vx,
                                              vy,
                                              vz]
                              )
    dumper.dump(control.time(), control.totalSteps)
