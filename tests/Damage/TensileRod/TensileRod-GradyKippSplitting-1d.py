#-------------------------------------------------------------------------------
# A rod of stainless steel undergoing tensile strain.  This is intended as a
# test of the cracking/failure models.
#
# This version imposes a constant velocity gradient across the rod, and uses
# a Grady-Kipp flaw distribution.
#
# See Benz & Asphaug (1994), Icarus, 107, 98
#-------------------------------------------------------------------------------
from Numeric import *
from SolidSpheral import *
from SpheralTestUtilities import *
from SpheralController import *
from findLastRestart import *
from SpheralVisitDump import dumpPhysicsState
from math import *

from SpheralMatplotlibUtilities import *
import pylab

# Load the mpi module if we're parallel.
import loadmpi
mpi, procID, numProcs = loadmpi.loadmpi()

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("1-D Tensile rod strength/damage model test")

#-------------------------------------------------------------------------------
# Generic problem parameters
# All CGS units.
#-------------------------------------------------------------------------------
seed = 'lattice'

xlength = 20.0
nx = 96
nPerh = 2.01

xmin = -0.5*xlength
xmax =  0.5*xlength

rho0 = 7.9
dx = xlength/nx

# Parameters for the time dependent strain and cracking.
v0 = 1.0e3
kWeibull = 57.7
mWeibull = 2.63
volume = xlength
randomSeed = 109482993

Qconstructor = MonaghanGingoldViscosity1d
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

neighborSearchType = Neighbor1d.NeighborSearchType.GatherScatter
numGridLevels = 20
topGridCellSize = 1.0
origin = Vector1d(-xlength)

goalTime = 500.0e-6
dtSample = 0.005*goalTime
dt = 1e-10
dtMin, dtMax = 1e-12, 1e-5
dtGrowth = 2.0
maxSteps = None
statsStep = 10
smoothIters = 0
HEvolution = Hydro1d.HEvolutionType.IdealH
sumForMassDensity = Hydro1d.MassDensityType.IntegrateDensity # HybridDensity # CorrectedSumDensity

restartStep = 1000
dataDir = "dumps-GradyKipp-1d-%i" % nx
restartDir = dataDir + "/restarts"
visitDir = dataDir + "/visit"
restartBaseName = restartDir + "/TensileRod-%i" % nx

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
# If we're restarting, find the set of most recent restart files.
#-------------------------------------------------------------------------------
restoreCycle = findLastRestart(restartBaseName)

#-------------------------------------------------------------------------------
# Stainless steel material properties.
#-------------------------------------------------------------------------------
eos = GruneisenEquationOfStateCGS1d(7.90,    # reference density  
                                    1e-4,    # etamin             
                                    6.0,     # etamax             
                                    0.457e6, # C0                 
                                    1.49,    # S1                 
                                    0.0,     # S2                 
                                    0.0,     # S3                 
                                    1.93,    # gamma0             
                                    0.5,     # b                  
                                    55.350)  # atomic weight
coldFit = NinthOrderPolynomialFit(-1.06797724e10,
                                  -2.06872020e10,
                                   8.24893246e11,
                                  -2.39505843e10,
                                  -2.44522017e10,
                                   5.38030101e10,
                                   0.0,
                                   0.0,
                                   0.0,
                                   0.0)
meltFit = NinthOrderPolynomialFit(7.40464217e10,
                                  2.49802214e11,
                                  1.00445029e12,
                                 -1.36451475e11,
                                  7.72897829e9,
                                  5.06390305e10,
                                  0.0,
                                  0.0,
                                  0.0,
                                  0.0)
strengthModel = SteinbergGuinanStrengthCGS1d(eos,
                                             7.700000e11,        # G0
                                             2.2600e-12,         # A
                                             4.5500-04,          # B
                                             3.4000e9,           # Y0
                                             2.5e10,             # Ymax
                                             1.0e-3,             # Yp
                                             43.0000,            # beta
                                             0.0,                # gamma0
                                             0.35,               # nhard
                                             coldFit,
                                             meltFit)

#-------------------------------------------------------------------------------
# Create the NodeLists.
#-------------------------------------------------------------------------------
nodes = SphSolidNodeList1d("Stainless steel", eos, strengthModel)
nodeSet = [nodes]
for n in nodeSet:
    n.nodesPerSmoothingScale = nPerh
    n.epsilonTensile = epsilonTensile
    n.nTensile = nTensile
    n.XSPH = XSPH
    output('n.name()')
    output('n.nodesPerSmoothingScale')
    output('n.epsilonTensile')
    output('n.nTensile')
    output('n.XSPH')
del n

#-------------------------------------------------------------------------------
# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
#-------------------------------------------------------------------------------
WT = TableKernel1d(BSplineKernel1d(), 1000)
WTPi = TableKernel1d(BSplineKernel1d(), 1000)
output('WT')
output('WTPi')
kernelExtent = WT.kernelExtent()

#-------------------------------------------------------------------------------
# Construct the neighbor objects and associate them with the node lists.
#-------------------------------------------------------------------------------
neighborTimer = SpheralTimer('Neighbor initialization.')
neighborTimer.start()
cache = []
for n in nodeSet:
    neighbor = NestedGridNeighbor1d(n,
                                    neighborSearchType,
                                    numGridLevels,
                                    topGridCellSize,
                                    origin,
                                    kernelExtent)
    n.registerNeighbor(neighbor)
    cache.append(neighbor)
neighborTimer.stop()
neighborTimer.printStatus()
del n

#-------------------------------------------------------------------------------
# Set node properties (positions, masses, H's, etc.)
#-------------------------------------------------------------------------------
if restoreCycle is None:
    print "Generating node distribution."
    from DistributeNodes import distributeNodesInRange1d
    distributeNodesInRange1d([(nodes, nx, rho0, (xmin, xmax))])
    output('mpi.reduce(nodes.numInternalNodes, mpi.MIN)')
    output('mpi.reduce(nodes.numInternalNodes, mpi.MAX)')
    output('mpi.reduce(nodes.numInternalNodes, mpi.SUM)')

    # Set node specific thermal energies
    u0 = eos.specificThermalEnergy(rho0, 300.0)
    nodes.setSpecificThermalEnergy(ScalarField1d("tmp", nodes, u0))

    # Set node velocites.
    for i in xrange(nodes.numInternalNodes):
        nodes.velocity()[i].x = nodes.positions()[i].x/(0.5*xlength)*v0

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase1d()
for n in nodeSet:
    db.appendNodeList(n)
del n
output('db')
output('db.numNodeLists')
output('db.numFluidNodeLists')

#-------------------------------------------------------------------------------
# Construct constant velocity boundary conditions to be applied to the rod ends.
#-------------------------------------------------------------------------------
xNodes = vector_of_int()
for i in xrange(nodes.numInternalNodes):
    if (nodes.positions()[i].x < -0.5*xlength + 4*dx or
        nodes.positions()[i].x >  0.5*xlength - 4*dx):
        xNodes.append(i)
print ("Selected %i constant velocity nodes." %
       (mpi.allreduce(len(xNodes), mpi.SUM)))
xbc = ConstantVelocityBoundary1d(nodes, xNodes)

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
hydro = Hydro1d(WT, WTPi, q)
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
strength = Strength1d()
output("strength")

#-------------------------------------------------------------------------------
# Construct a damage model.
#-------------------------------------------------------------------------------
damageModel = GradyKippParticleSplittingDamage1d(nodes,
                                                 kWeibull,
                                                 mWeibull,
                                                 volume,
                                                 randomSeed)
damageModel.excludeNodes = xNodes
output("damageModel")

#-------------------------------------------------------------------------------
# Construct a predictor corrector integrator.
#-------------------------------------------------------------------------------
integrator = PredictorCorrectorIntegrator1d(db)
integrator.appendPhysicsPackage(hydro)
integrator.appendPhysicsPackage(strength)
integrator.appendPhysicsPackage(damageModel)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
output('integrator')
output('integrator.havePhysicsPackage(hydro)')
output('integrator.havePhysicsPackage(strength)')
output('integrator.havePhysicsPackage(damageModel)')
output('integrator.valid()')
output('integrator.lastDt')
output('integrator.dtMin')
output('integrator.dtMax')
output('integrator.dtGrowth')

#-------------------------------------------------------------------------------
# Add the boundary conditions.
#-------------------------------------------------------------------------------
for package in integrator.physicsPackages():
    for bc in [xbc]:
        package.appendBoundary(bc)

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
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

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
while control.time() < goalTime:
    nextGoalTime = min(control.time() + dtSample, goalTime)
    control.advance(nextGoalTime, maxSteps)
    control.dropRestartFile()

# Plot the state.
rhoPlot = plotFieldList(db.fluidMassDensity,
                        plotStyle="o",
                        winTitle="rho @ %g %i" % (control.time(), mpi.procs))
velPlot = plotFieldList(db.fluidVelocity,
                        yFunction = "%s.x",
                        plotStyle="o",
                        winTitle="vel @ %g %i" % (control.time(), mpi.procs))
mPlot = plotFieldList(db.fluidMass,
                        plotStyle="o",
                      winTitle="mass @ %g %i" % (control.time(), mpi.procs))
PPlot = plotFieldList(db.fluidPressure,
                        plotStyle="o",
                      winTitle="pressure @ %g %i" % (control.time(), mpi.procs))

s = damageModel.strain()
sl = ScalarFieldList1d()
sl.appendField(s)
sPlot = plotFieldList(sl, winTitle="strain @ %g %i" % (control.time(), mpi.procs))
