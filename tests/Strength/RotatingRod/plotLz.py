#-------------------------------------------------------------------------------
# A rod of stainless steel undergoing rotation.
#-------------------------------------------------------------------------------
from Numeric import *
from SolidSpheral import *
from SpheralTestUtilities import *
from SpheralController import *
from findLastRestart import *
from SpheralVisitDump import dumpPhysicsState
from identifyFragments import identifyFragments, fragmentProperties
from math import *

# Load the mpi module if we're parallel.
import loadmpi
mpi, procID, numProcs = loadmpi.loadmpi()

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("2-D rotating steel rod strength test")

#-------------------------------------------------------------------------------
# Generic problem parameters
# All CGS units.
#-------------------------------------------------------------------------------
SolidNodeListConstructor = SphNodeList2d
seed = 'lattice'

xlength, ylength = 20.0, 5.0
nx, ny = 100, 25
nPerh = 3.01

xmin = (-0.5*xlength, -0.5*ylength)
xmax = ( 0.5*xlength,  0.5*ylength)

rho0 = 7.9
dx = xlength/nx
dy = ylength/ny

# Initial rotational rate (radians/sec).  Basically rotating 0.1 times every 
# millisec, or 1e2 rps = 6000 rpm.
omega0 = 2.0*pi/1.0e-2

Qconstructor = MonaghanGingoldViscosity2d
Cl, Cq = 1.0, 1.0
Qlimiter = False
balsaraCorrection = True
epsilon2 = 1e-4
negligibleSoundSpeed = 1e-5
csMultiplier = 1e-4
HsmoothMin, HsmoothMax = 1e-5, 5.0
cfl = 0.5
useVelocityMagnitudeForDt = False
XSPH = False
SPHInternalVelocityGradient = True
epsilonTensile = 0.0
nTensile = 4
hybridMassDensityThreshold = 0.01

goalTime = 10.0*2.0*pi/omega0
dtSample = 0.005*goalTime
dt = 1e-10
dtMin, dtMax = 1e-12, 1e-5
dtGrowth = 2.0
maxSteps = None
statsStep = 10
smoothIters = 0
HEvolution = Hydro2d.HEvolutionType.IdealH
sumForMassDensity = Hydro2d.MassDensityType.IntegrateDensity # HybridDensity # CorrectedSumDensity

restartStep = 1000
dataDir = "dumps-RotatingSteelRod-2d-%ix%i" % (nx, ny)
restartDir = dataDir + "/restarts"
visitDir = dataDir + "/visit"
restartBaseName = restartDir + "/RotatingSteelRod-%ix%i" % (nx, ny)

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
eos = GruneisenEquationOfStateCGS2d(7.90,    # reference density  
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
strengthModel = SteinbergGuinanStrengthCGS2d(eos,
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
nodes = SolidNodeListConstructor("Stainless steel", eos) # , strengthModel)
nodeSet = [nodes]
for n in nodeSet:
    n.nodesPerSmoothingScale = nPerh
    n.epsilonTensile = epsilonTensile
    n.nTensile = nTensile
    n.XSPH = XSPH
    n.SPHInternalVelocityGradient = SPHInternalVelocityGradient
    output('n.name()')
    output('  n.nodesPerSmoothingScale')
    output('  n.epsilonTensile')
    output('  n.nTensile')
    output('  n.XSPH')
    output("  n.SPHInternalVelocityGradient")
del n

#-------------------------------------------------------------------------------
# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
#-------------------------------------------------------------------------------
WT = TableKernel2d(BSplineKernel2d(), 1000)
WTPi = TableKernel2d(BSplineKernel2d(), 1000)
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
    neighbor = TreeNeighbor2d(n,
                              kernelExtent = kernelExtent)
    n.registerNeighbor(neighbor)
    cache.append(neighbor)
neighborTimer.stop()
neighborTimer.printStatus()
del n

#-------------------------------------------------------------------------------
# Set node properties (positions, masses, H's, etc.)
#-------------------------------------------------------------------------------
eps0 = 0.0
if restoreCycle is None:
    print "Generating node distribution."
    from GenerateNodeDistribution2d import *
    from DistributeNodes import distributeNodes2d
    generator = GenerateNodeDistribution2d(nx,
                                           ny,
                                           rho0,
                                           seed,
                                           xmin = xmin,
                                           xmax = xmax,
                                           nNodePerh = nPerh)
    n = generator.globalNumNodes()
    distributeNodes2d([(nodes, generator)])
    output('mpi.reduce(nodes.numInternalNodes, mpi.MIN)')
    output('mpi.reduce(nodes.numInternalNodes, mpi.MAX)')
    output('mpi.reduce(nodes.numInternalNodes, mpi.SUM)')

    # Set node specific thermal energies
    eps0 = eos.specificThermalEnergy(rho0, 300.0)
    nodes.setSpecificThermalEnergy(ScalarField2d("tmp", nodes, eps0))

    # Set node velocites.
    for i in xrange(nodes.numInternalNodes):
        xi = nodes.positions()[i]
        r = xi.magnitude()
        runit = xi.unitVector()
        vunit = Vector2d(-runit.y, runit.x)
        nodes.velocity()[i] = vunit*r*omega0

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase2d()
for n in nodeSet:
    db.appendNodeList(n)
del n
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
hydro = Hydro2d(WT, WTPi, q)
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
strength = Strength2d()
output("strength")

#-------------------------------------------------------------------------------
# Construct a predictor corrector integrator.
#-------------------------------------------------------------------------------
integrator = PredictorCorrectorIntegrator2d(db)
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
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            initializeMassDensity = False)
output('control')

#-------------------------------------------------------------------------------
# Plot the result with strength.
#-------------------------------------------------------------------------------
dataDir = "dumps-RotatingSteelRod-2d-%ix%i-nph=%4.2f" % (nx, ny, nPerh)
restartDir = dataDir + "/restarts"
restartBaseName = restartDir + "/RotatingSteelRod-%ix%i" % (nx, ny)
restoreCycle = findLastRestart(restartBaseName)
control.setRestartBaseName(restartBaseName)
control.loadRestartFile(restoreCycle)

Lz0 = [x.z for x in control.conserve.amomHistory]

import Gnuplot
from SpheralGnuPlotUtilities import generateNewGnuPlot
d0 = Gnuplot.Data(control.conserve.timeHistory, Lz0,
                  title = "With strength",
                  with = "linespoints",
                  inline = True)
p = generateNewGnuPlot()
p.plot(d0)

## import pylab
## if mpi.rank == 0:
##     pylab.plot(control.conserve.timeHistory, Lz0,
##                "b-",
##                label = "With strength")

#-------------------------------------------------------------------------------
# Plot the result without strength.
#-------------------------------------------------------------------------------
dataDir = "dumps-RotatingSteelRod-2d-%ix%i-nostrength" % (nx, ny)
restartDir = dataDir + "/restarts"
restartBaseName = restartDir + "/RotatingSteelRod-%ix%i" % (nx, ny)
restoreCycle = findLastRestart(restartBaseName)
control.setRestartBaseName(restartBaseName)
control.loadRestartFile(restoreCycle)

Lz1 = [x.z for x in control.conserve.amomHistory]

d1 = Gnuplot.Data(control.conserve.timeHistory, Lz1,
                  title = "Without strength",
                  with = "linespoints",
                  inline = True)
p.replot(d1)
p("set key bottom left")

## if mpi.rank == 0:
##     pylab.plot(control.conserve.timeHistory, Lz1,
##                "r-",
##                label = "Without strength")
##     pylab.legend(loc = "best")
