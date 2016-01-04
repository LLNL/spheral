#-------------------------------------------------------------------------------
# A rod of stainless steel undergoing rotation.
#-------------------------------------------------------------------------------
import mpi
import os
from SolidSpheral2d import *
from SpheralTestUtilities import *
from findLastRestart import *
from SpheralVisitDump import dumpPhysicsState
from identifyFragments import identifyFragments, fragmentProperties
from math import *

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("2-D rotating steel rod strength test")

#-------------------------------------------------------------------------------
# Generic problem parameters
# All CGS units.
#-------------------------------------------------------------------------------
commandLine(
    seed = "lattice",

    xlength = 3.0,
    ylength = 1.0,
    nx = 150,
    ny = 50,
    nPerh = 2.01,

    rho0 = 7.9,
    etamin = 0.6,
    etamax = 1.4,

    # Initial rotational rate (radians/sec).  Basically rotating once every
    # 0.01 seconds, or 1e2 rps = 6000 rpm.
    omega0 = 2.0*pi/1.0e-2,

    HydroConstructor = SolidASPHHydro,
    Qconstructor = MonaghanGingoldViscosity,
    Cl = 1.0,
    Cq = 1.0,
    Qlimiter = False,
    balsaraCorrection = True,
    epsilon2 = 1e-4,
    negligibleSoundSpeed = 1e-5,
    csMultiplier = 1e-4,
    hmin = 1e-5,
    hmax = 5.0,
    hminratio = 0.1,
    cfl = 0.5,
    useVelocityMagnitudeForDt = False,
    XSPH = False,
    epsilonTensile = 0.3,
    nTensile = 4,
    compatibleEnergyEvolution = True,
    gradhCorrection = True,
    correctVelocityGradient = True,

    goalTime = 0.01,
    steps = None,
    vizTime = 1e-5,
    vizStep = None,
    dt = 1e-10,
    dtMin = 1e-10,
    dtMax = 1e-5,
    dtGrowth = 2.0,
    maxSteps = None,
    statsStep = 10,
    smoothIters = 0,
    HEvolution = IdealH,
    densityUpdate = IntegrateDensity, # HybridDensity # CorrectedSumDensity

    restartStep = 1000,
    baseDir = "dumps-RotatingSteelRod-2d-%ix%i",
    )

dataDir = os.path.join(baseDir % (nx, ny),
                       "nPerh=%4.2f" % nPerh,
                       "correctVelocityGradient=%s" % correctVelocityGradient)
restartDir = os.path.join(dataDir, "restarts")
visitDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, "RotatingSteelRod-%ix%i" % (nx, ny))

xmin = (-0.5*xlength, -0.5*ylength)
xmax = ( 0.5*xlength,  0.5*ylength)

dx = xlength/nx
dy = ylength/ny

origin = Vector(-xlength, -ylength)

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
eos = GruneisenEquationOfStateCGS(rho0,    # reference density  
                                  etamin,  # etamin             
                                  etamax,  # etamax             
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
strengthModel = SteinbergGuinanStrengthCGS(eos,
                                           7.700000e11,        # G0
                                           2.2600e-12,         # A
                                           4.5500e-04,         # B
                                           3.4000e9,           # Y0
                                           2.5e10,             # Ymax
                                           1.0e-3,             # Yp
                                           43.0000,            # beta
                                           0.0,                # gamma0
                                           0.35,               # nhard
                                           coldFit,
                                           meltFit)

#-------------------------------------------------------------------------------
# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
#-------------------------------------------------------------------------------
WT = TableKernel(BSplineKernel(), 1000)
output("WT")
kernelExtent = WT.kernelExtent

#-------------------------------------------------------------------------------
# Create the NodeLists.
#-------------------------------------------------------------------------------
nodes = makeSolidNodeList("Stainless steel", eos, strengthModel, 
                          nPerh = nPerh,
                          hmin = hmin,
                          hmax = hmax,
                          hminratio = hminratio,
                          rhoMin = etamin*rho0,
                          rhoMax = etamax*rho0,
                          kernelExtent = kernelExtent,
                          xmin = Vector(-10*xlength, -10*ylength),  # Box size for neighbor selection
                          xmax = Vector( 10*xlength,  10*ylength))  # Box size for neighbor selection

output("nodes.name")
output("  nodes.nodesPerSmoothingScale")
output("  nodes.hmin")
output("  nodes.hmax")
output("  nodes.hminratio")
output("  nodes.rhoMin")
output("  nodes.rhoMax")

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
    distributeNodes2d((nodes, generator))
    output("mpi.reduce(nodes.numInternalNodes, mpi.MIN)")
    output("mpi.reduce(nodes.numInternalNodes, mpi.MAX)")
    output("mpi.reduce(nodes.numInternalNodes, mpi.SUM)")

    # Set node specific thermal energies
    nodes.specificThermalEnergy(ScalarField("tmp", nodes, eps0))

    # Set node velocites.
    for i in xrange(nodes.numInternalNodes):
        xi = nodes.positions()[i]
        r = xi.magnitude()
        runit = xi.unitVector()
        vunit = Vector(-runit.y, runit.x)
        nodes.velocity()[i] = vunit*r*omega0

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
db.appendNodeList(nodes)
output("db")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the artificial viscosities for the problem.
#-------------------------------------------------------------------------------
q = Qconstructor(Cl, Cq)
q.limiter = Qlimiter
q.balsaraShearCorrection = balsaraCorrection
q.epsilon2 = epsilon2
q.negligibleSoundSpeed = negligibleSoundSpeed
q.csMultiplier = csMultiplier
output("q")
output("q.Cl")
output("q.Cq")
output("q.limiter")
output("q.epsilon2")
output("q.negligibleSoundSpeed")
output("q.csMultiplier")
output("q.balsaraShearCorrection")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
hydro = HydroConstructor(W = WT,
                         Q = q,
                         cfl = cfl,
                         useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                         compatibleEnergyEvolution = compatibleEnergyEvolution,
                         gradhCorrection = False,
                         correctVelocityGradient = correctVelocityGradient,
                         densityUpdate = densityUpdate,
                         HUpdate = HEvolution,
                         XSPH = XSPH,
                         epsTensile = epsilonTensile,
                         nTensile = nTensile)
output("hydro")
output("hydro.cfl")
output("hydro.useVelocityMagnitudeForDt")
output("hydro.HEvolution")
output("hydro.sumForMassDensity")
output("hydro.compatibleEnergyEvolution")
output("hydro.gradhCorrection")
output("hydro.kernel()")
output("hydro.PiKernel()")

#-------------------------------------------------------------------------------
# Construct a predictor corrector integrator.
#-------------------------------------------------------------------------------
integrator = SynchronousRK2Integrator(db)
integrator.appendPhysicsPackage(hydro)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            restoreCycle = restoreCycle,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            vizDir = visitDir,
                            vizBaseName = "RotatingSteelRod-2d",
                            vizTime = vizTime,
                            vizStep = vizStep,
                            vizDerivs = True)
output("control")

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if steps is not None:
    control.step(steps)
else:
    control.advance(goalTime)
    control.conserve.writeHistory(os.path.join(dataDir, "conserve.txt"))
