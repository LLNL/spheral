#-------------------------------------------------------------------------------
# The Taylor anvil impact problem -- impact of a solid cylinder on an unyielding
# surface.
#
# This scenario is based on the v=205 m/sec example in
# Eakins & Thadhani, Journal of Applied Physics, 100, 073503 (2006)
#-------------------------------------------------------------------------------
#
# The following ATS setup is to generate reference data for the SpheralC tests.
#
#ATS:test(SELF, "--steps 100 --compatibleEnergy False --clearDirectories True --baseDir 2D_1proc_ref --siloSnapShotFile Spheral_state_snapshot_1proc", np=1, label="Generate 1 proc reference data")
#ATS:test(SELF, "--steps 100 --compatibleEnergy False --clearDirectories True --baseDir 2D_8proc_ref --siloSnapShotFile Spheral_state_snapshot_8proc", np=8, label="Generate 8 proc reference data")

import os, shutil
from math import *
import mpi

from SolidSpheral2d import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from SpheralController import *

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("2D Cu taylor anvil impact strength test")

#-------------------------------------------------------------------------------
# Generic problem parameters
# All (cm, gm, usec) units.
#-------------------------------------------------------------------------------
commandLine(seed = "lattice",

            # Geometry
            rlength = 0.945,
            zlength = 7.5,
            reflect = True,          # Use reflecting BC (True) or two rods (False)

            # Initial z velocity.
            vz0 = 2.05e-2,

            # Resolution
            nr = 10,
            nz = 80,
            nPerh = 2.01,
            
            # Material.
            etamin = 0.2,
            etamax = 4.0,

            # Artificial viscosity (and other numerical crap).
            CRKSPH = False,
            ASPH = True,     # Only for H evolution, not hydro algorithm
            Qconstructor = MonaghanGingoldViscosity,       # Artificial viscosity algorithm
            HUpdate = IdealH,
            densityUpdate = IntegrateDensity,
            compatibleEnergy = True,
            filter = 0.0,
            Cl = 1.0,                                      # Linear Q coefficient
            Cq = 1.0,                                      # Quadratic Q coefficient
            Qlimiter = False,                              # Q directional limiter switch
            balsaraCorrection = False,                     # Q shear switch
            epsilon2 = 1e-2,                               
            negligibleSoundSpeed = 1e-5,
            csMultiplier = 1e-4,
            hmin = 1e-5, 
            hmax = 1000.0, 
            hminratio = 0.1,
            limitIdealH = False,
            cfl = 0.4,
            useVelocityMagnitudeForDt = True,
            XSPH = True,
            epsilonTensile = 0.0,
            nTensile = 4,
            rigorousBoundaries = False,
            gradhCorrection = True,
            correctVelocityGradient = True,

            # Simulation control
            goalTime = 150.0,
            steps = None,
            dt = 1e-3,
            dtmin = 1.0e-3,
            dtmax = 100.0,
            dtGrowth = 2.0,
            maxSteps = None,
            statsStep = 10,
            restartStep = 1000,
            restoreCycle = -1,
            redistributeStep = 2000,
            vizCycle = 50,
            vizTime = 1.0,
            baseDir = "dumps-TaylorImpact-2d",
            verbosedt = False,
            clearDirectories = False,

            # Should we generate a state snapshot on completion?
            siloSnapShotFile = "",
            )

if CRKSPH:
    if ASPH:
        HydroConstructor = SolidACRKSPHHydro
    else:
        HydroConstructor = SolidCRKSPHHydro
    Qconstructor = CRKSPHMonaghanGingoldViscosity
else:
    if ASPH:
        HydroConstructor = SolidASPHHydro
    else:
        HydroConstructor = SolidSPHHydro

# Restart and output files.
dataDir = os.path.join(baseDir,
                       "reflect=%s" % reflect,
                       "%ix%i" % (nr, nz),
                       "XSPH=%s" % XSPH,
                       str(HydroConstructor).split("'")[1].split(".")[-1],
                       str(Qconstructor).split("'")[1].split(".")[-1])
restartDir = os.path.join(dataDir, "restarts", "proc-%04i" % mpi.rank)
vizDir = os.path.join(dataDir, "viz")
restartBaseName = os.path.join(restartDir, "TaylorImpact-%i-%i" % (nr, nz))

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if clearDirectories and os.path.exists(dataDir):
        shutil.rmtree(dataDir)
    if not os.path.exists(dataDir):
        os.makedirs(dataDir)
    if not os.path.exists(vizDir):
        os.makedirs(vizDir)
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
mpi.barrier()
if not os.path.exists(restartDir):
    os.makedirs(restartDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# Construct our base units.
#-------------------------------------------------------------------------------
units = PhysicalConstants(0.01,    # Unit length (m)
                          0.001,   # Unit mass (kg)
                          1.0e-6)  # Unit time (sec)

#-------------------------------------------------------------------------------
# Copper material parameters.
#-------------------------------------------------------------------------------
eosCu = TillotsonEquationOfState("copper",
                                 etamin,
                                 etamax,
                                 units)
rho0 = eosCu.referenceDensity
eps0 = 0.0

coldFitCu = NinthOrderPolynomialFit(-1.05111874e-02,
                                    -2.13429672e-02,
                                     6.92768584e-01,
                                    -2.45626513e-02,
                                    -2.48677403e-02,
                                     4.35373677e-02,
                                     0.00000000e+00,
                                     0.00000000e+00,
                                     0.00000000e+00,
                                     0.00000000e+00)

meltFitCu = NinthOrderPolynomialFit( 5.22055639e-02,
                                     1.90143176e-01,
                                     8.51351901e-01,
                                    -1.12049022e-01,
                                    -6.11436674e-03,
                                     4.36007831e-02,
                                     0.00000000e+00,
                                     0.00000000e+00,
                                     0.00000000e+00,
                                     0.00000000e+00)

strengthCu = SteinbergGuinanStrengthMKS(eosCu,
                                        4.770000e-01,   # G0 (Mb)
                                        2.8300e+00,     # A  (Mb^-1)
                                        3.7700e-04,     # B  (dimensionless)
                                        1.2000e-03,     # Y0 (Mb)
                                        6.4000e-03,     # Ymax (Mb)
                                        1.0000e-03,     # Yp (dimensionless)
                                        3.6000e+01,     # beta (dimensionless)
                                        0.0000e+00,     # gamma0 (dimensionless)
                                        4.5000e-01,     # nhard (dimensionless)
                                        coldFitCu,
                                        meltFitCu)

#-------------------------------------------------------------------------------
# Create the NodeLists.
#-------------------------------------------------------------------------------
nodes1 = makeSolidNodeList("Cylinder 1", eosCu, strengthCu,
                           nPerh = nPerh,
                           hmin = hmin,
                           hmax = hmax,
                           rhoMin = etamin*rho0,
                           rhoMax = etamax*rho0,
                           xmin = Vector(-10, -10, -10),          # Box size for neighbor finding
                           xmax = Vector( 10,  10,  10))          # Box size for neighbor finding
nodeSet = [nodes1]
if not reflect:
    nodes2 = makeSolidNodeList("Cylinder 2", eosCu, strengthCu,
                               nPerh = nPerh,
                               hmin = hmin,
                               hmax = hmax,
                               rhoMin = etamin*rho0,
                               rhoMax = etamax*rho0,
                               xmin = Vector(-10, -10, -10),          # Box size for neighbor finding
                               xmax = Vector( 10,  10,  10))          # Box size for neighbor finding
    nodeSet.append(nodes2)

for n in nodeSet:
    output("n.name")
    output("  n.nodesPerSmoothingScale")
    output("  n.hmin")
    output("  n.hmax")
    output("  n.rhoMin")
    output("  n.rhoMax")
del n

#-------------------------------------------------------------------------------
# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
#-------------------------------------------------------------------------------
WT = TableKernel(NBSplineKernel(3), 1000)
output('WT')

#-------------------------------------------------------------------------------
# Set node properties (positions, masses, H's, etc.)
#-------------------------------------------------------------------------------
from GenerateNodeDistribution2d import *
from VoronoiDistributeNodes import distributeNodes2d
print "Generating node distribution."
generator1 = GenerateNodeDistribution2d(2*nr, nz, 
                                        rho = rho0,
                                        distributionType = seed,
                                        xmin = (-rlength, 0.0),
                                        xmax = ( rlength, zlength),
                                        nNodePerh = nPerh)
stuff2distribute = [(nodes1, generator1)]
if not reflect:
    generator2 = GenerateNodeDistribution2d(2*nr, nz,
                                            rho = rho0,
                                            distributionType = seed,
                                            xmin = (-rlength, -zlength),
                                            xmax = ( rlength,  0.0),
                                            nNodePerh = nPerh)
    stuff2distribute.append((nodes2, generator2))
distributeNodes2d(*tuple(stuff2distribute))
for n in nodeSet:
    output('n.name')
    output('   mpi.reduce(n.numInternalNodes, mpi.MIN)')
    output('   mpi.reduce(n.numInternalNodes, mpi.MAX)')
    output('   mpi.reduce(n.numInternalNodes, mpi.SUM)')
del n

nodes1.specificThermalEnergy(ScalarField("tmp", nodes1, eps0))
nodes1.velocity(VectorField("tmp", nodes1, Vector(0.0, -vz0)))
if not reflect:
    nodes2.specificThermalEnergy(ScalarField("tmp", nodes2, eps0))
    nodes2.velocity(VectorField("tmp", nodes2, Vector(0.0, vz0)))

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
bcs = []
if reflect:
    yplane = Plane(Vector(0.0, 0.0), Vector(0.0, 1.0))
    bc = ReflectingBoundary(yplane)
    bcs.append(bc)

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
for n in nodeSet:
    db.appendNodeList(n)
output('db')
output('db.numNodeLists')
output('db.numFluidNodeLists')

#-------------------------------------------------------------------------------
# Construct the artificial viscosity.
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
if CRKSPH:
    hydro = HydroConstructor(W = WT,
                             Q = q,
                             filter = filter,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             XSPH = XSPH,
                             densityUpdate = densityUpdate,
                             HUpdate = HUpdate)
else:
    hydro = HydroConstructor(W = WT,
                             Q = q,
                             filter = filter,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             gradhCorrection = gradhCorrection,
                             correctVelocityGradient = correctVelocityGradient,
                             densityUpdate = densityUpdate,
                             HUpdate = HUpdate,
                             XSPH = XSPH,
                             epsTensile = epsilonTensile,
                             nTensile = nTensile)
for bc in bcs:
    hydro.appendBoundary(bc)
output("hydro")
output("hydro.cfl")
output("hydro.useVelocityMagnitudeForDt")
output("hydro.HEvolution")
output("hydro.densityUpdate")
output("hydro.compatibleEnergyEvolution")
output("hydro.kernel()")
output("hydro.PiKernel()")

#-------------------------------------------------------------------------------
# Construct a time integrator.
#-------------------------------------------------------------------------------
integrator = CheapSynchronousRK2Integrator(db)
integrator.appendPhysicsPackage(hydro)
integrator.lastDt = dt
integrator.verbose = verbosedt
if dtmin:
    integrator.dtMin = dtmin
if dtmax:
    integrator.dtMax = dtmax
integrator.dtGrowth = dtGrowth
integrator.rigorousBoundaries = rigorousBoundaries
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.rigorousBoundaries")
output("integrator.verbose")

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            redistributeStep = redistributeStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            vizBaseName = "TaylorImpact-2d",
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime)
output("control")

#-------------------------------------------------------------------------------
# In the two material case, it's useful to smooth the initial velocity field
# in order to avoid interpenetration at the interface.
#-------------------------------------------------------------------------------
if (not reflect) and control.totalSteps == 0:
    print "Smoothing initial velocity field."
    state = State(db, integrator.physicsPackages())
    derivs = StateDerivatives(db, integrator.physicsPackages())
    integrator.initialize(state, derivs)
    vel = db.fluidVelocity
    pos = db.fluidPosition
    H = db.fluidHfield
    m = db.fluidMass
    velsmooth = smoothVectorFieldsMash(vel, pos, m, H, WT)
    vel.assignFields(velsmooth)
    control.dropViz(control.totalSteps, 0.0, 0.0)

#-------------------------------------------------------------------------------
# Advance to completetion.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)
else:
    control.advance(goalTime, maxSteps)

#-------------------------------------------------------------------------------
# If requested, generate table output of the full results.
#-------------------------------------------------------------------------------
if siloSnapShotFile:
    from siloPointmeshDump import siloPointmeshDump
    print "Generating snapshot in silo files."

    # First generate the state and derivatives.
    state = State(db, integrator.physicsPackages())
    derivs = StateDerivatives(db, integrator.physicsPackages())
    derivs.Zero()
    integrator.preStepInitialize(state, derivs)
    dt = integrator.selectDt(dtmin, dtmax, state, derivs)
    integrator.initializeDerivatives(control.time() + dt, dt, state, derivs)
    integrator.evaluateDerivatives(control.time() + dt, dt, db, state, derivs)
    integrator.finalizeDerivatives(control.time() + dt, dt, db, state, derivs)

    # Grab the fields and their derivatives.
    mass = state.scalarFields(HydroFieldNames.mass)
    rho = state.scalarFields(HydroFieldNames.massDensity)
    pos = state.vectorFields(HydroFieldNames.position)
    eps = state.scalarFields(HydroFieldNames.specificThermalEnergy)
    vel = state.vectorFields(HydroFieldNames.velocity)
    H = state.symTensorFields(HydroFieldNames.H)
    P = state.scalarFields(HydroFieldNames.pressure)
    S = state.symTensorFields(SolidFieldNames.deviatoricStress)
    cs = state.scalarFields(HydroFieldNames.soundSpeed)
    K = state.scalarFields(SolidFieldNames.bulkModulus)
    mu = state.scalarFields(SolidFieldNames.shearModulus)
    Y = state.scalarFields(SolidFieldNames.yieldStrength)
    ps = state.scalarFields(SolidFieldNames.plasticStrain)
    massSum = derivs.scalarFields("new " + HydroFieldNames.massDensity)
    DrhoDt = derivs.scalarFields("delta " + HydroFieldNames.massDensity)
    DvelDt = derivs.vectorFields(HydroFieldNames.hydroAcceleration)
    DepsDt = derivs.scalarFields("delta " + HydroFieldNames.specificThermalEnergy)
    DvelDx = derivs.tensorFields(HydroFieldNames.velocityGradient)
    DHDt = derivs.symTensorFields("delta " + HydroFieldNames.H)
    Hideal = derivs.symTensorFields("new " + HydroFieldNames.H)
    DSDt = derivs.symTensorFields("delta " + SolidFieldNames.deviatoricStress)

    # Write the sucker.
    siloPointmeshDump(siloSnapShotFile, 
                      fieldLists = [mass, rho, pos, eps, vel, H, P, S, cs, K, mu, Y, ps,
                                    massSum, DrhoDt, DvelDt, DepsDt, DvelDx, DHDt, Hideal, DSDt],
                      baseDirectory = dataDir,
                      label = "Spheral++ snapshot of state and derivatives.",
                      time = control.time(),
                      cycle = control.totalSteps)
