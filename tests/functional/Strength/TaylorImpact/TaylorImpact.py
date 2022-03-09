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
# SPH 2D
#ATS:test(SELF, "--geometry 2d --crksph False --steps 100 --compatibleEnergy False --clearDirectories True --siloSnapShotFile Spheral_sph_2d_state_snapshot_1proc", np=1, label="Generate 1 proc SPH 2D reference data")
#ATS:test(SELF, "--geometry 2d --crksph False --steps 100 --compatibleEnergy False --clearDirectories True --siloSnapShotFile Spheral_sph_2d_state_snapshot_8proc", np=8, label="Generate 8 proc SPH 2D reference data")
#
# SPH RZ
#ATS:test(SELF, "--geometry RZ --crksph False --steps 100 --compatibleEnergy False --clearDirectories True --siloSnapShotFile Spheral_sph_rz_state_snapshot_1proc", np=1, label="Generate 1 proc SPH RZ reference data")
#ATS:test(SELF, "--geometry RZ --crksph False --steps 100 --compatibleEnergy False --clearDirectories True --siloSnapShotFile Spheral_sph_rz_state_snapshot_8proc", np=8, label="Generate 8 proc SPH RZ reference data")
#
# SPH 3D
#ATS:test(SELF, "--geometry 3d --crksph False --steps 100 --compatibleEnergy False --clearDirectories True --siloSnapShotFile Spheral_sph_3d_state_snapshot_8proc", np=8, label="Generate 8 proc SPH 3D reference data")
#
# SPH 2D (no grad h correction)
#ATS:test(SELF, "--geometry 2d --crksph False --steps 100 --compatibleEnergy False --clearDirectories True --gradhCorrection False --siloSnapShotFile Spheral_sph_nogradh_2d_state_snapshot_1proc", np=1, label="Generate 1 proc SPH 2D reference data (no grad h)")
#ATS:test(SELF, "--geometry 2d --crksph False --steps 100 --compatibleEnergy False --clearDirectories True --gradhCorrection False --siloSnapShotFile Spheral_sph_nogradh_2d_state_snapshot_8proc", np=8, label="Generate 8 proc SPH 2D reference data (no grad h)")
#
# SPH RZ (no grad h correction)
#ATS:test(SELF, "--geometry RZ --crksph False --steps 100 --compatibleEnergy False --clearDirectories True --gradhCorrection False --siloSnapShotFile Spheral_sph_nogradh_rz_state_snapshot_1proc", np=1, label="Generate 1 proc SPH RZ reference data (no grad h)")
#ATS:test(SELF, "--geometry RZ --crksph False --steps 100 --compatibleEnergy False --clearDirectories True --gradhCorrection False --siloSnapShotFile Spheral_sph_nogradh_rz_state_snapshot_8proc", np=8, label="Generate 8 proc SPH RZ reference data (no grad h)")
#
# SPH 3D (no grad h correction)
#ATS:test(SELF, "--geometry 3d --crksph False --steps 100 --compatibleEnergy False --clearDirectories True --gradhCorrection False --siloSnapShotFile Spheral_sph_nogradh_3d_state_snapshot_8proc", np=8, label="Generate 8 proc SPH 3D reference data (no grad h)")
#
# CRK 2D
#ATS:test(SELF, "--geometry 2d --crksph True  --steps 100 --compatibleEnergy False --densityUpdate RigorousSumDensity --clearDirectories True --siloSnapShotFile Spheral_crk_2d_state_snapshot_1proc", np=1, label="Generate 1 proc CRK 2D reference data")
#ATS:test(SELF, "--geometry 2d --crksph True  --steps 100 --compatibleEnergy False --densityUpdate RigorousSumDensity --clearDirectories True --siloSnapShotFile Spheral_crk_2d_state_snapshot_8proc", np=8, label="Generate 8 proc CRK 2D reference data")
#
# CRK RZ
#ATS:test(SELF, "--geometry RZ --crksph True  --steps 100 --compatibleEnergy False --densityUpdate RigorousSumDensity --clearDirectories True --siloSnapShotFile Spheral_crk_rz_state_snapshot_1proc", np=1, label="Generate 1 proc CRK RZ reference data")
#ATS:test(SELF, "--geometry RZ --crksph True  --steps 100 --compatibleEnergy False --densityUpdate RigorousSumDensity --clearDirectories True --siloSnapShotFile Spheral_crk_rz_state_snapshot_8proc", np=8, label="Generate 8 proc CRK RZ reference data")
#
# CRK 3D
#ATS:test(SELF, "--geometry 3d --crksph True  --steps 100 --compatibleEnergy False --densityUpdate RigorousSumDensity --clearDirectories True --siloSnapShotFile Spheral_crk_3d_state_snapshot_8proc", np=8, label="Generate 8 proc CRK 3D reference data")

import os, shutil
from math import *
import mpi

from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from SpheralController import *

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("Cu taylor anvil impact strength test")

#-------------------------------------------------------------------------------
# Generic problem parameters
# All (cm, gm, usec) units.
#-------------------------------------------------------------------------------
commandLine(geometry = "2d",         # one of (2d, 3d, RZ)

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
            crksph = False,
            asph = True,     # Only for H evolution, not hydro algorithm
            HUpdate = IdealH,
            densityUpdate = IntegrateDensity,
            compatibleEnergy = True,
            evolveTotalEnergy = False,
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
            correctionOrder = LinearOrder,                # CRKSPH
            volumeType = RKVoronoiVolume,                 # CRKSPH

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
            baseDir = "dumps-TaylorImpact",
            verbosedt = False,
            clearDirectories = False,

            # Should we generate a state snapshot on completion?
            siloSnapShotFile = "",
            )

assert geometry in ("2d", "3d", "RZ")
assert not (compatibleEnergy and evolveTotalEnergy)

exec("from Spheral%s import *" % geometry)

if crksph:
    hydroname = os.path.join("CRKSPH",
                             str(correctionOrder),
                             str(volumeType))
else:
    hydroname = os.path.join("SPH", "gradh=%s" % gradhCorrection)
if asph:
    hydroname = "A" + hydroname

# Restart and output files.
dataDir = os.path.join(baseDir,
                       geometry,
                       hydroname,
                       "XSPH=%s" % XSPH,
                       "reflect=%s" % reflect,
                       "%ix%i" % (nr, nz),
                       "procs=%i" % mpi.procs)
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
print "Generating node distribution."
#...............................................................................
# 2D
if geometry == "2d":
    from GenerateNodeDistribution2d import *
    from VoronoiDistributeNodes import distributeNodes2d as distributeNodes
    generator1 = GenerateNodeDistribution2d(nz, 2*nr, 
                                            rho = rho0,
                                            distributionType = "lattice",
                                            xmin = (0.0,    -rlength),
                                            xmax = (zlength, rlength),
                                            nNodePerh = nPerh,
                                            SPH = not asph)
    stuff2distribute = [(nodes1, generator1)]
    if not reflect:
        generator2 = GenerateNodeDistribution2d(nz, 2*nr,
                                                rho = rho0,
                                                distributionType = "lattice",
                                                xmin = (-zlength, -rlength),
                                                xmax = ( 0.0,      rlength),
                                                nNodePerh = nPerh,
                                                SPH = not asph)
        stuff2distribute.append((nodes2, generator2))

#...............................................................................
# RZ
elif geometry == "RZ":
    from GenerateNodeDistribution2d import *
    from VoronoiDistributeNodes import distributeNodes2d as distributeNodes
    generator1 = RZGenerator(GenerateNodeDistribution2d(nz, nr,
                                                        rho = rho0,
                                                        distributionType = "lattice",
                                                        xmin = (0.0,     0.0),
                                                        xmax = (zlength, rlength),
                                                        nNodePerh = nPerh,
                                                        SPH = not asph))
    stuff2distribute = [(nodes1, generator1)]
    if not reflect:
        generator2 = RZGenerator(GenerateNodeDistribution2d(nz, nr,
                                                            rho = rho0,
                                                            distributionType = "lattice",
                                                            xmin = (-zlength, 0.0),
                                                            xmax = ( 0.0,     rlength),
                                                            nNodePerh = nPerh,
                                                            SPH = not asph))
        stuff2distribute.append((nodes2, generator2))

#...............................................................................
# 3D
else:
    from GenerateNodeDistribution3d import *
    from VoronoiDistributeNodes import distributeNodes3d as distributeNodes
    rmin = 0.0
    rmax = rlength
    zmin = 0.0
    zmax = zlength
    generator1 = GenerateNodeDistribution3d(nr, nz, 0,
                                            rho = rho0,
                                            distributionType = "cylindrical",
                                            rmin = rmin,
                                            rmax = rmax,
                                            thetamin = 0.0,
                                            thetamax = 2.0*pi,
                                            zmin = zmin,
                                            zmax = zmax,
                                            nNodePerh = nPerh,
                                            SPH = not asph)
    stuff2distribute = [(nodes1, generator1)]
    if not reflect:
        generator2 = GenerateNodeDistribution3d(nr, nz, 0,
                                                rho = rho0,
                                                distributionType = "cylindrical",
                                                rmin = rmin,
                                                rmax = rmax,
                                                thetamin = 0.0,
                                                thetamax = 2.0*pi,
                                                zmin = -zmax,
                                                zmax = -zmin,
                                                nNodePerh = nPerh,
                                                SPH = not asph)
        stuff2distribute.append((nodes2, generator2))

#...............................................................................
distributeNodes(*tuple(stuff2distribute))
for n in nodeSet:
    output('n.name')
    output('   mpi.reduce(n.numInternalNodes, mpi.MIN)')
    output('   mpi.reduce(n.numInternalNodes, mpi.MAX)')
    output('   mpi.reduce(n.numInternalNodes, mpi.SUM)')
del n

#-------------------------------------------------------------------------------
# Set initial conditions
#-------------------------------------------------------------------------------
if geometry in ("2d", "RZ"):
    v0 = Vector(-vz0, 0.0)
else:
    v0 = Vector(0.0, 0.0, -vz0)
nodes1.specificThermalEnergy(ScalarField("tmp", nodes1, eps0))
nodes1.velocity(VectorField("tmp", nodes1, v0))
if not reflect:
    nodes2.specificThermalEnergy(ScalarField("tmp", nodes2, eps0))
    nodes2.velocity(VectorField("tmp", nodes2, -v0))

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
bcs = []
if reflect:
    if geometry == "3d":
        zplane = Plane(Vector(0.0, 0.0, 0.0), Vector(0.0, 0.0, 1.0))
        bc = ReflectingBoundary(zplane)
        bcs = [bc]
    else:
        xplane = Plane(Vector(0.0, 0.0), Vector(1.0, 0.0))
        bc = ReflectingBoundary(xplane)
        bcs = [bc]

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
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if crksph:
    hydro = CRKSPH(dataBase = db,
                   order = correctionOrder,
                   filter = filter,
                   cfl = cfl,
                   compatibleEnergyEvolution = compatibleEnergy,
                   evolveTotalEnergy = evolveTotalEnergy,
                   XSPH = XSPH,
                   densityUpdate = densityUpdate,
                   HUpdate = HUpdate,
                   ASPH = asph,
                   RZ = (geometry == "RZ"))
else:
    hydro = SPH(dataBase = db,
                W = WT,
                filter = filter,
                cfl = cfl,
                compatibleEnergyEvolution = compatibleEnergy,
                evolveTotalEnergy = evolveTotalEnergy,
                gradhCorrection = gradhCorrection,
                correctVelocityGradient = correctVelocityGradient,
                densityUpdate = densityUpdate,
                HUpdate = HUpdate,
                XSPH = XSPH,
                epsTensile = epsilonTensile,
                nTensile = nTensile,
                ASPH = asph,
                RZ = (geometry == "RZ"))

for bc in bcs:
    hydro.appendBoundary(bc)
output("hydro")
output("hydro.cfl")
output("hydro.useVelocityMagnitudeForDt")
output("hydro.HEvolution")
output("hydro.densityUpdate")
output("hydro.compatibleEnergyEvolution")

#-------------------------------------------------------------------------------
# Set the artificial viscosity parameters.
#-------------------------------------------------------------------------------
q = hydro.Q
if Cl:
    q.Cl = Cl
if Cq:
    q.Cq = Cq
if epsilon2:
    q.epsilon2 = epsilon2
if Qlimiter:
    q.limiter = Qlimiter
if balsaraCorrection:
    q.balsaraShearCorrection = balsaraCorrection
output("q")
output("q.Cl")
output("q.Cq")
output("q.epsilon2")
output("q.limiter")
output("q.balsaraShearCorrection")
try:
    output("q.linearInExpansion")
    output("q.quadraticInExpansion")
except:
    pass

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
from SpheralPointmeshSiloDump import dumpPhysicsState
control = SpheralController(integrator, WT,
                            volumeType = volumeType,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            redistributeStep = redistributeStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            vizMethod = dumpPhysicsState,
                            vizBaseName = "TaylorImpact",
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
    state0 = State(db, integrator.physicsPackages())
    state0.copyState()
    derivs = StateDerivatives(db, integrator.physicsPackages())
    derivs.Zero()
    integrator.preStepInitialize(state, derivs)
    dt = integrator.selectDt(dtmin, dtmax, state, derivs)
    integrator.initializeDerivatives(control.time() + dt, dt, state, derivs)
    integrator.evaluateDerivatives(control.time() + dt, dt, db, state, derivs)
    integrator.finalizeDerivatives(control.time() + dt, dt, db, state, derivs)

    # Grab the fields and their derivatives.
    mass = state0.scalarFields(HydroFieldNames.mass)
    rho = state0.scalarFields(HydroFieldNames.massDensity)
    pos = state0.vectorFields(HydroFieldNames.position)
    eps = state0.scalarFields(HydroFieldNames.specificThermalEnergy)
    vel = state0.vectorFields(HydroFieldNames.velocity)
    H = state0.symTensorFields(HydroFieldNames.H)
    P = state0.scalarFields(HydroFieldNames.pressure)
    S = state0.symTensorFields(SolidFieldNames.deviatoricStress)
    cs = state0.scalarFields(HydroFieldNames.soundSpeed)
    K = state0.scalarFields(SolidFieldNames.bulkModulus)
    mu = state0.scalarFields(SolidFieldNames.shearModulus)
    Y = state0.scalarFields(SolidFieldNames.yieldStrength)
    ps = state0.scalarFields(SolidFieldNames.plasticStrain)
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
