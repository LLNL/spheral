#-------------------------------------------------------------------------------
# The Taylor anvil impact problem -- impact of a solid cylinder on an unyielding
# surface.
#-------------------------------------------------------------------------------
from math import *
import mpi

from SolidSpheral3d import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from SpheralController import *
from findLastRestart import *

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("3D Cu taylor anvil impact strength test")

#-------------------------------------------------------------------------------
# Generic problem parameters
# All (cm, gm, usec) units.
#-------------------------------------------------------------------------------
units = PhysicalConstants

commandLine(seed = "cylindrical",

            # Geometry
            rlength = 0.25,
            zlength = 2.0,
            reflect = True,                                # Use reflecting BC

            # Initial z velocity.
            vz0 = -2.5e-2,

            # Resolution
            nr = 10,
            nz = 80,
            nPerh = 2.01,
            
            # Material.
            etamin = 0.2,
            etamax = 4.0,

            # Artificial viscosity (and other numerical crap).
            HydroConstructor = SolidASPHHydro,             # Hydro algorithm
            Qconstructor = MonaghanGingoldViscosity,       # Artificial viscosity algorithm
            HEvolution = IdealH,
            densityUpdate = IntegrateDensity,
            compatibleEnergyEvolution = True,
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
            restoreCycle = None,
            redistributeStep = 2000,
            vizCycle = 50,
            vizTime = 1.0,
            baseDir = "dumps-TaylorImpact-3d",
            verbosedt = False,
            )

rmin = 0.0
rmax = rlength
zmin = 0.0
zmax = zlength

# Restart and output files.
dataDir = os.path.join(baseDir,
                       "reflect=%s" % reflect,
                       "%ix%i" % (nr, nz))
restartDir = os.path.join(dataDir, "restarts", "proc-%04i" % mpi.rank)
vizDir = os.path.join(dataDir, "viz")
restartBaseName = os.path.join(restartDir, "TaylorImpact-%i-%i" % (nr, nz))

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
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
# If we're restarting, find the set of most recent restart files.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    restoreCycle = findLastRestart(restartBaseName)

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
    nodes2 = makeSolidNodeList("Cylinder 1", eosCu, strengthCu,
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
WT = TableKernel(BSplineKernel(), 1000)
WTPi = WT
output('WT')
output('WTPi')

#-------------------------------------------------------------------------------
# Set node properties (positions, masses, H's, etc.)
#-------------------------------------------------------------------------------
if restoreCycle is None:
    from GenerateNodeDistribution3d import *
    from VoronoiDistributeNodes import distributeNodes3d
    print "Generating node distribution."
    generator1 = GenerateNodeDistribution3d(nr, nz, 0,
                                            rho = rho0,
                                            distributionType = seed,
                                            rmin = rmin,
                                            rmax = rmax,
                                            thetamin = 0.0,
                                            thetamax = 2.0*pi,
                                            zmin = zmin,
                                            zmax = zmax,
                                            nNodePerh = nPerh)
    stuff2distribute = [(nodes1, generator1)]
    if not reflect:
        generator2 = GenerateNodeDistribution3d(nr, nz, 0,
                                                rho = rho0,
                                                distributionType = seed,
                                                rmin = rmin,
                                                rmax = rmax,
                                                thetamin = 0.0,
                                                thetamax = 2.0*pi,
                                                zmin = -zmax,
                                                zmax = -zmin,
                                                nNodePerh = nPerh)
        stuff2distribute.append((nodes2, generator2))
    distributeNodes3d(*tuple(stuff2distribute))
    for n in nodeSet:
        output('n.name')
        output('   mpi.reduce(n.numInternalNodes, mpi.MIN)')
        output('   mpi.reduce(n.numInternalNodes, mpi.MAX)')
        output('   mpi.reduce(n.numInternalNodes, mpi.SUM)')
    del n

    nodes1.specificThermalEnergy(ScalarField("tmp", nodes1, eps0))
    nodes1.velocity(VectorField("tmp", nodes1, Vector(0.0, 0.0, vz0)))
    if not reflect:
        nodes2.specificThermalEnergy(ScalarField("tmp", nodes2, eps0))
        nodes2.velocity(VectorField("tmp", nodes2, Vector(0.0, 0.0, -vz0)))

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
bcs = []
if reflect:
    zplane = Plane(Vector(0.0, 0.0, 0.0), Vector(0.0, 0.0, 1.0))
    bc = ReflectingBoundary(zplane)
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
hydro = HydroConstructor(WT,
                         WTPi,
                         q,
                         cfl = cfl,
                         useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                         compatibleEnergyEvolution = compatibleEnergyEvolution,
                         gradhCorrection = False,
                         densityUpdate = densityUpdate,
                         HUpdate = HEvolution,
                         XSPH = XSPH,
                         epsTensile = epsilonTensile,
                         nTensile = nTensile)
for bc in bcs:
    hydro.appendBoundary(bc)
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
from SpheralVisitDump import dumpPhysicsState
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            redistributeStep = redistributeStep,
                            restartBaseName = restartBaseName,
                            vizBaseName = "RoundAsteroid-KineticImpactor-3d",
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime,
                            vizMethod = dumpPhysicsState)
output("control")

#-------------------------------------------------------------------------------
# Advance to completetion.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)
    raise ValueError, ("Completed %i steps." % steps)

else:
    control.advance(goalTime, maxSteps)
