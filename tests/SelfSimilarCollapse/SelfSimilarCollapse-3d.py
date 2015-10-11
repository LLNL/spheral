#-------------------------------------------------------------------------------
# Spheral++ script to model the gravitationally driven collapse of an isothermal
# sphere of gas.  The fluid is assumed to obey a polytropic equation of state:
#           P = K \rho^\gamma
# and there is an analytic self-similar solution for this case.
# 
# See P. Blottiau, S. Bouquet, & J.P. Cheize, Astron. Astrophs., 207, 24-36 (1988)
#-------------------------------------------------------------------------------
import os, sys, shutil, mpi
from math import *
from Spheral3d import *
from SpheralTestUtilities import *
from findLastRestart import *
from GenerateNodeDistribution3d import *

from VoronoiDistributeNodes import distributeNodes3d as distributeNodes

title("3-D gravitational hydro test -- self-similar collapse of an isothermal sphere.")

#-------------------------------------------------------------------------------
# Problem parameters.
#-------------------------------------------------------------------------------
commandLine(nx = 100,         # Number of points across diameter of sphere.
            M0 = 1.0,         # Initial mass of sphere (solar masses)
            R0 = 1.0,         # Initial radius (AU)

            # Material properties.
            K = 1.0,          # Polytropic constant
            gamma = 1.0,      # Polytropic index (solution applies to \gamma \in [1, 4/3]
            mu = 2.0,         # Molecular weight

            # Gravity parameters.
            plummerLength = 1.0e-2, # (AU) Plummer softening scale
            opening = 1.0,          # (dimensionless, OctTreeGravity) opening parameter for tree walk
            fdt = 0.1,              # (dimensionless, OctTreeGravity) timestep multiplier

            # Hydro parameters.
            nPerh = 2.01,
            HydroConstructor = SPHHydro,
            Qconstructor = MonaghanGingoldViscosity,
            Cl = 1.0, 
            Cq = 0.75,
            Qlimiter = False,
            balsaraCorrection = False,
            hmin = 1.0e-10,
            hmax = 10.0,
            hminratio = 0.1,
            cfl = 0.5,
            XSPH = True,
            epsilonTensile = 0.0,
            nTensile = 4.0,
            densityUpdate = RigorousSumDensity,
            compatibleEnergy = False,
            gradhCorrection = False,

            # Timestepping/advancement
            IntegratorConstructor = CheapSynchronousRK2Integrator,
            goalTime = 1.0,
            steps = None,
            vizCycle = 20,
            vizTime = 0.01,
            dt = 1.0e-5,
            dtMin = 1.0e-8, 
            dtMax = 1.0e-2,
            dtGrowth = 2.0,
            maxSteps = None,
            statsStep = 10,
            smoothIters = 0,
            HEvolution = IdealH,
            domainIndependent = False,
            rigorousBoundaries = False,
            dtverbose = False,

            clearDirectories = False,
            restoreCycle = None,
            restartStep = 50,
            checkRestart = False,
            dataDir = "dumps-SphericalCollapse-3d",
            vizName = "IsothermalSphericalCollapse",
            )

G = Solar().G

restartDir = os.path.join(dataDir, "restarts")
vizDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, vizName)

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
if mpi.rank == 0:
    for DIR in (restartDir, vizDir):
        if clearDirectories and os.path.exists(DIR):
            shutil.rmtree(DIR)
        if not os.path.exists(DIR):
            os.makedirs(DIR)
mpi.barrier()

#-------------------------------------------------------------------------------
# If we're restarting, find the set of most recent restart files.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    restoreCycle = findLastRestart(restartBaseName)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(BSplineKernel(), 100000)
output("WT")

#===============================================================================
# Material properties
#===============================================================================
eos = PolytropicEquationOfStateSolar3d(K, gamma, mu)

#===============================================================================
# Create the NodeLists.
#===============================================================================
nodes = makeFluidNodeList("gas", eos,
                          hmin = hmin,
                          hmax = hmax,
                          hminratio = hminratio,
                          nPerh = nPerh)

#===============================================================================
# Create the generators.
#===============================================================================
if restoreCycle is None:
    rho0 = M0/(4.0/3.0*pi*R0**3)
    generator = GenerateNodeDistribution3d(nx, nx, nx,
                                           rho = rho0,
                                           distributionType = "lattice",
                                           xmin = (-R0, -R0, -R0),
                                           xmax = ( R0,  R0,  R0),
                                           rmin = 0.0,
                                           rmax = R0,
                                           nNodePerh = nPerh,
                                           SPH = (HydroConstructor == SPHHydro))

    distributeNodes((nodes, generator))
    print "Num internal nodes for ", nodes.name, " : ", mpi.allreduce(nodes.numInternalNodes, mpi.SUM)

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
output("db.appendNodeList(nodes)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the artificial viscosity.
#-------------------------------------------------------------------------------
q = Qconstructor(Cl, Cq)
q.limiter = Qlimiter
q.balsaraShearCorrection = balsaraCorrection
output("q")
output("q.Cl")
output("q.Cq")
output("q.epsilon2")
output("q.limiter")
output("q.balsaraShearCorrection")

#-------------------------------------------------------------------------------
# Gravity baby!
#-------------------------------------------------------------------------------
gravity = OctTreeGravity(G = G,
                         softeningLength = plummerLength,
                         opening = opening,
                         ftimestep = fdt)

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
hydro = HydroConstructor(W = WT,
                         Q = q,
                         cfl = cfl,
                         compatibleEnergyEvolution = compatibleEnergy,
                         gradhCorrection = gradhCorrection,
                         XSPH = XSPH,
                         densityUpdate = densityUpdate,
                         HUpdate = HEvolution,
                         epsTensile = epsilonTensile,
                         nTensile = nTensile)
output("hydro")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.gradhCorrection")
output("hydro.XSPH")
output("hydro.densityUpdate")
output("hydro.HEvolution")
output("hydro.epsilonTensile")
output("hydro.nTensile")

#-------------------------------------------------------------------------------
# Construct a time integrator, and add the physics packages.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
integrator.appendPhysicsPackage(gravity)
integrator.appendPhysicsPackage(hydro)
integrator.lastDt = dt
integrator.dtMin = dtMin
integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.domainDecompositionIndependent = domainIndependent
integrator.verbose = dtverbose
integrator.rigorousBoundaries = rigorousBoundaries
output("integrator")
output("integrator.havePhysicsPackage(gravity)")
output("integrator.havePhysicsPackage(hydro)")
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
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            vizBaseName = vizName,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime)
output("control")

# Smooth the initial conditions.
if restoreCycle is None:
    control.iterateIdealH(hydro)
    control.smoothState(smoothIters)
    if densityUpdate in (VoronoiCellDensity, SumVoronoiCellDensity):
        print "Reinitializing node masses."
        control.voronoiInitializeMass()
    control.dropRestartFile()
    control.dropViz()

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)

else:
    control.advance(goalTime, maxSteps)
    control.updateViz(control.totalSteps, integrator.currentTime, 0.0)
    control.dropRestartFile()
