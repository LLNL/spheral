#-------------------------------------------------------------------------------
# This is our shearing variant of the planar Noh problem.
#-------------------------------------------------------------------------------
from math import *
import shutil
import mpi
from Spheral2d import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from findLastRestart import *
from SpheralVisitDump import dumpPhysicsState

from GenerateNodeDistribution2d import *

title("2-D integrated hydro test -- Shearing Planar Noh problem")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(seed = "lattice",

            nx = 20,
            ny = 100,
            nPerh = 2.01,

            x0 = 0.0,
            x1 = 0.2,
            y0 = 0.0,
            y1 = 1.0,

            rho1 = 1.0,
            eps1 = 0.0,
            vshear = 1.0,
            vy = -1.0,

            gamma = 5.0/3.0,
            mu = 1.0,

            HydroConstructor = SPHHydro,
            Qconstructor = MonaghanGingoldViscosity,
            #Qconstructor = TensorMonaghanGingoldViscosity,
            Cl = 1.0, 
            Cq = 0.75,
            Qlimiter = False,
            balsaraCorrection = False,
            epsilon2 = 1e-2,
            hmin = 0.0001, 
            hmax = 0.5,
            hminratio = 0.1,
            cfl = 0.5,
            XSPH = True,
            epsilonTensile = 0.0,
            nTensile = 8,
            hourglass = None,
            hourglassOrder = 0,
            hourglassLimiter = 0,
            hourglassFraction = 0.5,

            IntegratorConstructor = CheapSynchronousRK2Integrator,
            goalTime = 0.6,
            steps = None,
            vizCycle = None,
            vizTime = 0.05,
            dt = 0.0001,
            dtMin = 1.0e-5, 
            dtMax = 0.1,
            dtGrowth = 2.0,
            maxSteps = None,
            statsStep = 10,
            smoothIters = 0,
            HEvolution = IdealH,
            domainIndependent = False,
            rigorousBoundaries = True,
            dtverbose = False,

            densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
            compatibleEnergy = False,
            gradhCorrection = False,

            clearDirectories = False,
            restoreCycle = None,
            restartStep = 1000,
            dataRoot = "dumps-shearingNoh-2d",

            graphics = True,
            )

dataDir = os.path.join(dataRoot,
                       str(HydroConstructor).split()[1].replace(">", "").replace("'",""),
                       str(Qconstructor).split()[1].replace(">", "").replace("'",""),
                       "basaraShearCorrection=%s_Qlimiter=%s" % (balsaraCorrection, Qlimiter),
                       "nperh=%4.2f" % nPerh,
                       "XSPH=%s" % XSPH,
                       "densityUpdate=%s" % densityUpdate,
                       "compatibleEnergy=%s" % compatibleEnergy,
                       "gradhCorrection=%s" % gradhCorrection,
                       "nx=%i_ny=%i" % (nx, ny))
restartDir = os.path.join(dataDir, "restarts")
vizDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, "Noh-shear-2d-%ix%i" % (nx, ny))

xmin = (x0, y0)
xmax = (x1, y1)

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
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
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS(gamma, mu)

#-------------------------------------------------------------------------------
# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
#-------------------------------------------------------------------------------
WT = TableKernel(BSplineKernel(), 1000)
WTPi = WT
output("WT")
output("WTPi")

#-------------------------------------------------------------------------------
# Create a NodeList and associated Neighbor object.
#-------------------------------------------------------------------------------
nodes1 = makeFluidNodeList("nodes1", eos, 
                           hmin = hmin,
                           hmax = hmax,
                           nPerh = nPerh)

#-------------------------------------------------------------------------------
# Set node properties.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    from DistributeNodes import distributeNodes2d
    print "Generating node distribution."
    generator1 = GenerateNodeDistribution2d(nx, ny, rho1, seed,
                                            xmin = xmin,
                                            xmax = xmax,
                                            nNodePerh = nPerh)
    n1 = generator1.globalNumNodes()

    if mpi.procs > 1:
        from VoronoiDistributeNodes import distributeNodes2d
        #from PeanoHilbertDistributeNodes import distributeNodes2d
    else:
        from DistributeNodes import distributeNodes2d
    distributeNodes2d((nodes1, generator1))
    output("mpi.reduce(nodes1.numInternalNodes, mpi.MIN)")
    output("mpi.reduce(nodes1.numInternalNodes, mpi.MAX)")
    output("mpi.reduce(nodes1.numInternalNodes, mpi.SUM)")

    # Set node specific thermal energies
    nodes1.specificThermalEnergy(ScalarField("tmp", nodes1, eps1))

    # Set node velocities
    pos = nodes1.positions()
    vel = nodes1.velocity()
    for i in xrange(nodes1.numInternalNodes):
        vel[i] = Vector(vshear*cos(2.0*pi*pos[i].y), vy)

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
output("db.appendNodeList(nodes1)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the artificial viscosity.
#-------------------------------------------------------------------------------
q = Qconstructor(Cl, Cq)
q.epsilon2 = epsilon2
q.limiter = Qlimiter
q.balsaraShearCorrection = balsaraCorrection
output("q")
output("q.Cl")
output("q.Cq")
output("q.epsilon2")
output("q.limiter")
output("q.balsaraShearCorrection")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
hydro = HydroConstructor(WT,
                         WTPi,
                         q,
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

packages = [hydro]

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane0 = Plane(Vector(*xmin), Vector( 1.0, 0.0))
xPlane1 = Plane(Vector(*xmax), Vector(-1.0, 0.0))
yPlane0 = Plane(Vector(*xmin), Vector( 0.0, 1.0))
yPlane1 = Plane(Vector(*xmax), Vector( 0.0,-1.0))
xbc = PeriodicBoundary(xPlane0, xPlane1)
ybc0 = ReflectingBoundary(yPlane0)
ybc1 = ReflectingBoundary(yPlane1)
for p in packages:
    for bc in (xbc, ybc0, ybc1):
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
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
if hourglass:
    output("integrator.havePhysicsPackage(hg)")
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
                            vizBaseName = "Noh-shear-2d-%ix%i" % (nx, ny),
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime)
output("control")

# Smooth the initial conditions.
if restoreCycle is not None:
    control.loadRestartFile(restoreCycle)
else:
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

#-------------------------------------------------------------------------------
# Plot the results.
#-------------------------------------------------------------------------------
if graphics:

    # Plot the node positions.
    import Gnuplot
    rPlot = plotNodePositions2d(db, colorNodeLists=0, colorDomains=1)

    # Plot the final state.
    rhoPlot, vrPlot, epsPlot, PPlot, HPlot = plotState(db,
                                                       xFunction = "%s.y",
                                                       vecyFunction = "%s.y",
                                                       tenyFunction = "%s.yy ** -1")

    # Overplot the analytic solution.
    import NohAnalyticSolution
    answer = NohAnalyticSolution.NohSolution(1,
                                             h0 = nPerh*y1/ny)
    plotAnswer(answer, control.time(),
               rhoPlot = rhoPlot,
               velPlot = vrPlot,
               epsPlot = epsPlot,
               PPlot = PPlot,
               HPlot = HPlot)

    # Report the error norms.
    rmin, rmax = 0.05, 0.35
    r = mpi.allreduce([x.y for x in nodes1.positions().internalValues()], mpi.SUM)
    rho = mpi.allreduce(list(nodes1.massDensity().internalValues()), mpi.SUM)
    v = mpi.allreduce([x.y for x in nodes1.velocity().internalValues()], mpi.SUM)
    eps = mpi.allreduce(list(nodes1.specificThermalEnergy().internalValues()), mpi.SUM)
    Pf = ScalarField("pressure", nodes1)
    nodes1.pressure(Pf)
    P = mpi.allreduce(list(Pf.internalValues()), mpi.SUM)
    if mpi.rank == 0:
        from SpheralGnuPlotUtilities import multiSort
        import Pnorm
        multiSort(r, rho, v, eps, P)
        rans, vans, epsans, rhoans, Pans, hans = answer.solution(control.time(), r)
        print "\tQuantity \t\tL1 \t\t\tL2 \t\t\tLinf"
        for (name, data, ans) in [("Mass Density", rho, rhoans),
                                  ("Pressure", P, Pans),
                                  ("Velocity", v, vans),
                                  ("Thermal E", eps, epsans)]:
            assert len(data) == len(ans)
            error = [data[i] - ans[i] for i in xrange(len(data))]
            Pn = Pnorm.Pnorm(error, r)
            L1 = Pn.gridpnorm(1, rmin, rmax)
            L2 = Pn.gridpnorm(2, rmin, rmax)
            Linf = Pn.gridpnorm("inf", rmin, rmax)
            print "\t%s \t\t%g \t\t%g \t\t%g" % (name, L1, L2, Linf)
