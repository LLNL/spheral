#-------------------------------------------------------------------------------
# 3-D stirred box test
#
#
#-------------------------------------------------------------------------------
from math import *
from Spheral3d import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from findLastRestart import *
from SpheralVisitDump import dumpPhysicsState
import mpi

from GenerateNodeDistribution3d import *

title("3-D stirred box turbulence test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(seed = "lattice",

            nx = 100,
            ny = 100,
            nz = 100,

            rho1 = 1.0,
            eps1 = 0.0,
            nPerh = 2.01,

            gamma = 5.0/3.0,
            mu = 1.0,
            
            SVPH = False,
            CRKSPH = False,
            ASPH = False,
            SPH = True,   # This just chooses the H algorithm -- you can use this with CRKSPH for instance.
            filter = 0.0,   # CRKSPH filtering
            Qconstructor = MonaghanGingoldViscosity,
            #Qconstructor = TensorMonaghanGingoldViscosity,
            linearConsistent = False,
            fcentroidal = 0.0,
            fcellPressure = 0.0,
            boolReduceViscosity = False,
            nh = 5.0,
            aMin = 0.1,
            aMax = 2.0,
            Qhmult = 1.0,
            Cl = 1.0,
            Cq = 1.0,
            linearInExpansion = False,
            Qlimiter = False,
            balsaraCorrection = False,
            epsilon2 = 1e-2,
            hmin = 0.0001,
            hmax = 0.5,
            hminratio = 0.1,
            cfl = 0.5,
            useVelocityMagnitudeForDt = False,
            XSPH = False,
            epsilonTensile = 0.0,
            nTensile = 8,


            IntegratorConstructor = CheapSynchronousRK2Integrator,
            goalTime = 2.0,
            steps = None,
            vizCycle = None,
            vizTime = 0.1,
            dt = 0.0001,
            dtMin = 1.0e-8,
            dtMax = 0.1,
            dtGrowth = 2.0,
            maxSteps = None,
            statsStep = 10,
            smoothIters = 0,
            HUpdate = IdealH,
            domainIndependent = False,
            rigorousBoundaries = False,
            dtverbose = False,
            
            densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
            compatibleEnergy = True,            # <--- Important!  rigorousBoundaries does not work with the compatibleEnergy algorithm currently.
            gradhCorrection = False,
            
            useVoronoiOutput = False,
            clearDirectories = False,
            restoreCycle = None,
            restartStep = 100,
            redistributeStep = 500,
            checkRestart = False,
            dataDir = "stir-3d",
            outputFile = "None",
            comparisonFile = "None",
            
            serialDump = False, #whether to dump a serial ascii file at the end for viz
            )

# Decide on our hydro algorithm.
if SVPH:
  if ASPH:
    HydroConstructor = ASVPHFacetedHydro
  else:
    HydroConstructor = SVPHFacetedHydro
elif CRKSPH:
  Qconstructor = CRKSPHMonaghanGingoldViscosity
  if ASPH:
    HydroConstructor = ACRKSPHHydro
  else:
    HydroConstructor = CRKSPHHydro
else:
  if ASPH:
    HydroConstructor = ASPHHydro
  else:
    HydroConstructor = SPHHydro

xmin = (0.0, 0.0, 0.0)
xmax = (1.0, 1.0, 1.0)

dataDir = os.path.join(dataDir,
                       "rho1=%g" % rho1,
                       str(HydroConstructor).split("'")[1].split(".")[-1],
                       "densityUpdate=%s" % (densityUpdate),
                       "compatibleEnergy=%s" % (compatibleEnergy),
                       "XSPH=%s" % XSPH,
                       "filter=%s" % filter,
                       "%s-Cl=%g-Cq=%g" % (str(Qconstructor).split("'")[1].split(".")[-1], Cl, Cq),
                       "%ix%ix%i" % (nx, ny, nz),
                       "nPerh=%g-Qhmult=%g" % (nPerh, Qhmult))
restartDir = os.path.join(dataDir, "restarts")
vizDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, "stir-3d")
vizBaseName = "stir-3d"


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
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(BSplineKernel(), 1000)
WTPi = TableKernel(BSplineKernel(), 1000, Qhmult)
output("WT")
output("WTPi")
kernelExtent = WT.kernelExtent

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
nodes1 = makeFluidNodeList("High density gas", eos,
                           hmin = hmin,
                           hmax = hmax,
                           hminratio = hminratio,
                           nPerh = nPerh)
output("nodes1.nodesPerSmoothingScale")
output("nodes1.hmin")
output("nodes1.hmax")
output("nodes1.hminratio")

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    generator1 = GenerateNodeDistribution3d(nx, ny, nz, rho1, seed,
                                            xmin = xmin,
                                            xmax = xmax,
                                            nNodePerh = nPerh)

    if mpi.procs > 1:
      from VoronoiDistributeNodes import distributeNodes3d
    else:
      from DistributeNodes import distributeNodes3d
    distributeNodes3d((nodes1, generator1))
    output("mpi.reduce(nodes1.numInternalNodes, mpi.MIN)")
    output("mpi.reduce(nodes1.numInternalNodes, mpi.MAX)")
    output("mpi.reduce(nodes1.numInternalNodes, mpi.SUM)")

    # Set node specific thermal energies
    nodes1.specificThermalEnergy(ScalarField("tmp", nodes1, eps1))

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase3d()
output("db")
db.appendNodeList(nodes1)
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the artificial viscosity.
#-------------------------------------------------------------------------------
q = Qconstructor(Cl, Cq, linearInExpansion)
q.epsilon2 = epsilon2
q.limiter = Qlimiter
q.balsaraShearCorrection = balsaraCorrection
output("q")
output("q.Cl")
output("q.Cq")
output("q.epsilon2")
output("q.limiter")
output("q.balsaraShearCorrection")
output("q.linearInExpansion")
output("q.quadraticInExpansion")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if SVPH:
  hydro = HydroConstructor(WT, q,
                           cfl = cfl,
                           useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                           compatibleEnergyEvolution = compatibleEnergy,
                           densityUpdate = densityUpdate,
                           XSVPH = XSPH,
                           linearConsistent = linearConsistent,
                           generateVoid = False,
                           HUpdate = HUpdate,
                           fcentroidal = fcentroidal,
                           fcellPressure = fcellPressure,
                           xmin = Vector(-2.0, -2.0, -2.0),
                           xmax = Vector(3.0, 3.0, 3.0))
# xmin = Vector(x0 - 0.5*(x2 - x0), y0 - 0.5*(y2 - y0)),
# xmax = Vector(x2 + 0.5*(x2 - x0), y2 + 0.5*(y2 - y0)))
elif CRKSPH:
  hydro = HydroConstructor(WT, WTPi, q,
                           filter = filter,
                           cfl = cfl,
                           useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                           compatibleEnergyEvolution = compatibleEnergy,
                           XSPH = XSPH,
                           densityUpdate = densityUpdate,
                           HUpdate = HUpdate)
else:
  hydro = HydroConstructor(WT,
                           WTPi,
                           q,
                           cfl = cfl,
                           useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                           compatibleEnergyEvolution = compatibleEnergy,
                           gradhCorrection = gradhCorrection,
                           XSPH = XSPH,
                           densityUpdate = densityUpdate,
                           HUpdate = HUpdate,
                           epsTensile = epsilonTensile,
                           nTensile = nTensile)
output("hydro")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.densityUpdate")
output("hydro.HEvolution")

packages = [hydro]

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane0 = Plane(Vector(*xmin), Vector(1.0, 0.0, 0.0))
xPlane1 = Plane(Vector(*xmax), Vector(-1.0,0.0, 0.0))

yPlane0 = Plane(Vector(*xmin), Vector(0.0, 1.0, 0.0))
yPlane1 = Plane(Vector(*xmax), Vector(0.0, -1.0, 0.0))

zPlane0 = Plane(Vector(*xmin), Vector(0.0, 0.0, 1.0))
zPlane1 = Plane(Vector(*xmax), Vector(0.0, 0.0, -1.0))

xbc = PeriodicBoundary(xPlane0, xPlane1)
ybc = PeriodicBoundary(yPlane0, yPlane1)
zbc = PeriodicBoundary(zPlane0, zPlane1)

for p in packages:
    p.appendBoundary(xbc)
    p.appendBoundary(ybc)
    p.appendBoundary(zbc)

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
if useVoronoiOutput:
  import SpheralVoronoiSiloDump
  vizMethod = SpheralVoronoiSiloDump.dumpPhysicsState
else:
  import SpheralPointmeshSiloDump
  vizMethod = SpheralPointmeshSiloDump.dumpPhysicsState
control = SpheralController(integrator, WT,
                            initializeDerivatives = True,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            redistributeStep = redistributeStep,
                            vizMethod = vizMethod,
                            vizBaseName = vizBaseName,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime,
                            SPH = SPH)
output("control")

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
hstats([nodes1])
while control.time() < goalTime:
    nextGoalTime = min(control.time() + dtSample, goalTime)
    control.advance(nextGoalTime, maxSteps)
    control.dropRestartFile()
    dumpPhysicsState(integrator,
                     "Noh-planar-3d-%ix%ix%i" % (nx, ny, nz),
                     visitDir)
hstats([nodes1])

#-------------------------------------------------------------------------------
# Plot the elongation (h1/h2) for the H tensors.
#-------------------------------------------------------------------------------
import Gnuplot
hratio0 = [H.eigenValues().minElement()/H.eigenValues().maxElement() for H in nodes1.Hfield().internalValues()]
x0 = [ri.x for ri in nodes1.positions().internalValues()]
hratio = mpi.allreduce(hratio0, mpi.SUM)
x = mpi.allreduce(x0, mpi.SUM)
hratioPlot = Gnuplot.Gnuplot()
cache = []
if mpi.rank == 0:
    hratioData = Gnuplot.Data(x, hratio,
                              title = "hmin/hmax ratio",
                              inline = True)
    hratioPlot.plot(hratioData)
    cache.extend([hratioData, hratioPlot])

# Plot the final state.
vxPlot = plotFieldList(db.fluidVelocity, xFunction="%s.x", yFunction="%s.x", plotStyle="points", winTitle="x velocity")
vyPlot = plotFieldList(db.fluidVelocity, xFunction="%s.x", yFunction="%s.y", plotStyle="points", winTitle="y velocity")
vzPlot = plotFieldList(db.fluidVelocity, xFunction="%s.x", yFunction="%s.z", plotStyle="points", winTitle="z velocity")
rhoPlot = plotFieldList(db.fluidMassDensity, xFunction="%s.y", plotStyle="points", winTitle="Mass Density")
P = db.fluidPressure
PPlot = plotFieldList(P, xFunction="%s.x", plotStyle="points", winTitle="Pressure")
Hinverse = db.fluidHinverse
hx = db.newFluidScalarFieldList()
hy = db.newFluidScalarFieldList()
hz = db.newFluidScalarFieldList()
xunit = Vector3d(1, 0, 0)
yunit = Vector3d(0, 1, 0)
zunit = Vector3d(0, 0, 1)
for Hfield, hxfield, hyfield, hzfield in zip(Hinverse.fields(),
                                             hx.fields(),
                                             hy.fields(),
                                             hz.fields()):
    n = Hfield.numElements()
    assert hxfield.numElements() == n and hyfield.numElements() == n and hzfield.numElements() == n
    for i in xrange(n):
        hxfield[i] = (Hfield[i]*xunit).magnitude()
        hyfield[i] = (Hfield[i]*yunit).magnitude()
        hzfield[i] = (Hfield[i]*zunit).magnitude()
hxPlot = plotFieldList(hx, xFunction="%s.x", plotStyle="points", winTitle="h_x")
hyPlot = plotFieldList(hy, xFunction="%s.x", plotStyle="points", winTitle="h_y")
hzPlot = plotFieldList(hy, xFunction="%s.z", plotStyle="points", winTitle="h_z")

# Overplot the analytic solution.
from NohAnalyticSolution import *
answer = NohSolution(1,
                     h0 = nPerh*1.0/nx)
plotAnswer(answer, control.time(),
           rhoPlot = rhoPlot,
           velPlot = vyPlot,
           PPlot = PPlot,
           HPlot = hxPlot)
