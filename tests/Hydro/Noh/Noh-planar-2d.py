#-------------------------------------------------------------------------------
# The Planar Noh test case run in 1-D.
#
# W.F. Noh 1987, JCP, 72, 78-120.
#-------------------------------------------------------------------------------
from math import *
from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from findLastRestart import *
from SpheralVisitDump import dumpPhysicsState
import mpi

from GenerateNodeDistribution2d import *

title("2-D integrated hydro test -- Planar Noh problem")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(NodeListConstructor = AsphNodeList2d,

            seed = "lattice",

            nx = 20,
            ny = 100,

            rho1 = 1.0,
            eps1 = 0.0,
            vy = -1.0,
            vshear = 0.0,
            nPerh = 2.01,

            gamma = 5.0/3.0,
            mu = 1.0,
            Qconstructor = MonaghanGingoldViscosity2d,
            #Qconstructor = TensorMonaghanGingoldViscosity2d,
            Cl = 1.0,
            Cq = 1.5,
            Qlimiter = False,
            balsaraCorrection = False,
            epsilon2 = 1e-2,
            negligibleSoundSpeed = 1e-5,
            csMultiplier = 1e-4,
            energyMultiplier = 0.1,
            hmin = 1e-5,
            hmax = 1.0,
            hminratio = 0.05,
            HsmoothFraction = 0.0,
            cfl = 0.5,
            XSPH = True,
            epsilonTensile = 0.0,
            nTensile = 8,

            hourglassMultiplier = 0.0,
            hourglassAccelerationFactor = 0.01,

            useIdealHUpdate = False,
            idealHnPerh = 2.01,
            snapRound = True,

            HEvolution = Hydro2d.HEvolutionType.IdealH,
            limitIdealH = False,

            neighborSearchType = Neighbor2d.NeighborSearchType.GatherScatter,
            numGridLevels = 20,
            topGridCellSize = 2.0,
            origin = Vector2d(0.0, 0.0),

            goalTime = 0.6,
            dt = 0.0001,
            dtMin = 1.0e-5, 
            dtMax = None,
            dtGrowth = 2.0,
            dtSample = 0.1,
            rigorousBoundaries = False,
            maxSteps = None,
            statsStep = 1,
            smoothIters = 0,
            sumForMassDensity = Hydro2d.MassDensityType.RigorousSumDensity, # VolumeScaledDensity,
            compatibleEnergy = True,

            restoreCycle = None,
            restartStep = 1000,

            # Parameters for the test acceptance.,
            L1rho = 0.702502,
            L2rho = 1.01547,
            Linfrho = 3.07796,
            
            L1P = 0.194579,
            L2P = 0.341756,
            LinfP = 1.3412,
            
            L1v = 0.339706,
            L2v = 0.478579,
            Linfv = 1.10007,
            
            L1eps = 0.111857,
            L2eps = 0.170687,
            Linfeps = 0.493335,

            L1h = 0.00698501,
            L2h = 0.00828045,
            Linfh = 0.0149375,

            tol = 1.0e-5,

            graphics = "gnu",
            )

xmin = (0.0, 0.0)
xmax = (1.0, 1.0)
#xmax = (0.2, 1.0)

dataDir = "dumps-planar-%ix%i" % (nx, ny)
restartDir = dataDir + "/restarts"
visitDir = dataDir + "/visit"
restartBaseName = restartDir + "/Noh-planar-2d-%ix%i" % (nx, ny)

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
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS2d(gamma, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel2d(BSplineKernel2d(), 100)
WTPi = TableKernel2d(BSplineKernel2d(), 100)
output("WT")
output("WTPi")
kernelExtent = WT.kernelExtent()

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
nodes1 = NodeListConstructor("nodes", eos, WT, WTPi)
output("nodes1")
nodes1.HsmoothFraction = HsmoothFraction
nodes1.XSPH = XSPH
nodes1.nodesPerSmoothingScale = nPerh
nodes1.epsilonTensile = epsilonTensile
nodes1.nTensile = nTensile
nodes1.hmin = hmin
nodes1.hmax = hmax
nodes1.hminratio = hminratio
output("nodes1.HsmoothFraction")
output("nodes1.nodesPerSmoothingScale")
output("nodes1.epsilonTensile")
output("nodes1.nTensile")
output("nodes1.XSPH")
output("nodes1.hmin")
output("nodes1.hmax")
output("nodes1.hminratio")

#-------------------------------------------------------------------------------
# Construct the neighbor object.
#-------------------------------------------------------------------------------
neighbor1 = NestedGridNeighbor2d(nodes1,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes1.registerNeighbor(neighbor1)

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    #from ParMETISDistributeNodes import distributeNodes2d
    from PeanoHilbertDistributeNodes import distributeNodes2d
    #from MortonOrderDistributeNodes import distributeNodes2d
    generator1 = GenerateNodeDistribution2d(nx, ny, rho1, seed,
                                            xmin = xmin,
                                            xmax = xmax,
                                            nNodePerh = nPerh)
    distributeNodes2d((nodes1, generator1))
    output("mpi.reduce(nodes1.numInternalNodes, mpi.MIN)")
    output("mpi.reduce(nodes1.numInternalNodes, mpi.MAX)")
    output("mpi.reduce(nodes1.numInternalNodes, mpi.SUM)")

    # Set node specific thermal energies
    nodes1.specificThermalEnergy(ScalarField2d("tmp", nodes1, eps1))

    # Set node velocities
    for i in xrange(nodes1.numInternalNodes):
        y = nodes1.positions()[i].y
        nodes1.velocity()[i] = Vector2d(vshear*cos(2.0*pi*y), vy)

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase2d()
output("db")
output("db.appendNodeList(nodes1)")
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
hydro = Hydro2d(WT, WTPi, q, compatibleEnergy)
hydro.cfl = cfl
hydro.HEvolution = HEvolution
hydro.limitIdealH = limitIdealH
hydro.sumForMassDensity = sumForMassDensity
hydro.HsmoothMin = hmin
hydro.HsmoothMax = hmax
output("hydro")
output("hydro.cfl")
output("hydro.HEvolution")
output("hydro.sumForMassDensity")
output("hydro.HsmoothMin")
output("hydro.HsmoothMax")
output("hydro.compatibleEnergyEvolution")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.valid()")

packages = [hydro]

#-------------------------------------------------------------------------------
# Optionally construct an hour glass control object.
#-------------------------------------------------------------------------------
if hourglassMultiplier > 0.0:
    hourglass = SecondMomentHourglassControl2d(hourglassMultiplier,
                                               hourglassAccelerationFactor)
    output("hourglass")
    output("hourglass.multiplier")
    output("hourglass.maxAccelerationFactor")

    packages += [hourglass]

#-------------------------------------------------------------------------------
# Optionally construct an ideal H update object.
#-------------------------------------------------------------------------------
if useIdealHUpdate:
    idealHUpdate = IdealHUpdate2d(idealHnPerh, snapRound)
    output("idealHUpdate")
    output("idealHUpdate.updateNodesPerSmoothingScale")
    output("idealHUpdate.snapRound")
    packages += [idealHUpdate]

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane0 = Plane2d(Vector2d(*xmin), Vector2d(1.0, 0.0))
xPlane1 = Plane2d(Vector2d(*xmax), Vector2d(-1.0, 0.0))
yPlane0 = Plane2d(Vector2d(*xmin), Vector2d(0.0, 1.0))
xbc = PeriodicBoundary2d(xPlane0, xPlane1)
ybc = ReflectingBoundary2d(yPlane0)
for p in packages:
    p.appendBoundary(xbc)
    p.appendBoundary(ybc)

#-------------------------------------------------------------------------------
# Construct a predictor corrector integrator.
#-------------------------------------------------------------------------------
integrator = SynchronousRK2Integrator2d(db)
for p in packages:
    integrator.appendPhysicsPackage(p)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.rigorousBoundaries = rigorousBoundaries
output("integrator")
output("integrator.valid()")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.rigorousBoundaries")

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName)
output("control")

# Smooth the initial conditions.
if restoreCycle is not None:
    control.loadRestartFile(restoreCycle)
else:
    control.iterateIdealH()
    control.smoothState(smoothIters)
    dumpPhysicsState(integrator,
                     "Noh-planar-2d-%ix%i" % (nx, ny),
                     visitDir)

control.step(2)

## #-------------------------------------------------------------------------------
## # Advance to the end time.
## #-------------------------------------------------------------------------------
## while control.time() < goalTime:
##     nextGoalTime = min(control.time() + dtSample, goalTime)
##     control.advance(nextGoalTime, maxSteps)
##     control.dropRestartFile()
##     dumpPhysicsState(integrator,
##                      "Noh-planar-2d-%ix%i" % (nx, ny),
##                      visitDir)

## #-------------------------------------------------------------------------------
## # Plot the elongation (h1/h2) for the H tensors.
## #-------------------------------------------------------------------------------
## import Gnuplot
## hratio0 = [H.eigenValues().minElement()/H.eigenValues().maxElement() for H in nodes1.Hfield().internalValues()]
## y0 = [ri.y for ri in nodes1.positions().internalValues()]
## hratio = mpi.allreduce(hratio0, mpi.SUM)
## y = mpi.allreduce(y0, mpi.SUM)
## hratioPlot = Gnuplot.Gnuplot()
## cache = []
## if mpi.rank == 0:
##     hratioData = Gnuplot.Data(y, hratio,
##                               title = "hmin/hmax ratio",
##                               inline = True)
##     hratioPlot.plot(hratioData)
##     cache.extend([hratioData, hratioPlot])

## rPlot = plotNodePositions2d(db, colorNodeLists=0, colorDomains=1)

## # Plot the final state.
## vxPlot = plotFieldList(db.fluidVelocity, xFunction="%s.y", yFunction="%s.x", plotStyle="points", winTitle="x velocity")
## vyPlot = plotFieldList(db.fluidVelocity, xFunction="%s.y", yFunction="%s.y", plotStyle="points", winTitle="y velocity")
## rhoPlot = plotFieldList(db.fluidMassDensity, xFunction="%s.y", plotStyle="points", winTitle="Mass Density")
## P = db.fluidPressure
## PPlot = plotFieldList(P, xFunction="%s.y", plotStyle="points", winTitle="Pressure")
## Hinverse = db.fluidHinverse
## hx = db.newFluidScalarFieldList()
## hy = db.newFluidScalarFieldList()
## xunit = Vector2d(1, 0)
## yunit = Vector2d(0, 1)
## for Hfield, hxfield, hyfield in zip(Hinverse.fields(),
##                                     hx.fields(),
##                                     hy.fields()):
##     n = Hfield.numElements()
##     assert hxfield.numElements() == n and hyfield.numElements() == n
##     for i in xrange(n):
##         hxfield[i] = (Hfield[i]*xunit).magnitude()
##         hyfield[i] = (Hfield[i]*yunit).magnitude()
## hxPlot = plotFieldList(hx, xFunction="%s.y", plotStyle="points", winTitle="h_x")
## hyPlot = plotFieldList(hy, xFunction="%s.y", plotStyle="points", winTitle="h_y")

## # Overplot the analytic solution.
## from NohAnalyticSolution import *
## answer = NohSolution(1,
##                      h0 = nPerh*1.0/ny)
## plotAnswer(answer, control.time(),
##            rhoPlot = rhoPlot,
##            velPlot = vyPlot,
##            PPlot = PPlot,
##            HPlot = hyPlot)


## # Connect the dots!
## xx = [x.x for x in nodes1.positions().internalValues()]
## yy = [x.y for x in nodes1.positions().internalValues()]
## d = Gnuplot.Data(xx, yy,
##                  with = "lines")
## p = Gnuplot.Gnuplot()
## p.plot(d)
## p("set xrange [0:1]; set yrange [0:1]"); p.refresh()

