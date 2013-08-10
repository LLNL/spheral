from Spheral1d import *
from SpheralTestUtilities import *
from SodAnalyticSolution import *

title("1-D integrated hydro test -- planar Sod problem")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(nx1 = 400,
            nx2 = 100,
            rho1 = 1.0,
            rho2 = 0.25,
            P1 = 1.0,
            P2 = 0.1795,

            numNodeLists = 1,

            x0 = -0.5,
            x1 = 0.0,
            x2 = 0.5,

            smoothDiscontinuity = False,

            nPerh = 2.01,

            gammaGas = 5.0/3.0,
            mu = 1.0,
            Qconstructor = MonaghanGingoldViscosity,
            Cl = 1.0,
            Cq = 1.5,
            Qlimiter = False,
            epsilon2 = 1e-4,
            hmin = 1e-10,
            hmax = 1.0,
            cfl = 0.5,
            XSPH = False,
            epsilonTensile = 0.0,
            nTensile = 8,
            rhoMin = 0.01,
            hourglass = None,
            hourglassOrder = 1,
            hourglassLimiter = 1,

            SVPH = False,
            IntegratorConstructor = CheapSynchronousRK2Integrator,
            steps = None,
            goalTime = 0.15,
            dt = 1e-4,
            dtMin = 1.0e-5,
            dtMax = 0.1,
            dtGrowth = 2.0,
            rigorousBoundaries = False,
            maxSteps = None,
            statsStep = 10,
            HEvolution = IdealH,
            densityUpdate = RigorousSumDensity,
            compatibleEnergy = True,
            gradhCorrection = True,

            useRefinement = False,

            restoreCycle = None,
            restartStep = 200,
            dataDirBase = "Sod-planar-1d",
            restartBaseName = "Sod-planar-1d-restart",

            graphics = "gnu",
            )

dataDir = dataDirBase + ("/%i" % (nx1 + nx2))
restartDir = dataDir + "/restarts"
restartBaseName = restartDir + "/Sod-planar-1d-%i" % (nx1 + nx2)

assert numNodeLists in (1, 2)

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS(gammaGas, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(BSplineKernel(), 1000)
WTPi = WT
output("WT")
output("WTPi")

#-------------------------------------------------------------------------------
# Make the NodeLists.
#-------------------------------------------------------------------------------
nodes1 = makeFluidNodeList("nodes1", eos, 
                           hmin = hmin,
                           hmax = hmax,
                           nPerh = nPerh,
                           rhoMin = rhoMin)
nodes2 = makeFluidNodeList("nodes2", eos, 
                           hmin = hmin,
                           hmax = hmax,
                           nPerh = nPerh,
                           rhoMin = rhoMin)
nodeSet = [nodes1, nodes2]

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
from DistributeNodes import distributeNodesInRange1d
if numNodeLists == 1:
    distributeNodesInRange1d([(nodes1, [(nx1, rho1, (x0, x1)), (nx2, rho2, (x1, x2))])])
else:
    distributeNodesInRange1d([(nodes1, [(nx1, rho1, (x0, x1))]),
                              (nodes2, [(nx2, rho2, (x1, x2))])])
output("nodes1.numNodes")
output("nodes2.numNodes")

# Set node specific thermal energies
def specificEnergy(x):
    if x < x1:
        return P1/((gammaGas - 1.0)*rho1)
    else:
        return P2/((gammaGas - 1.0)*rho2)
for nodes in nodeSet:
    pos = nodes.positions()
    eps = nodes.specificThermalEnergy()
    for i in xrange(nodes.numInternalNodes):
        eps[i] = specificEnergy(pos[i].x)

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
for nodes in nodeSet:
    output("db.appendNodeList(nodes)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the artificial viscosity.
#-------------------------------------------------------------------------------
q = Qconstructor(Cl, Cq)
q.limiter = Qlimiter
q.epsilon2 = epsilon2
output("q")
output("q.Cl")
output("q.Cq")
output("q.limiter")
output("q.epsilon2")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if SVPH:
    hydro = SVPHHydro(WT, q,
                      cfl = cfl,
                      compatibleEnergyEvolution = compatibleEnergy,
                      XSVPH = XSPH,
                      densityUpdate = densityUpdate,
                      HUpdate = HEvolution,
                      xmin = Vector(-100.0),
                      xmax = Vector( 100.0))
else:
    hydro = SPHHydro(WT,
                     WTPi,
                     q,
                     cfl = cfl,
                     compatibleEnergyEvolution = compatibleEnergy,
                     gradhCorrection = gradhCorrection,
                     densityUpdate = densityUpdate,
                     HUpdate = HEvolution,
                     XSPH = XSPH,
                     epsTensile = epsilonTensile,
                     nTensile = nTensile)
output("hydro")

packages = [hydro]

#-------------------------------------------------------------------------------
# Optionally construct an hourglass control object.
#-------------------------------------------------------------------------------
if hourglass:
    hg = hourglass(WT, hourglassOrder, hourglassLimiter)
    output("hg")
    output("hg.order")
    output("hg.limiter")
    packages.append(hg)

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane0 = Plane(Vector(x0), Vector( 1.0))
xPlane1 = Plane(Vector(x2), Vector(-1.0))
xbc0 = ReflectingBoundary(xPlane0)
xbc1 = ReflectingBoundary(xPlane1)

for p in packages:
    p.appendBoundary(xbc0)
    p.appendBoundary(xbc1)

#-------------------------------------------------------------------------------
# Construct an integrator, and add the one physics package.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
for p in packages:
    integrator.appendPhysicsPackage(p)
integrator.lastDt = dt
integrator.dtGrowth = dtGrowth
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.rigorousBoundaries = rigorousBoundaries
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
if hourglass:
    output("integrator.havePhysicsPackage(hg)")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.rigorousBoundaries")

#-------------------------------------------------------------------------------
# If requested, smooth the state at the discontinuity.
#-------------------------------------------------------------------------------
if smoothDiscontinuity:
    W0 = WT.kernelValue(0.0, 1.0)
    m1 = (x1 - x0)*rho1/nx1
    m2 = (x2 - x1)*rho2/nx2
    dx1 = (x1 - x0)/nx1
    dx2 = (x2 - x1)/nx2
    h1 = 2.0*dx1 * nPerh
    h2 = 2.0*dx2 * nPerh
    H1 = SymTensor(1.0/(nPerh * dx1))
    H2 = SymTensor(1.0/(nPerh * dx2))
    A1 = P1/(rho1**gammaGas)
    A2 = P2/(rho2**gammaGas)

    for i in xrange(nodes1.numInternalNodes):
        fi = WT.kernelValue(abs(nodes1.positions()[i].x - x1)/h1, 1.0) / W0
        fi = 1.0 - min(1.0, abs(nodes1.positions()[i].x - x1)/h1)
        assert fi >= 0.0 and fi <= 1.0
        nodes1.mass()[i] = (1.0 - fi)*m1 + 0.5*fi*(m1 + m2)
    for i in xrange(nodes2.numInternalNodes):
        fi = WT.kernelValue(abs(nodes2.positions()[i].x - x1)/h2, 1.0) / W0
        fi = 1.0 - min(1.0, abs(nodes2.positions()[i].x - x1)/h2)
        assert fi >= 0.0 and fi <= 1.0
        nodes2.mass()[i] = (1.0 - fi)*m2 + 0.5*fi*(m1 + m2)

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle)
output("control")

#-------------------------------------------------------------------------------
# If we want to use refinement, build the refinemnt algorithm.
#-------------------------------------------------------------------------------
class SelectNodes:
    def __init__(self):
        return
    def prepareForSelection(self, dataBase):
        return
    def selectNodes(self, nodeList):
        if control.totalSteps == 50 and nodeList.name == nodes1.name:
            return range(nodes1.numInternalNodes)
        else:
            return []

class SelectNodes2:
    def __init__(self):
        return
    def prepareForSelection(self, dataBase):
        return
    def selectNodes(self, nodeList):
        if control.totalSteps == 20:
            return range(nodeList.numInternalNodes)
        else:
            return []

if useRefinement:
    from AdaptiveRefinement import *
    select = SelectNodes()
    refine = SplitNodes1d()
    package = AdaptiveRefinement(db, select, refine, control)
    control.appendPeriodicWork(package.refineNodes, 1)

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)
else:
    control.advance(goalTime, maxSteps)
    control.dropRestartFile()

#-------------------------------------------------------------------------------
# Plot the final state.
#-------------------------------------------------------------------------------
if graphics in ("gnu", "matplot"):
    if graphics == "gnu":
        from SpheralGnuPlotUtilities import *
    elif graphics == "matplot":
        from SpheralMatplotlibUtilities import *

    rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db)

    # Now overplot the analytic solution.
    dx1 = (x1 - x0)/nx1
    dx2 = (x2 - x1)/nx2
    h1 = 1.0/(nPerh*dx1)
    h2 = 1.0/(nPerh*dx2)
    answer = SodSolution(nPoints=nx1 + nx2,
                         gamma = gammaGas,
                         rho1 = rho1,
                         P1 = P1,
                         rho2 = rho2,
                         P2 = P2,
                         x0 = x0,
                         x1 = x1,
                         x2 = x2,
                         h1 = 1.0/h1,
                         h2 = 1.0/h2)

    plotAnswer(answer, control.time(),
               rhoPlot, velPlot, epsPlot, PPlot, HPlot)
    pE = plotEHistory(control.conserve)

    def createList(x):
        xx = x
        if xx == []:
            xx = [-1e50,]
        tmp = mpi.allreduce(xx, mpi.SUM)
        result = [y for y in tmp if y != -1e50]
        return result

    # Compute the simulated specific entropy.
    A = []
    for nodes in nodeSet:
        rho = createList(nodes.massDensity().internalValues())
        pressure = ScalarField("pressure", nodes)
        nodes.pressure(pressure)
        P = createList(pressure.internalValues())
        A += [Pi/rhoi**gammaGas for (Pi, rhoi) in zip(P, rho)]

    # The analytic solution for the simulated entropy.
    xprof = createList([x.x for x in nodes1.positions().internalValues()] +
                       [x.x for x in nodes2.positions().internalValues()])
    xans, vans, uans, rhoans, Pans, hans = answer.solution(control.time(), xprof)
    Aans = [Pi/rhoi**gammaGas for (Pi, rhoi) in zip(Pans,  rhoans)]

    # # Plot the specific entropy.
    # if mpi.rank == 0:
    #     ll = zip(xprof, A, Aans)
    #     ll.sort()
    #     AsimData = Gnuplot.Data(xprof, A,
    #                             with_ = "points",
    #                             title = "Simulation",
    #                             inline = True)
    #     AansData = Gnuplot.Data(xprof, Aans,
    #                             with_ = "lines",
    #                             title = "Analytic",
    #                             inline = True)
    #     Aplot = Gnuplot.Gnuplot()
    #     Aplot.plot(AsimData)
    #     Aplot.replot(AansData)
    # else:
    #     Aplot = fakeGnuplot()

    # # Plot the grad h correction term (omega)
    # omegaPlot = plotFieldList(db.fluidOmegaGradh,
    #                           winTitle = "grad h correction",
    #                           colorNodeLists = False)

print "Energy conservation: original=%g, final=%g, error=%g" % (control.conserve.EHistory[0],
                                                                control.conserve.EHistory[-1],
                                                                (control.conserve.EHistory[-1] - control.conserve.EHistory[0])/control.conserve.EHistory[0])
