#-------------------------------------------------------------------------------
# A mock ICF kind 'o problem.
#-------------------------------------------------------------------------------
from math import *
from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from SpheralVisitDump import dumpPhysicsState
from findLastRestart import findLastRestart
from GzipFileNodeGenerator import *

# Load the mpi module if we're parallel.
import loadmpi
mpi, rank, procs = loadmpi.loadmpi()

title("2-D ICF test problem")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(NodeListConstructor = AsphNodeList2d,

            rhoAir = 1.4,
            rhoDrive = 0.1,
            rhoShell = 1.0,
            PAir = 1.0,
            PDrive = 100.0,
            PShell = 1.0,
            gammaAir = 1.4,
            gammaDrive = 1.4,
            gammaShell = 1.6,
            mu = 1.0,

            Qconstructor = MonaghanGingoldViscosity2d,
            #Qconstructor = TensorMonaghanGingoldViscosity2d,
            Cl = 1.0, 
            Cq = 0.75,
            Qlimiter = False,
            balsaraCorrection = False,
            epsilon2 = 1e-2,
            hmin = 1e-5,
            hmax = 0.5,
            hminratio = 0.1,
            nPerh = 2.01,
            cfl = 0.5,
            XSPH = True,
            epsilonTensile = 0.0,
            nTensile = 8,

            goalTime = 0.06,
            dtSample = 0.01,
            dt = 0.0001,
            dtMin = 1.0e-6, 
            dtMax = 0.1,
            dtGrowth = 2.0,
            maxSteps = None,
            statsStep = 10,
            smoothIters = 0,
            HEvolution = Hydro2d.HEvolutionType.IdealH,
            sumForMassDensity = Hydro2d.MassDensityType.RigorousSumDensity,
            compatibleEnergy = True,

            restoreCycle = None,
            restartStep = 1000,
            dataDirBase = "icf-2d",

            graphics = True,
            )

epsAir = PAir/((gammaAir - 1.0)*rhoAir)
epsDrive = PDrive/((gammaDrive - 1.0)*rhoDrive)
epsShell = PShell/((gammaShell - 1.0)*rhoShell)

dataDir = dataDirBase
restartDir = dataDir + "/restarts"
visitDir = dataDir + "/visit"
restartBaseName = restartDir + "/icf-2d"

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
eosAir = GammaLawGasCGS2d(gammaAir, mu)
eosDrive = GammaLawGasCGS2d(gammaDrive, mu)
eosShell = GammaLawGasCGS2d(gammaShell, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel2d(BSplineKernel2d(), 1000)
WTPi = TableKernel2d(BSplineKernel2d(), 1000)
output("WT")
output("WTPi")
kernelExtent = WT.kernelExtent()

#-------------------------------------------------------------------------------
# Make the NodeLists.
#-------------------------------------------------------------------------------
nodesAir = NodeListConstructor("air", eosAir, WT, WTPi)
nodesDrive = NodeListConstructor("drive", eosAir, WT, WTPi)
nodesShell = NodeListConstructor("shell", eosAir, WT, WTPi)
nodeSet = [nodesAir, nodesDrive, nodesShell]
for nodes in nodeSet:
    nodes.XSPH = XSPH
    nodes.hmin = hmin
    nodes.hmax = hmax
    nodes.hminratio = hminratio
    nodes.nodesPerSmoothingScale = nPerh
    nodes.epsilonTensile = epsilonTensile
    nodes.nTensile = nTensile
    output("nodes.name()")
    output("    nodes.hmin")
    output("    nodes.hmax")
    output("    nodes.hminratio")
    output("    nodes.nodesPerSmoothingScale")
    output("    nodes.epsilonTensile")
    output("    nodes.nTensile")
    output("    nodes.XSPH")

#-------------------------------------------------------------------------------
# Construct the neighbor objects.
#-------------------------------------------------------------------------------
_cache = []
for nodes in nodeSet:
    neighbor = TreeNeighbor2d(nodes,
                              kernelExtent = kernelExtent)
    nodes.registerNeighbor(neighbor)
    _cache.append(neighbor)

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
outerNodeFlags = IntField2d("outer node flags", nodesDrive)
if restoreCycle is None:
    filename = "icf-10-20-8x90.gz"
    generatorAir = GzipFileNodeGeneratorRZto2D(filename, "Inner", nPerh,
                                               SPH = (type(nodesAir) == SphNodeList2d))
    generatorDrive = GzipFileNodeGeneratorRZto2D(filename, "Driver", nPerh,
                                                 SPH = (type(nodesDrive) == SphNodeList2d),
                                                 extraFields = ["Driver"])
    generatorShell = GzipFileNodeGeneratorRZto2D(filename, "Shell", nPerh,
                                                 SPH = (type(nodesShell) == SphNodeList2d))

    # Get the outer node flags.
    nodesDrive.numInternalNodes = generatorDrive.numLocalNodes()
    for i in xrange(generatorDrive.numLocalNodes()):
        outerNodeFlags[i] = int(generatorDrive.outerNodes[i] + 0.1)

    from ParMETISDistributeNodes import distributeNodes2d
    distributeNodes2d((nodesAir, generatorAir),
                      (nodesDrive, generatorDrive),
                      (nodesShell, generatorShell))
    for nodes in nodeSet:
        output("nodes.name()")
        output("    mpi.reduce(nodes.numInternalNodes, mpi.MIN)")
        output("    mpi.reduce(nodes.numInternalNodes, mpi.MAX)")
        output("    mpi.reduce(nodes.numInternalNodes, mpi.SUM)")

    # Set node specific thermal energies
    nodesAir.specificThermalEnergy(ScalarField2d("tmp", nodesAir, epsAir))
    nodesDrive.specificThermalEnergy(ScalarField2d("tmp", nodesDrive, epsAir))
    nodesShell.specificThermalEnergy(ScalarField2d("tmp", nodesShell, epsAir))

#-------------------------------------------------------------------------------
# Construct a DataBase.
#-------------------------------------------------------------------------------
db = DataBase2d()
for nodes in nodeSet:
    db.appendNodeList(nodes)
output("db")
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
hydro = Hydro2d(WT, WTPi, q, compatibleEnergy)
hydro.cfl = cfl
hydro.HEvolution = HEvolution
hydro.sumForMassDensity = sumForMassDensity
hydro.HsmoothMin = hmin
hydro.HsmoothMax = hmax
hydro.HratioMin = hminratio
output("hydro")
output("hydro.cfl")
output("hydro.HEvolution")
output("hydro.sumForMassDensity")
output("hydro.HsmoothMin")
output("hydro.HsmoothMax")
output("hydro.compatibleEnergyEvolution")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.HratioMin")
output("hydro.valid()")

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane0 = Plane2d(Vector2d(0.0, 0.0), Vector2d(1.0, 0.0))
yPlane0 = Plane2d(Vector2d(0.0, 0.0), Vector2d(0.0, 1.0))
xbc0 = ReflectingBoundary2d(xPlane0)
ybc0 = ReflectingBoundary2d(yPlane0)
hydro.appendBoundary(xbc0)
hydro.appendBoundary(ybc0)
output("hydro.haveBoundary(xbc0)")
output("hydro.haveBoundary(ybc0)")

#-------------------------------------------------------------------------------
# Construct an integrator, and add the physics packages.
#-------------------------------------------------------------------------------
integrator = SynchronousRK2Integrator2d(db)
integrator.appendPhysicsPackage(hydro)
integrator.lastDt = dt
integrator.dtMin = dtMin
integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.valid()")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")

#-------------------------------------------------------------------------------
# Make the problem controller.
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
    control.dropRestartFile()
    dumpPhysicsState(integrator,
                     "icf-2d",
                     visitDir)

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
hstats(nodeSet)
while control.time() < goalTime:
    nextGoalTime = min(control.time() + dtSample, goalTime)
    control.advance(nextGoalTime, maxSteps)
    control.dropRestartFile()
    dumpPhysicsState(integrator,
                     "icf-2d",
                     visitDir)

#-------------------------------------------------------------------------------
# Plot the results.
#-------------------------------------------------------------------------------
if graphics:

    # Plot the elongation (h1/h2) for the H tensors.
    import Gnuplot
    rPlot = plotNodePositions2d(db, colorNodeLists=True, colorDomains=False)

    # Plot the final state.
    rhoPlot, vrPlot, epsPlot, PPlot, HPlot = plotRadialState(db)
    del HPlot
    Hinverse = db.fluidHinverse
    hr = db.newFluidScalarFieldList()
    ht = db.newFluidScalarFieldList()
    for Hfield, hrfield, htfield in zip(Hinverse.fields(),
                                        hr.fields(),
                                        ht.fields()):
        n = Hfield.numElements()
        assert hrfield.numElements() == n
        assert htfield.numElements() == n
        positions = Hfield.nodeList().positions()
        for i in xrange(n):
            runit = positions[i].unitVector()
            tunit = Vector2d(-(positions[i].y), positions[i].x).unitVector()
            hrfield[i] = (Hfield[i]*runit).magnitude()
            htfield[i] = (Hfield[i]*tunit).magnitude()
    hrPlot = plotFieldList(hr, xFunction="%s.magnitude()", plotStyle="points", winTitle="h_r")
    htPlot = plotFieldList(ht, xFunction="%s.magnitude()", plotStyle="points", winTitle="h_t")
