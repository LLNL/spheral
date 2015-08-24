#-------------------------------------------------------------------------------
# 2-D test of the Kelvin-Helmholtz instability.
# Designed to match the setup described in
# Agertz et al. 2006 (arXiv:astro-ph/0610051)
#-------------------------------------------------------------------------------
from math import *
from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from SpheralVisitDump import dumpPhysicsState
from findLastRestart import *
from GenerateNodeDistribution2d import *

import loadmpi
mpi, rank, procs = loadmpi.loadmpi()

title("2-D Kelvin-Helmholtz shearing test.")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(NodeListConstructor = SphNodeList2d,
            seed = "lattice",
            nx1 = 100,
            nx2 = 100,
            ny1 = 50,
            ny2 = 50,
            nPerh = 2.01,

            gamma = 5.0/3.0,
            mu = 1.0,
            rho1 = 1.0,
            rho2 = 2.0,
            P1 = 2.5,
            v1 = -0.5,
            v2 = 0.5,
            perturbWaves = 2.0,
            perturbHeight = 0.05,
            perturbMag = 0.025,

            Qconstructor = TensorMonaghanGingoldViscosity2d, #MonaghanGingoldViscosity2d
            Cl = 1.0,
            Cq = 0.75,
            epsilon2 = 1e-4,
            Qlimiter = True,
            balsaraCorrection = False,
            hmin = 1e-5,
            hmax = 0.2,
            hminratio = 0.1,
            cfl = 0.5,
            sumForMassDensity = Hydro2d.MassDensityType.RigorousSumDensity, #IntegrateDensity, # 
            HEvolution = Hydro2d.HEvolutionType.IdealH,
            XSPH = False,
            compatibleEnergy = False,
            rigorousBoundaries = True,

            IntegratorConstructor = SynchronousRK2Integrator2d,

            neighborSearchType = Neighbor2d.NeighborSearchType.GatherScatter,
            numGridLevels = 25,
            topGridCellSize = 1.0,

            goalTime = 5.0,
            nSample = 20,
            dtMin = 1.0e-5,
            dtMax = 0.1,
            dtGrowth = 2.0,
            maxSteps = None,
            statsStep = 10,
            smoothIters = 0,
            restartStep = 20,
            dataDirBase = "dumps-KelvinHelmholtz-2d",

            graphics = True,
            )

xmin1, xmax1 = Vector2d(0.0, -0.5), Vector2d(1.0, 0.0)
xmin2, xmax2 = Vector2d(0.0,  0.0), Vector2d(1.0, 0.5)
origin = Vector2d(0.0, 0.0)

dataDir = "%s/%s/Cl=%3.1f_Cq=%3.1f/nPerh=%4.2f/KelvinHelmholtz-2d-%ix%i" % (dataDirBase,
                                                                            str(Qconstructor).split("'")[1],
                                                                            Cl, Cq,
                                                                            nPerh,
                                                                            nx1, ny1 + ny2)
restartDir = dataDir + "/restarts"
visitDir = dataDir + "/visit"
restartBaseName = restartDir + "/KelvinHelmholtz-2d"

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
#eos = PolytropicEquationOfStateMKS2d(1.0, 1.0, mu)

eps1 = P1/((gamma - 1.0)*rho1)
eps2 = rho1/rho2*eps1
P2 = eos.pressure(rho2, eps2)
c1 = eos.soundSpeed(rho1, eps1)
c2 = eos.soundSpeed(rho2, eps2)
assert fuzzyEqual(P1, P2)

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
nodes1 = NodeListConstructor("Bottom nodes", eos, WT, WTPi)
nodes2 = NodeListConstructor("Top nodes", eos, WT, WTPi)
nodeSet = [nodes1, nodes2]
for nodes in nodeSet:
    nodes.XSPH = XSPH
    nodes.hmin = hmin
    nodes.hmax = hmax
    nodes.hminratio = hminratio
    nodes.nodesPerSmoothingScale = nPerh
    output("nodes")
    output("nodes.XSPH")
    output("nodes.hmin")
    output("nodes.hmax")
    output("nodes.hminratio")
    output("nodes.nodesPerSmoothingScale")
del nodes

#-------------------------------------------------------------------------------
# Construct the neighbor objects.
#-------------------------------------------------------------------------------
neighbor1 = NestedGridNeighbor2d(nodes1,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes1.registerNeighbor(neighbor1)
neighbor2 = NestedGridNeighbor2d(nodes2,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes2.registerNeighbor(neighbor2)

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    import ParMETISDistributeNodes
    generator1 = GenerateNodeDistribution2d(nx1, ny1, rho1, seed,
                                            xmin = (xmin1.x, xmin1.y),
                                            xmax = (xmax1.x, xmax1.y),
                                            nNodePerh = nPerh,
                                            SPH = (NodeListConstructor == SphNodeList2d))
    generator2 = GenerateNodeDistribution2d(nx2, ny2, rho2, seed,
                                            xmin = (xmin2.x, xmin2.y),
                                            xmax = (xmax2.x, xmax2.y),
                                            nNodePerh = nPerh,
                                            SPH = (NodeListConstructor == SphNodeList2d))
    ParMETISDistributeNodes.distributeNodes2d((nodes1, generator1),
                                              (nodes2, generator2))
    for nodes in nodeSet:
        output("nodes.name()")
        output("mpi.reduce(nodes.numInternalNodes, mpi.MIN)")
        output("mpi.reduce(nodes.numInternalNodes, mpi.MAX)")
        output("mpi.reduce(nodes.numInternalNodes, mpi.SUM)")
    del nodes

    # Set the node densities.
    nodes1.massDensity(ScalarField2d("tmp", nodes1, rho1))
    nodes2.massDensity(ScalarField2d("tmp", nodes2, rho2))

    # Set the node specific energies.
    nodes1.specificThermalEnergy(ScalarField2d("tmp", nodes1, eps1))
    nodes2.specificThermalEnergy(ScalarField2d("tmp", nodes2, eps2))

    # Set node velocities
    for (nodes, vx) in ((nodes1, v1),
                        (nodes2, v2)):
        for i in xrange(nodes.numInternalNodes):
            r = nodes.positions()[i]
            vy = perturbMag*sin(2*pi*perturbWaves*r.x) * max(0.0, 1.0 - 2.0/perturbHeight*abs(r.y))
            nodes.velocity()[i] = Vector2d(vx, vy)
    del nodes

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase2d()
output("db")
output("db.appendNodeList(nodes1)")
output("db.appendNodeList(nodes2)")
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
xPlane0 = Plane2d(xmin1, Vector2d( 1.0, 0.0))
xPlane1 = Plane2d(xmax1, Vector2d(-1.0, 0.0))
yPlane0 = Plane2d(xmin1, Vector2d(0.0,  1.0))
yPlane1 = Plane2d(xmax2, Vector2d(0.0, -1.0))
xbc0 = PeriodicBoundary2d(xPlane0, xPlane1)
xbc1 = ReflectingBoundary2d(xPlane0)
xbc2 = ReflectingBoundary2d(xPlane1)
ybc0 = ReflectingBoundary2d(yPlane0)
ybc1 = ReflectingBoundary2d(yPlane1)
for bc in (xbc0, ybc0, ybc1):
    hydro.appendBoundary(bc)
del bc

#-------------------------------------------------------------------------------
# Construct a time integrator, and add the physics packages.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
integrator.appendPhysicsPackage(hydro)
integrator.lastDt = dtMin
integrator.dtMin = dtMin
integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.rigorousBoundaries = rigorousBoundaries
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.valid()")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.rigorousBoundaries")

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            initializeMassDensity = (sumForMassDensity != Hydro2d.MassDensityType.IntegrateDensity))
output("control")

# Smooth the initial conditions.
if restoreCycle is not None:
    control.loadRestartFile(restoreCycle)
else:
    control.iterateIdealH()
    control.smoothState(smoothIters)

    # Just in case the intended pressure distribution got screwed up by summing the
    # mass density, reset it.
    for (nodes, P) in ((nodes1, P1),
                       (nodes2, P2)):
        for i in xrange(nodes.numInternalNodes):
            nodes.specificThermalEnergy()[i] = P/((gamma - 1.0)*nodes.massDensity()[i])
    del nodes

    control.dropRestartFile()
    dumpPhysicsState(integrator,
                     "KelvinHelmholtz-2d",
                     visitDir)
hstats(nodeSet)

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
dtSample = goalTime / nSample
while control.time() < goalTime:
    nextGoalTime = min(control.time() + dtSample, goalTime)
    control.advance(nextGoalTime, maxSteps)
    control.dropRestartFile()
    dumpPhysicsState(integrator,
                     "KelvinHelmholtz-2d",
                     visitDir)

#-------------------------------------------------------------------------------
# Plot graphics if requested.
#-------------------------------------------------------------------------------
if graphics:
    # Plot the postion and velocity fields.
    rPlot = plotNodePositions2d(db,
                                colorNodeLists = False,
                                colorDomains = False)

    rhoPlot = plotFieldList(db.fluidMassDensity,
                            xFunction = "%s.y",
                            plotStyle = "points",
                            winTitle = "Density")
    PPlot = plotFieldList(db.fluidPressure,
                          xFunction = "%s.y",
                          plotStyle = "points",
                          winTitle = "Pressure")
    vxPlot = plotFieldList(db.fluidVelocity,
                           xFunction = "%s.y",
                           yFunction = "%s.x",
                           plotStyle = "points",
                           winTitle = "v_x")
    vyPlot = plotFieldList(db.fluidVelocity,
                           xFunction = "%s.y",
                           yFunction = "%s.y",
                           plotStyle = "points",
                           winTitle = "v_y")
