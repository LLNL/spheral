#-------------------------------------------------------------------------------
# The Hydrostatic Equilibrium/Surface Tension Test
#-------------------------------------------------------------------------------
import shutil
from math import *
from Spheral1d import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from findLastRestart import *

import mpi
from DistributeNodes import distributeNodesInRange1d

title("1-D integrated hydro test --  Hydrostatic Equilibrium/Surface Tension Test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(
    # Outer state.
    rho1 = 1.0,
    P1 = 1.0,
    gamma1 = 1.5,

    # Inner state
    rho2 = 4.0,
    P2 = 1.0,
    gamma2 = 1.5,

    # Geometry 
    x0 = 0.0,
    x1 = 0.25,
    x2 = 0.75,
    x3 = 1.0,

    #Translation
    vel=0.0,

    # Resolution and node seeding.
    nx1 = 100,
    ny1 = 100,

    nx2 = 50,
    ny2 = 50,

    nPerh = 1.51,

    SVPH = False,
    CRKSPH = False,
    ASPH = False,
    SPH = True,   # This just chooses the H algorithm -- you can use this with CRKSPH for instance.
    filter = 0.0,  # For CRKSPH
    Qconstructor = MonaghanGingoldViscosity,
    #Qconstructor = TensorMonaghanGingoldViscosity,
    boolReduceViscosity = False,
    nh = 5.0,
    aMin = 0.1,
    aMax = 2.0,
    linearConsistent = False,
    fcentroidal = 0.0,
    fcellPressure = 0.0,
    Cl = 1.0, 
    Cq = 0.75,
    Qlimiter = False,
    balsaraCorrection = True,
    epsilon2 = 1e-2,
    hmin = 1e-5,
    hmax = 0.5,
    hminratio = 0.1,
    cfl = 0.5,
    XSPH = True,
    epsilonTensile = 0.0,
    nTensile = 8,

    IntegratorConstructor = CheapSynchronousRK2Integrator,
    goalTime = 7.0,
    steps = None,
    vizCycle = 5,
    vizTime = 0.1,
    dt = 0.0001,
    dtMin = 1.0e-5, 
    dtMax = 0.1,
    dtGrowth = 2.0,
    maxSteps = None,
    statsStep = 10,
    HUpdate = IdealH,
    domainIndependent = False,
    rigorousBoundaries = False,
    dtverbose = False,

    densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
    compatibleEnergy = True,
    gradhCorrection = False,

    clearDirectories = False,
    restoreCycle = None,
    restartStep = 200,
    dataDir = "dumps-boxtension-1d",
    graphics = True,
    )

# Decide on our hydro algorithm.
if SVPH:
    if ASPH:
        HydroConstructor = ASVPHFacetedHydro
    else:
        HydroConstructor = SVPHFacetedHydro
elif CRKSPH:
    if ASPH:
        HydroConstructor = ACRKSPHHydro
    else:
        HydroConstructor = CRKSPHHydro
else:
    if ASPH:
        HydroConstructor = ASPHHydro
    else:
        HydroConstructor = SPHHydro

# Build our directory paths.
densityUpdateLabel = {IntegrateDensity : "IntegrateDensity",
                      SumDensity : "SumDensity",
                      RigorousSumDensity : "RigorousSumDensity",
                      SumVoronoiCellDensity : "SumVoronoiCellDensity"}
baseDir = os.path.join(dataDir,
                       HydroConstructor.__name__,
                       Qconstructor.__name__,
                       densityUpdateLabel[densityUpdate],
                       "linearConsistent=%s" % linearConsistent,
                       "XSPH=%s" % XSPH,
                       "nPerh=%3.1f" % nPerh,
                       "fcentroidal=%1.3f" % fcentroidal,
                       "fcellPressure = %1.3f" % fcellPressure,
                       "%ix%i" % (nx1 + nx2, ny1 + ny2))
restartDir = os.path.join(baseDir, "restarts")
restartBaseName = os.path.join(restartDir, "boxtension-1d-%i" % (nx1 + nx2))

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if clearDirectories and os.path.exists(baseDir):
        shutil.rmtree(baseDir)
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# If we're restarting, find the set of most recent restart files.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    restoreCycle = findLastRestart(restartBaseName)

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
mu = 1.0
eos1 = GammaLawGasMKS(gamma1, mu)
eos2 = GammaLawGasMKS(gamma1, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(BSplineKernel(), 1000)
WTPi = TableKernel(BSplineKernel(), 1000)
output("WT")
output("WTPi")
kernelExtent = WT.kernelExtent

#-------------------------------------------------------------------------------
# Make the NodeLists.
#-------------------------------------------------------------------------------
outerNodes1 = makeFluidNodeList("outer1", eos1,
                                hmin = hmin,
                                hmax = hmax,
                                hminratio = hminratio,
                                nPerh = nPerh)
outerNodes2 = makeFluidNodeList("outer2", eos1,
                                hmin = hmin,
                                hmax = hmax,
                                hminratio = hminratio,
                                nPerh = nPerh)
innerNodes = makeFluidNodeList("inner", eos2,
                               hmin = hmin,
                               hmax = hmax,
                               hminratio = hminratio,
                               nPerh = nPerh)
nodeSet = (outerNodes1, outerNodes2, innerNodes)
for nodes in nodeSet:
    output("nodes.name")
    output("    nodes.hmin")
    output("    nodes.hmax")
    output("    nodes.hminratio")
    output("    nodes.nodesPerSmoothingScale")
del nodes

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    distributeNodesInRange1d([(outerNodes1, [(nx1/4, rho1, (x0, x1))]),
                              (innerNodes,  [(nx2, rho2, (x1, x2))]),
                              (outerNodes2, [(nx1/4, rho1, (x2, x3))])])
    for nodes in nodeSet:
        print nodes.name, ":"
        output("    mpi.reduce(nodes.numInternalNodes, mpi.MIN)")
        output("    mpi.reduce(nodes.numInternalNodes, mpi.MAX)")
        output("    mpi.reduce(nodes.numInternalNodes, mpi.SUM)")
    del nodes

    # Set node specific thermal energies
    for (nodes, gamma, rho, P) in ((outerNodes1, gamma1, rho1, P1),
                                   (outerNodes2, gamma1, rho1, P1),
                                   (innerNodes, gamma2, rho2, P2)):
        eps0 = P/((gamma - 1.0)*rho)
        nodes.specificThermalEnergy(ScalarField("tmp", nodes, eps0))
        vels = nodes.velocity()
        for i in xrange(nodes.numInternalNodes):
           vels[i]=Vector(vel)
    del nodes

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node lists
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
for nodes in nodeSet:
    db.appendNodeList(nodes)
del nodes
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
if SVPH:
    hydro = HydroConstructor(WT, q,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             densityUpdate = densityUpdate,
                             XSVPH = XSPH,
                             linearConsistent = linearConsistent,
                             generateVoid = False,
                             HUpdate = HUpdate,
                             fcentroidal = fcentroidal,
                             fcellPressure = fcellPressure,
                             xmin = Vector(x0 - (x3 - x0), y0 - (y3 - y0)),
                             xmax = Vector(x3 + (x3 - x0), y3 + (y3 - y0)))
elif CRKSPH:
    hydro = HydroConstructor(WT, WTPi, q,
                             filter = filter,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             XSPH = XSPH,
                             densityUpdate = densityUpdate,
                             HUpdate = HUpdate)
else:
    hydro = HydroConstructor(WT,
                             WTPi,
                             q,
                             cfl = cfl,
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
# Construct the MMRV physics object.
#-------------------------------------------------------------------------------

if boolReduceViscosity:
    #q.reducingViscosityCorrection = True
    evolveReducingViscosityMultiplier = MorrisMonaghanReducingViscosity(q,nh,aMin,aMax)
    
    packages.append(evolveReducingViscosityMultiplier)


#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane0 = Plane(Vector(x0), Vector( 1.0))
xPlane1 = Plane(Vector(x3), Vector(-1.0))

xbc = PeriodicBoundary(xPlane0, xPlane1)

xbc0 = ReflectingBoundary(xPlane0)
xbc1 = ReflectingBoundary(xPlane1)

bcSet = [xbc]
#bcSet = [xbc0, xbc1]

for p in packages:
    for bc in bcSet:
        p.appendBoundary(bc)

#-------------------------------------------------------------------------------
# Construct a time integrator, and add the physics packages.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
for p in packages:
    integrator.appendPhysicsPackage(p)
integrator.cullGhostNodes = False
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
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            skipInitialPeriodicWork = (HydroConstructor in (SVPHFacetedHydro, ASVPHFacetedHydro)),
                            SPH = SPH)
output("control")

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)
else:
    control.advance(goalTime, maxSteps)

print "Energy conservation: original=%g, final=%g, error=%g" % (control.conserve.EHistory[0],
                                                                control.conserve.EHistory[-1],
                                                                (control.conserve.EHistory[-1] - control.conserve.EHistory[0])/control.conserve.EHistory[0])

#-------------------------------------------------------------------------------
# Plot the results.
#-------------------------------------------------------------------------------
if graphics:
    from SpheralGnuPlotUtilities import *

    rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db, plotGhosts=True)
    pE = plotEHistory(control.conserve)

    if CRKSPH:
        APlot = plotFieldList(hydro.A(),
                              winTitle = "A",
                              colorNodeLists = False)
        BPlot = plotFieldList(hydro.B(),
                              yFunction = "%s.x",
                              winTitle = "B",
                              colorNodeLists = False)

    cs = db.newFluidScalarFieldList(0.0, "sound speed")
    db.fluidSoundSpeed(cs)
    csPlot = plotFieldList(cs, winTitle="Sound speed")
