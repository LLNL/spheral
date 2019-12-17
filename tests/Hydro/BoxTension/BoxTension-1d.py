#-------------------------------------------------------------------------------
# The Hydrostatic Equilibrium/Surface Tension Test
#-------------------------------------------------------------------------------
import shutil
from math import *
from SolidSpheral1d import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *

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
    nx1 = 50,
    nx2 = 50,

    nPerh = 1.51,
    order = 5,

    svph = False,
    crksph = False,
    psph = False,
    asph = False,
    filter = 0.0,   # For CRKSPH
    evolveTotalEnergy = False,
    HopkinsConductivity = False,
    solid = False,  # Try out the solid form of the hydro

    linearConsistent = False,
    fcentroidal = 0.0,
    fcellPressure = 0.0,
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

    volumeType = RKVoronoiVolume,
    densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
    correctionOrder = LinearOrder,
    compatibleEnergy = True,
    gradhCorrection = True,
    correctVelocityGradient = True,

    clearDirectories = False,
    restoreCycle = None,
    restartStep = 200,
    dataDir = "dumps-boxtension-1d",
    graphics = True,
    serialDump = False,
    )

# Decide on our hydro algorithm.
if svph:
    hydroname = "SVPH"
elif crksph:
    hydroname = "CRKSPH"
elif psph:
    hydroname = "PSPH"
else:
    hydroname = "SPH"
if asph:
    hydroname = "A" + hydroname

# Build our directory paths.
densityUpdateLabel = {IntegrateDensity : "IntegrateDensity",
                      SumDensity : "SumDensity",
                      RigorousSumDensity : "RigorousSumDensity",
                      SumVoronoiCellDensity : "SumVoronoiCellDensity"}
baseDir = os.path.join(dataDir,
                       hydroname,
                       densityUpdateLabel[densityUpdate],
                       "linearConsistent=%s" % linearConsistent,
                       "XSPH=%s" % XSPH,
                       "nPerh=%3.1f" % nPerh,
                       "fcentroidal=%1.3f" % fcentroidal,
                       "fcellPressure = %1.3f" % fcellPressure,
                       "%i" % (nx1 + nx2))
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
# Material properties.
#-------------------------------------------------------------------------------
mu = 1.0
eos1 = GammaLawGasMKS(gamma1, mu)
eos2 = GammaLawGasMKS(gamma1, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(NBSplineKernel(order), 1000)
WTPi = WT
output("WT")

#-------------------------------------------------------------------------------
# Make the NodeLists.
#-------------------------------------------------------------------------------
if solid:
    MNL = makeSolidNodeList
else:
    MNL = makeFluidNodeList
outerNodes1 = MNL("outer1", eos1,
                  hmin = hmin,
                  hmax = hmax,
                  hminratio = hminratio,
                  kernelExtent = WT.kernelExtent,
                  nPerh = nPerh)
outerNodes2 = MNL("outer2", eos1,
                  hmin = hmin,
                  hmax = hmax,
                  hminratio = hminratio,
                  kernelExtent = WT.kernelExtent,
                  nPerh = nPerh)
innerNodes = MNL("inner", eos2,
                 hmin = hmin,
                 hmax = hmax,
                 hminratio = hminratio,
                 kernelExtent = WT.kernelExtent,
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
distributeNodesInRange1d([(outerNodes1, [(nx1/2, rho1, (x0, x1))]),
                          (innerNodes,  [(nx2, rho2, (x1, x2))]),
                          (outerNodes2, [(nx1/2, rho1, (x2, x3))])])
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
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if svph:
    hydro = SVPH(dataBase = db,
                 W = WT, 
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
                 xmax = Vector(x3 + (x3 - x0), y3 + (y3 - y0)),
                 ASPH = asph)
elif crksph:
    hydro = CRKSPH(dataBase = db,
                   W = WT, 
                   filter = filter,
                   epsTensile = epsilonTensile,
                   nTensile = nTensile,
                   cfl = cfl,
                   compatibleEnergyEvolution = compatibleEnergy,
                   XSPH = XSPH,
                   correctionOrder = correctionOrder,
                   volumeType = volumeType,
                   densityUpdate = densityUpdate,
                   HUpdate = HUpdate,
                   ASPH = asph)
elif psph:
    hydro = PSPH(dataBase = db,
                 W = WT,
                 Q = q,
                 filter = filter,
                 cfl = cfl,
                 compatibleEnergyEvolution = compatibleEnergy,
                 evolveTotalEnergy = evolveTotalEnergy,
                 HopkinsConductivity = HopkinsConductivity,
                 densityUpdate = densityUpdate,
                 correctVelocityGradient = correctVelocityGradient,
                 HUpdate = HUpdate,
                 XSPH = XSPH,
                 ASPH = asph)

else:
    hydro = SPH(dataBase = db,
                W = WT,
                cfl = cfl,
                compatibleEnergyEvolution = compatibleEnergy,
                gradhCorrection = gradhCorrection,
                correctVelocityGradient = correctVelocityGradient,
                XSPH = XSPH,
                densityUpdate = densityUpdate,
                HUpdate = HUpdate,
                epsTensile = epsilonTensile,
                nTensile = nTensile,
                ASPH = asph)

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
                            skipInitialPeriodicWork = svph,
                            SPH = (not asph))
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
procs = mpi.procs
rank = mpi.rank
serialData = []
pError = []
rhoError = []
xdata = []
velError = []
i,j = 0,0
from Pnorm import Pnorm
for i in xrange(procs):
    for (nodeL, gamma, rho, P) in ((outerNodes1, gamma1, rho1, P1),
                                   (outerNodes2, gamma1, rho1, P1),
                                   (innerNodes, gamma2, rho2, P2)):
      if rank == i:
        Pfield = ScalarField("pressure", nodeL)
        nodeL.pressure(Pfield)
        for j in xrange(nodeL.numInternalNodes):
          pError.append(Pfield[j]-P)
          rhoError.append(nodeL.massDensity()[j] - rho)
          velError.append(nodeL.velocity()[j][0]) #Velcoity is supposed to be zero
          xdata.append(nodeL.positions()[j][0])
          if serialDump:
            serialData.append([nodeL.positions()[j],3.0/(nodeL.Hfield()[j].Trace()),nodeL.mass()[j],nodeL.massDensity()[j],nodeL.specificThermalEnergy()[j],Pfield[j]])
serialData = mpi.reduce(serialData,mpi.SUM)
pError = mpi.reduce(pError,mpi.SUM)
rhoError = mpi.reduce(rhoError,mpi.SUM)
velError = mpi.reduce(velError,mpi.SUM)
xdata = mpi.reduce(xdata,mpi.SUM)
if rank == 0 and serialDump:
    f = open(os.path.join(dataDir, "./serialDump.ascii"),'w')
    for i in xrange(len(serialData)):
      f.write("{0} {1} {2} {3} {4} {5} {6}\n".format(i,serialData[i][0][0],serialData[i][1],serialData[i][2],serialData[i][3],serialData[i][4], serialData[i][5]))
    f.close()
if rank == 0:
 print len(pError)
 print "Pressure results: L1 error = %g, L2 Error = %g, L inf Error = %g \n" % (Pnorm(pError, xdata).pnorm(1), Pnorm(pError, xdata).pnorm(2),Pnorm(pError, xdata).pnorm("inf"))
 print "Density results: L1 error = %g, L2 Error = %g, L inf Error = %g \n" % (Pnorm(rhoError, xdata).pnorm(1), Pnorm(rhoError, xdata).pnorm(2),Pnorm(rhoError, xdata).pnorm("inf"))
 print "Velocity results: L1 error = %g, L2 Error = %g, L inf Error = %g \n" % (Pnorm(velError, xdata).pnorm(1), Pnorm(velError, xdata).pnorm(2),Pnorm(velError, xdata).pnorm("inf"))
#-------------------------------------------------------------------------------
# Plot the results.
#-------------------------------------------------------------------------------
if graphics:
    from SpheralGnuPlotUtilities import *

    rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db, plotGhosts=True)
    pE = plotEHistory(control.conserve)
    massPlot = plotFieldList(db.fluidMass,
                             winTitle = "mass",
                             plotStyle = "points",
                             colorNodeLists = False)

    if crksph:
        APlot = plotFieldList(hydro.A(),
                              winTitle = "A",
                              plotStyle = "points",
                              colorNodeLists = False)
        BPlot = plotFieldList(hydro.B(),
                              yFunction = "%s.x",
                              winTitle = "B",
                              plotStyle = "points",
                              colorNodeLists = False)
        volPlot = plotFieldList(hydro.volume(),
                                winTitle = "volume",
                                plotStyle = "points",
                                colorNodeLists = False)
        spPlot = plotFieldList(hydro.surfacePoint(),
                              winTitle = "surface point",
                              plotStyle = "points",
                              colorNodeLists = False)

    cs = db.newFluidScalarFieldList(0.0, "sound speed")
    db.fluidSoundSpeed(cs)
    csPlot = plotFieldList(cs,
                           plotStyle = "points",
                           winTitle="Sound speed")
