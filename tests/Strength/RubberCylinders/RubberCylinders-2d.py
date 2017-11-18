#-------------------------------------------------------------------------------
# A pair of rubber cylinders which collide, bounce, and rebound.
# This is a test for the tensile instability.
#
# See Monaghan 2000, JCP, 159, 290.
#-------------------------------------------------------------------------------
from SolidSpheral import *
from SpheralTestUtilities import *
from findLastRestart import *
from SpheralVisitDump import dumpPhysicsState
from math import *

# Load the mpi module if we're parallel.
import loadmpi
mpi, procID, numProcs = loadmpi.loadmpi()

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("2-D bouncing rubber cylinders test.")

#-------------------------------------------------------------------------------
# Generic problem parameters
# All CGS units.
#-------------------------------------------------------------------------------
commandLine(
    SolidNodeListConstructor = SphSolidNodeList2d,

    # Should we do one or two distinct cylinders?
    twoCylinders = False,

    # Geometry
    r0 = 3.0,
    r1 = 4.0,
    x0 = 5.0,
    y0 = 0.0,

    # Numbers of nodes.
    nr = 10,

    # Inital velocity of the cylinder.
    vx0 = -0.059*8.52e4,

    # Node seeding stuff.
    seed = "lattice",
    nPerh = 2.01,

    # Material properties.
    rho0 = 1.01,
    c0 = 8.52e4,
    mu0 = 0.22 * 1.01 * 8.52e4**2,
    Y0 = 1.0e100,

    # Material specific bounds on the mass density.
    etamin = 0.5,
    etamax = 1.5,

    # Material specific limits on the smoothing scale.
    hminratio = 0.1,
    hmin = 1e-5,
    hmax = 0.5,

    # Hydro parameters.
    Qconstructor = MonaghanGingoldViscosity2d,
    Cl = 1.5,
    Cq = 1.5,
    Qlimiter = False,
    balsaraCorrection = False,
    epsilon2 = 1e-2,
    cfl = 0.5,
    XSPH = True,
    epsilonTensile = 0.0,
    nTensile = 4,
    HEvolution = Hydro2d.HEvolutionType.IdealH,
    sumForMassDensity = Hydro2d.MassDensityType.IntegrateDensity,
    compatibleEnergy = True,

    # Times, and simulation control.
    goalTime = 4000.0e-6,
    dtSample = 5e-6,
    dt = 1e-10,
    dtMin = 1e-10,
    dtMax = 1e-3,
    dtGrowth = 2.0,
    maxSteps = 200,
    statsStep = 10,
    redistributeStep = None,
    smoothIters = 0,

    # Restart and output files.
    restoreCycle = None,
    restartStep = 200,
    baseDir = "dumps-rubberCylinders-2d",
    )

# Derive some node number parameters.
dr = r1 - r0
nx = int(2.0*r1/dr * nr + 0.5)
ntheta = int(pi*(r0 + r1)/dr * nr + 0.5)

# Restart and output files.
dataDir = "%s/%s/XSPH=%s/compatibleEnergy=%s/nr=%i/nperh=%4.2f" % (baseDir,
                                                                   str(SolidNodeListConstructor).split("'")[1],
                                                                   XSPH,
                                                                   compatibleEnergy,
                                                                   nr,
                                                                   nPerh)
restartDir = dataDir + "/restarts/proc-%04i" % mpi.rank
visitDir = dataDir + "/visit"
restartBaseName = restartDir + "/RubberCylinders"

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if not os.path.exists(dataDir):
        os.makedirs(dataDir)
    if not os.path.exists(visitDir):
        os.makedirs(visitDir)
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
# Material properties.
#-------------------------------------------------------------------------------
ack = rho0*c0*c0
eos = LinearPolynomialEquationOfStateCGS2d(rho0,    # reference density  
                                           etamin,  # etamin             
                                           etamax,  # etamax             
                                           0.0,     # A0
                                           ack,     # A1
                                           0.0,     # A2
                                           0.0,     # A3
                                           0.0,     # B0
                                           0.0,     # B1
                                           0.0,     # B2
                                           55.350)  # atomic weight

strengthModel = ConstantStrength(mu0,
                                 Y0)

#-------------------------------------------------------------------------------
# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
#-------------------------------------------------------------------------------
WT = TableKernel2d(BSplineKernel2d(), 1000)
WTPi = TableKernel2d(BSplineKernel2d(), 1000)
output("WT")
output("WTPi")
kernelExtent = WT.kernelExtent()

#-------------------------------------------------------------------------------
# Create the NodeLists.
#-------------------------------------------------------------------------------
nodes1 = SolidNodeListConstructor("Rubber 1", eos, strengthModel, WT, WTPi)
nodes2 = SolidNodeListConstructor("Rubber 2", eos, strengthModel, WT, WTPi)
nodeSet = [nodes1]
if twoCylinders:
    nodeSet.append(nodes2)
for nodes in nodeSet:
    nodes.nodesPerSmoothingScale = nPerh
    nodes.epsilonTensile = epsilonTensile
    nodes.nTensile = nTensile
    nodes.hmin = hmin
    nodes.hmax = hmax
    nodes.hminratio = hminratio
    nodes.XSPH = XSPH
    nodes.rhoMin = etamin*rho0
    nodes.rhoMax = etamax*rho0
    output("nodes.name()")
    output("  nodes.nodesPerSmoothingScale")
    output("  nodes.epsilonTensile")
    output("  nodes.nTensile")
    output("  nodes.hmin")
    output("  nodes.hmax")
    output("  nodes.hminratio")
    output("  nodes.XSPH")
    output("  nodes.rhoMin")
    output("  nodes.rhoMax")
del nodes

#-------------------------------------------------------------------------------
# Construct the neighbor objects and associate them with the node lists.
#-------------------------------------------------------------------------------
neighbor1 = TreeNeighbor2d(nodes1,
                           kernelExtent = kernelExtent)
neighbor2 = TreeNeighbor2d(nodes2,
                           kernelExtent = kernelExtent)
nodes1.registerNeighbor(neighbor1)
nodes2.registerNeighbor(neighbor2)

#-------------------------------------------------------------------------------
# Set node properties (positions, velocites, etc.)
#-------------------------------------------------------------------------------
if restoreCycle is None:
    print "Generating node distribution."
    from GenerateNodeDistribution2d import *
    from ParMETISDistributeNodes import distributeNodes2d
    if seed == "lattice":
        ny = nx
    else:
        ny = ntheta
    generator1 = GenerateNodeDistribution2d(nx,
                                            ny,
                                            rho = rho0,
                                            distributionType = seed,
                                            xmin = (-r1, -r1),
                                            xmax = (r1, r1),
                                            rmin = r0,
                                            rmax = r1,
                                            nNodePerh = nPerh,
                                            theta = 2.0*pi,
                                            SPH = not isinstance(nodes1, AsphSolidNodeList2d))
    generator2 = GenerateNodeDistribution2d(nx,
                                            ny,
                                            rho = rho0,
                                            distributionType = seed,
                                            xmin = (-r1, -r1),
                                            xmax = (r1, r1),
                                            rmin = r0,
                                            rmax = r1,
                                            nNodePerh = nPerh,
                                            theta = 2.0*pi,
                                            SPH = not isinstance(nodes2, AsphSolidNodeList2d))

    # Displace the nodes to the correct centering.
    assert generator1.localNumNodes() == generator2.localNumNodes()
    for i in xrange(generator1.localNumNodes()):
        generator1.x[i] += x0
        generator1.y[i] += y0
        generator2.x[i] -= x0
        generator2.y[i] += y0

    print "Starting node distribution..."
    if twoCylinders:
        distributeNodes2d((nodes1, generator1),
                          (nodes2, generator2))
    else:
        distributeNodes2d((nodes1, generator1))
    for nodes in nodeSet:
        output("nodes.name()")
        output("    mpi.allreduce(nodes.numInternalNodes, mpi.MIN)")
        output("    mpi.allreduce(nodes.numInternalNodes, mpi.MAX)")
        output("    mpi.allreduce(nodes.numInternalNodes, mpi.SUM)")
    del nodes

    # Set node specific thermal energies
    print "Initial pressure for %s: %g" % (nodes1.name(),
                                           nodes1.equationOfState().pressure(rho0, 0.0))
    print "Initial pressure for %s: %g" % (nodes2.name(),
                                           nodes2.equationOfState().pressure(rho0, 0.0))

    # Set the projectile velocities.
    nodes1.velocity(VectorField2d("tmp", nodes1, Vector2d(vx0, 0.0)))
    nodes2.velocity(VectorField2d("tmp", nodes2, Vector2d(-vx0, 0.0)))

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node lists.
#-------------------------------------------------------------------------------
db = DataBase2d()
for nodes in nodeSet:
    db.appendNodeList(nodes)
del nodes
output("db")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the artificial viscosity.
#-------------------------------------------------------------------------------
q = Qconstructor(Cl, Cq)
q.limiter = Qlimiter
q.balsaraShearCorrection = balsaraCorrection
q.epsilon2 = epsilon2
output("q")
output("q.Cl")
output("q.Cq")
output("q.limiter")
output("q.epsilon2")
output("q.balsaraShearCorrection")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
hydro = Hydro2d(WT, WTPi, q, compatibleEnergy)
hydro.cfl = cfl
hydro.useVelocityMagnitudeForDt = True
hydro.HEvolution = HEvolution
hydro.sumForMassDensity = sumForMassDensity
hydro.HsmoothMin = hmin
hydro.HsmoothMax = hmin
output("hydro")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.useVelocityMagnitudeForDt")
output("hydro.HEvolution")
output("hydro.sumForMassDensity")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.valid()")

#-------------------------------------------------------------------------------
# Construct a strength physics object.
#-------------------------------------------------------------------------------
strength = Strength2d()
output("strength")

#-------------------------------------------------------------------------------
# Construct a predictor corrector integrator.
#-------------------------------------------------------------------------------
integrator = PredictorCorrectorIntegrator2d(db)
integrator.appendPhysicsPackage(hydro)
integrator.appendPhysicsPackage(strength)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.havePhysicsPackage(strength)")
output("integrator.valid()")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")

#-------------------------------------------------------------------------------
# Construct boundary conditions, and add them to our physics packages.
#-------------------------------------------------------------------------------
if not twoCylinders:
    xbcPlane = Plane2d(Vector2d(0.0, 0.0, 0.0), Vector2d(1.0, 0.0, 0.0))
    xbc = ReflectingBoundary2d(xbcPlane)
    for package in integrator.physicsPackages():
        package.appendBoundary(xbc)

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            redistributeStep = redistributeStep,
                            restartBaseName = restartBaseName,
                            initializeMassDensity = False)
output("control")

#-------------------------------------------------------------------------------
# Drop visualization files.
#-------------------------------------------------------------------------------
def viz(fields = [],
        filename = "RubberCylinders-2d"):
    dumpPhysicsState(integrator,
                     filename,
                     visitDir,
                     fields = fields,
                     )

#-------------------------------------------------------------------------------
# Smooth the initial conditions/restore state.
#-------------------------------------------------------------------------------
if restoreCycle is not None:
    control.loadRestartFile(restoreCycle)
    control.setRestartBaseName(restartBaseName)
    control.setFrequency(control.updateDomainDistribution, redistributeStep)

else:
    control.iterateIdealH()
    control.smoothState(smoothIters)
    control.dropRestartFile()
    viz()

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
while control.time() < goalTime:
    nextGoalTime = min(control.time() + dtSample, goalTime)
    control.advance(nextGoalTime, maxSteps)
    control.dropRestartFile()
    viz()
