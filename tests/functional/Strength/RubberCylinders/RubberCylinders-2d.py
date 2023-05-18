#-------------------------------------------------------------------------------
# A pair of rubber cylinders which collide, bounce, and rebound.
# This is a test for the tensile instability.
#
# See Monaghan 2000, JCP, 159, 290.
#-------------------------------------------------------------------------------
from Spheral2d import *
from SpheralTestUtilities import *
from GenerateNodeDistribution2d import *
from PeanoHilbertDistributeNodes import distributeNodes2d

import mpi, os, sys, shutil
from math import *


#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("2-D bouncing rubber cylinders test.")

#-------------------------------------------------------------------------------
# Generic problem parameters
# All CGS units.
#-------------------------------------------------------------------------------
commandLine(

    # Should we do one or two distinct cylinders?
    twoCylinders = True,

    # Geometry
    r0 = 3.0,
    r1 = 4.0,
    x0 = 5.0,
    y0 = 0.0,

    # Numbers of nodes.
    nr = 10,
    seed = "lattice",

    # Inital velocity of the cylinder.
    vx0 = -0.059*8.52e4,

    # Node seeding stuff.
    nPerh = 4.01,
    hminratio = 0.1,
    hmin = 1e-5,
    hmax = 0.5,
    HUpdate = IdealH,

    # Material properties.
    rho0 = 1.01,
    c0 = 8.52e4,
    mu0 = 0.22 * 1.01 * 8.52e4**2,
    Y0 = 1.0e100,

    # Material specific bounds on the mass density.
    etamin = 0.5,
    etamax = 1.5,

    # Hydro
    fsisph=False,
    crksph=False,

    #hydro options
    xsph = False,
    asph = False,
    densityUpdate = IntegrateDensity,
    compatibleEnergy = True,
    evolveTotalEnergy = False,
    gradhCorrection = True,
    correctVelocityGradient = True,
    epsilonTensile = 0.0,
    nTensile = 4,

    # fsi options 
    fsiRhoStabCoeff = 0.0, 
    fsiEpsDiffuseCoeff = 0.0, 
    fsiXSPHCoeff = 1.0,

    #crk options
    correctionOrder = LinearOrder,

    # artificial viscosity
    Qconstructor = LimitedMonaghanGingoldViscosity,
    Cl = 1.0,
    Cq = 2.0,
    Qlimiter = False,
    balsaraCorrection = False,
    epsilon2 = 1e-2,
    
    # Times, and simulation control.
    cfl = 0.25,
    goalTime = 6000.0e-6,
    dtSample = 5e-6,
    dt = 1e-10,
    dtMin = 1e-10,
    dtMax = 1e-3,
    dtGrowth = 2.0,
    steps = None,
    maxSteps = None,
    statsStep = 10,
    redistributeStep = None,
    smoothIters = 0,
    domainIndependent=False,
    dtverbose = False,

    # Outputs.
    clearDirectories=False,
    vizTime = 100.0e-6,
    vizCycle = None,
    vizDerivs = False,
    restoreCycle = None,
    restartStep = 200,
    baseDir = "dumps-rubberCylinders-2d",
    )

assert not (compatibleEnergy and evolveTotalEnergy)
assert not (fsisph and crksph)

# Derive some node number parameters.
dr = r1 - r0
nx = int(2.0*r1/dr * nr + 0.5)
ntheta = int(pi*(r0 + r1)/dr * nr + 0.5)

hydroname = "SPH"
if crksph:
    hydroname = "CRK"+hydroname
elif fsisph:
    hydroname = "FSI"+hydroname
if asph:
    hydroname = "A"+hydroname
# Restart and output files.
dataDir = os.path.join(baseDir,
                       hydroname,
                        "ntensile=%s" % nTensile,
                        "epstensile=%s" % epsilonTensile,
                        "nr=%s" % nr,
                        "nperh=%s" % nPerh)

restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "RubberCylinders-2d-%i" % (nr))
vizDir = os.path.join(dataDir, "visit")
if vizTime is None and vizCycle is None:
    vizBaseName = None
else:
    vizBaseName =  "RubberCylinders-2d-%i" % (nr)
#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
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
ack = rho0*c0*c0
eos = LinearPolynomialEquationOfStateCGS(rho0,    # reference density  
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
WT = TableKernel2d(WendlandC2Kernel(), 1000)
WTPi = TableKernel2d(WendlandC2Kernel(), 1000)
output("WT")
output("WTPi")
kernelExtent = WT.kernelExtent

#-------------------------------------------------------------------------------
# Create the NodeLists.
#-------------------------------------------------------------------------------
nodes1 = makeSolidNodeList("Rubber 1", eos, strengthModel, 
                               hmin = hmin,
                               hmax = hmax,
                               kernelExtent = kernelExtent,
                               hminratio = hminratio,
                               nPerh = nPerh)

nodes2 = makeSolidNodeList("Rubber 2", eos, strengthModel, 
                               hmin = hmin,
                               hmax = hmax,
                               kernelExtent = kernelExtent,
                               hminratio = hminratio,
                               nPerh = nPerh)
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
    nodes.rhoMin = etamin*rho0
    nodes.rhoMax = etamax*rho0
    output("nodes.name")
    output("  nodes.nodesPerSmoothingScale")
    output("  nodes.hmin")
    output("  nodes.hmax")
    output("  nodes.hminratio")
    output("  nodes.rhoMin")
    output("  nodes.rhoMax")

#-------------------------------------------------------------------------------
# Set node properties (positions, velocites, etc.)
#-------------------------------------------------------------------------------
if restoreCycle is None:
    print("Generating node distribution.")
    
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
                                            SPH = not asph)
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
                                            SPH = not asph)

    # Displace the nodes to the correct centering.
    assert generator1.localNumNodes() == generator2.localNumNodes()
    for i in range(generator1.localNumNodes()):
        generator1.x[i] += x0
        generator1.y[i] += y0
        generator2.x[i] -= x0
        generator2.y[i] += y0

    print("Starting node distribution...")
    if twoCylinders:
        distributeNodes2d((nodes1, generator1),
                          (nodes2, generator2))
    else:
        distributeNodes2d((nodes1, generator1))
    for nodes in nodeSet:
        output("nodes.name")
        output("    mpi.allreduce(nodes.numInternalNodes, mpi.MIN)")
        output("    mpi.allreduce(nodes.numInternalNodes, mpi.MAX)")
        output("    mpi.allreduce(nodes.numInternalNodes, mpi.SUM)")

    # Set node specific thermal energies
    print("Initial pressure for %s: %g" % (nodes1.name,
                                           nodes1.eos.pressure(rho0, 0.0)))
    print("Initial pressure for %s: %g" % (nodes2.name,
                                           nodes2.eos.pressure(rho0, 0.0)))

    # Set the projectile velocities.
    nodes1.velocity(VectorField2d("tmp", nodes1, Vector2d(vx0, 0.0)))
    nodes2.velocity(VectorField2d("tmp", nodes2, Vector2d(-vx0, 0.0)))

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node lists.
#-------------------------------------------------------------------------------
db = DataBase()
for nodes in nodeSet:
    db.appendNodeList(nodes)
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
if crksph:
    hydro = CRKSPH(dataBase = db,
                   order = correctionOrder,
                   cfl = cfl,
                   compatibleEnergyEvolution = compatibleEnergy,
                   XSPH = xsph,
                   densityUpdate = densityUpdate,
                   HUpdate = HUpdate,
                   ASPH = asph)
elif fsisph:
    hydro = FSISPH(dataBase = db,
                Q=q,
                W = WT,
                cfl = cfl,
                densityStabilizationCoefficient = fsiRhoStabCoeff, 
                specificThermalEnergyDiffusionCoefficient = fsiEpsDiffuseCoeff,  
                correctVelocityGradient = correctVelocityGradient,
                compatibleEnergyEvolution = compatibleEnergy,
                evolveTotalEnergy = evolveTotalEnergy,
                HUpdate = HUpdate,
                ASPH = asph,
                xsphCoefficient = fsiXSPHCoeff,
                epsTensile = epsilonTensile,
                nTensile = nTensile)
else:
    hydro = SPH(dataBase = db,
                Q=q,
                W = WT,
                cfl = cfl,
                compatibleEnergyEvolution = compatibleEnergy,
                evolveTotalEnergy = evolveTotalEnergy,
                gradhCorrection = gradhCorrection,
                correctVelocityGradient = correctVelocityGradient,
                densityUpdate = densityUpdate,
                HUpdate = HUpdate,
                XSPH = xsph,
                epsTensile = epsilonTensile,
                nTensile = nTensile,
                ASPH = asph)

output("hydro")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.useVelocityMagnitudeForDt")
output("hydro.HEvolution")
output("hydro.densityUpdate")


#-------------------------------------------------------------------------------
# Construct a predictor corrector integrator.
#-------------------------------------------------------------------------------
integrator = CheapSynchronousRK2Integrator(db)
integrator.appendPhysicsPackage(hydro)
integrator.lastDt = dt
integrator.dtMin = dtMin
integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.domainDecompositionIndependent = domainIndependent
integrator.verbose = dtverbose
integrator.cullGhostNodes = False
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
# Construct boundary conditions, and add them to our physics packages.
#-------------------------------------------------------------------------------
if not twoCylinders:
    xbcPlane = Plane(Vector2d(0.0, 0.0, 0.0), Vector2d(1.0, 0.0, 0.0))
    xbc = ReflectingBoundary(xbcPlane)
    for package in integrator.physicsPackages():
        package.appendBoundary(xbc)

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            iterateInitialH=True,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            redistributeStep = redistributeStep,
                            restartBaseName = restartBaseName,
                            vizBaseName = vizBaseName,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime,
                            vizDerivs = vizDerivs)
output("control")

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)
else:
    control.advance(goalTime, maxSteps)
    control.dropRestartFile()

