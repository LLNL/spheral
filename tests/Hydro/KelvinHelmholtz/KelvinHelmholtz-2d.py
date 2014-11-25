#-------------------------------------------------------------------------------
# This is the basic Kelvin-Helmholtz problem as discussed in
# Springel 2010, MNRAS, 401, 791-851.
#-------------------------------------------------------------------------------
import shutil
from math import *
from Spheral2d import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from findLastRestart import *
from GenerateNodeDistribution2d import *
from CompositeNodeDistribution import *
from CentroidalVoronoiRelaxation import *

import mpi
import DistributeNodes

title("Kelvin-Helmholtz test problem in 2D")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(nx1 = 100,
            ny1 = 50,
            nx2 = 100,
            ny2 = 50,
            
            rho1 = 2.0,
            rho2 = 1.0,
            P1 = 2.5,
            P2 = 2.5,
            vx1 = 0.5,
            vx2 = -0.5,
            freq = 4.0,
            w0 = 0.1,
            sigma = 0.05/sqrt(2.0),

            numNodeLists = 2,  # If 2, makes this a two material problem.

            gamma = 5.0/3.0,
            mu = 1.0,

            nPerh = 1.51,

            SVPH = False,
            CSPH = False,
            ASPH = False,
            SPH = True,   # This just chooses the H algorithm -- you can use this with CSPH for instance.
            filter = 0.1,   # CSPH filtering
            momentumConserving = True, # For CSPH
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
            Cq = 0.75,
            Qlimiter = False,
            balsaraCorrection = True,
            epsilon2 = 1e-2,
            hmin = 0.0001, 
            hmax = 0.5,
            hminratio = 0.1,
            cfl = 0.5,
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
            rigorousBoundaries = True,
            dtverbose = False,

            densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
            compatibleEnergy = False,           # <--- Important!  rigorousBoundaries does not work with the compatibleEnergy algorithm currently.
            gradhCorrection = False,

            useVoronoiOutput = True,
            clearDirectories = False,
            restoreCycle = None,
            restartStep = 100,
            redistributeStep = 500,
            checkRestart = False,
            dataDir = "dumps-KelvinHelmholtz-2d",
            outputFile = "None",
            comparisonFile = "None",
            )

assert numNodeLists in (1, 2)

# Decide on our hydro algorithm.
if SVPH:
    if ASPH:
        HydroConstructor = ASVPHFacetedHydro
    else:
        HydroConstructor = SVPHFacetedHydro
elif CSPH:
    if ASPH:
        HydroConstructor = ACSPHHydro
    else:
        HydroConstructor = CSPHHydro
else:
    if ASPH:
        HydroConstructor = ASPHHydro
    else:
        HydroConstructor = SPHHydro

dataDir = os.path.join(dataDir,
                       "rho1=%g-rho2=%g" % (rho1, rho2),
                       "vx1=%g-vx2=%g" % (abs(vx1), abs(vx2)),
                       str(HydroConstructor).split("'")[1].split(".")[-1],
                       "densityUpdate=%s" % (densityUpdate),
                       "XSPH=%s" % XSPH,
                       "filter=%s" % filter,
                       "%s-Cl=%g-Cq=%g" % (str(Qconstructor).split("'")[1].split(".")[-1], Cl, Cq),
                       "%ix%i" % (nx1, ny1 + ny2),
                       "nPerh=%g-Qhmult=%g" % (nPerh, Qhmult))
restartDir = os.path.join(dataDir, "restarts")
vizDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, "KelvinHelmholtz-2d")
vizBaseName = "KelvinHelmholtz-2d"

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
nodes2 = makeFluidNodeList("Low density gas", eos,
                           hmin = hmin,
                           hmax = hmax,
                           hminratio = hminratio,
                           nPerh = nPerh)
nodeSet = [nodes1, nodes2]
for nodes in nodeSet:
    output("nodes.name")
    output("nodes.hmin")
    output("nodes.hmax")
    output("nodes.hminratio")
    output("nodes.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    generator1 = GenerateNodeDistribution2d(nx1, ny1,
                                            rho = rho1,
                                            distributionType = "lattice",
                                            xmin = (0.0,  0.25),
                                            xmax = (1.0,  0.75),
                                            nNodePerh = nPerh,
                                            SPH = SPH)
    generator21 = GenerateNodeDistribution2d(nx2, int(0.5*ny2 + 0.5),
                                             rho = rho2,
                                             distributionType = "lattice",
                                             xmin = (0.0, 0.0),
                                             xmax = (1.0, 0.25),
                                             nNodePerh = nPerh,
                                             SPH = SPH)
    generator22 = GenerateNodeDistribution2d(nx2, int(0.5*ny2 + 0.5),
                                             rho = rho2,
                                             distributionType = "lattice",
                                             xmin = (0.0, 0.75),
                                             xmax = (1.0, 1.0),
                                             nNodePerh = nPerh,
                                             SPH = SPH)
    generator2 = CompositeNodeDistribution(generator21, generator22)

    if mpi.procs > 1:
        from VoronoiDistributeNodes import distributeNodes2d
    else:
        from DistributeNodes import distributeNodes2d

    if numNodeLists == 2:
        distributeNodes2d((nodes1, generator1),
                          (nodes2, generator2))
    else:
        gen = CompositeNodeDistribution(generator1, generator2)
        distributeNodes2d((nodes1, gen))
        

    # A helpful method for setting y velocities.
    def vy(ri):
        thpt = 1.0/(2.0*sigma*sigma)
        return (w0*sin(freq*pi*ri.x) *
                (exp(-((ri.y - 0.25)**2 * thpt)) +
                 exp(-((ri.y - 0.75)**2 * thpt))))*abs(0.5 - ri.y)

    # Finish initial conditions.
    eps1 = P1/((gamma - 1.0)*rho1)
    eps2 = P2/((gamma - 1.0)*rho2)
    if numNodeLists == 2:
        nodes1.specificThermalEnergy(ScalarField("tmp", nodes1, eps1))
        nodes2.specificThermalEnergy(ScalarField("tmp", nodes2, eps2))
        for (nodes, vx) in ((nodes1, vx1),
                            (nodes2, vx2)):
            pos = nodes.positions()
            vel = nodes.velocity()
            for i in xrange(nodes.numInternalNodes):
                vel[i] = Vector(vx, vy(pos[i]))
    else:
        pos = nodes1.positions()
        vel = nodes1.velocity()
        eps = nodes1.specificThermalEnergy()
        for i in xrange(nodes1.numInternalNodes):
            if pos[i].y > 0.25 and pos[i].y < 0.75:
                eps[i] = eps1
                vel[i] = Vector(vx1, vy(pos[i]))
            else:
                eps[i] = eps2
                vel[i] = Vector(vx2, vy(pos[i]))

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
for nodes in nodeSet:
    db.appendNodeList(nodes)
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
                             xmin = Vector(-2.0, -2.0),
                             xmax = Vector(3.0, 3.0))
                             # xmin = Vector(x0 - 0.5*(x2 - x0), y0 - 0.5*(y2 - y0)),
                             # xmax = Vector(x2 + 0.5*(x2 - x0), y2 + 0.5*(y2 - y0)))
elif CSPH:
    hydro = HydroConstructor(WT, WTPi, q,
                             filter = filter,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             XSPH = XSPH,
                             densityUpdate = densityUpdate,
                             HUpdate = HUpdate,
                             momentumConserving = momentumConserving)
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
    evolveReducingViscosityMultiplier = MorrisMonaghanReducingViscosity(q,nh,aMin,aMax)
    
    packages.append(evolveReducingViscosityMultiplier)

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xp1 = Plane(Vector(0.0, 0.0), Vector( 1.0, 0.0))
xp2 = Plane(Vector(1.0, 0.0), Vector(-1.0, 0.0))
yp1 = Plane(Vector(0.0, 0.0), Vector(0.0,  1.0))
yp2 = Plane(Vector(0.0, 1.0), Vector(0.0, -1.0))
xbc = PeriodicBoundary(xp1, xp2)
ybc = PeriodicBoundary(yp1, yp2)
ybc1 = ReflectingBoundary(yp1)
ybc2 = ReflectingBoundary(yp2)
bcSet = [xbc, ybc1, ybc2]

for p in packages:
    for bc in bcSet:
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

# Blago!  Currently a problem with periodic boundaries.
# integrator.cullGhostNodes = False

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
    import SpheralVisitDump
    vizMethod = SpheralVisitDump.dumpPhysicsState
control = SpheralController(integrator, WT,
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
if not steps is None:
    control.step(steps)

else:
    control.advance(goalTime, maxSteps)
    control.updateViz(control.totalSteps, integrator.currentTime, 0.0)
    control.dropRestartFile()
