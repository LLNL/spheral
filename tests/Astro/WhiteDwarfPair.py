from math import *
from Spheral import *
from Spheral3d import *
from AsciiFileNodeGenerator import *
from VoronoiDistributeNodes import distributeNodes3d as distributeNodes
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from findLastRestart import *
from SpheralVisitDump import *
import mpi

title("White Dwarf pair test from ic readin file")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(HydroConstructor = SPHHydro,
            Qconstructor = MonaghanGingoldViscosity3d,
            #Qconstructor = TensorMonaghanGingoldViscosity3d,
            Cl = 1.0,
            Cq = 1.5,
            Qlimiter = False,
            balsaraCorrection = False,
            epsilon2 = 1e-2,
            epsilonTensile = 0.0,
            nTensile = 8,
            negligibleSoundSpeed = 1e-5,
            csMultiplier = 1e-4,
            energyMultiplier = 0.1,
            hmin = 1e-5,
            hmax = 1e20,
            hminratio = 0.05,
            HsmoothFraction = 0.0,
            cfl = 0.5,
            XSPH = True,
            
            plummerLength = 0.1e3,        # (cm) Plummer softening scale
            opening = 1.0,                 # (dimensionless, OctTreeGravity) opening parameter for tree walk
            fdt = 0.1,                     # (dimensionless, OctTreeGravity) timestep multiplier
            
            timeStepChoice = AccelerationRatio,
            myIntegrator = CheapSynchronousRK2Integrator3d,
            
            #neighborSearchType = Neighbor3d.NeighborSearchType.GatherScatter,
            numGridLevels = 20,
            topGridCellSize = 2.0,
            origin = Vector3d(0.0, 0.0),
            
            goalTime = 20.0,
            realGoalTime = 30.0,
            dt = 0.0001,
            dtMin = 1.0e-5,
            dtMax = None,
            dtGrowth = 2.0,
            dtSample = 1,
            rigorousBoundaries = False,
            maxSteps = None,
            statsStep = 1,
            smoothIters = 0,
            HEvolution = IdealH,
            densityUpdate = RigorousSumDensity,
            compatibleEnergy = True,
            gradhCorrection = False,
            verbosedt = False,
            restoreCycle = None,
            restartStep = 1000
            )

units = CGS()

dataDir = "dumps-wd"
restartDir = dataDir + "/restarts"
visitDir = dataDir + "/visit"
restartBaseName = restartDir + "/wd-pair-3d"

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
#restoreCycle = findLastRestart(restartBaseName)
restoreCycle = None

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
units = PhysicalConstants(0.01,0.001,1.0e-6)
Pmin = 1e-6
Pmax = 1e35
Tmin = 100.0
eos = HelmholtzEquationOfState(units,Pmin,Pmax,Tmin)
#eos = GammaLawMKS3d(4.0/3.0,1.0)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel3d(BSplineKernel3d(), 100)
WTPi = TableKernel3d(BSplineKernel3d(), 100)
output("WT")
output("WTPi")
kernelExtent = WT.kernelExtent

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
nPerh = 1.51
nodes = makeFluidNodeList("nodes", eos,
                          hmin = hmin,
                          hmax = hmax,
                          nPerh = nPerh,
                          xmin = Vector.one * -1e20,
                          xmax = Vector.one * 1e20)
output("nodes")
output("nodes.hmin")
output("nodes.hmax")
output("nodes.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Generate them nodes.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    generator = AsciiFileNodeGenerator3D(filename = "ic_1.0.sdf.ascii",
                                         materialName = "Default",
                                         nNodePerh = nPerh)
    nodes.numInternalNodes = generator.localNumNodes()
    vel = nodes.velocity()
    eps = nodes.specificThermalEnergy()
    abund = []
    vel1 = VectorField("thisField",nodes)

    for i in xrange(nodes.numInternalNodes):
        nodes.velocity()[i].x = generator.vx[i]
        if(generator.vx[i] != 0):
            print "velx = " + str(generator.vx[i])
            print "stored velx = " + str(nodes.velocity()[i].x)
        vel[i].y = generator.vy[i]
        vel[i].z = generator.vz[i]
        eps[i] = generator.eps[i]
        vel1[i] = vel[i]

    distributeNodes((nodes, generator),)
    for i in xrange(nodes.numInternalNodes):
        nodes.velocity()[i] = vel1[i]

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase3d()
output("db")
output("db.appendNodeList(nodes)")
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
hydro = HydroConstructor(WT,
                         WTPi,
                         q,
                         cfl = cfl,
                         compatibleEnergyEvolution = compatibleEnergy,
                         gradhCorrection = gradhCorrection,
                         XSPH = XSPH,
                         densityUpdate = densityUpdate,
                         HUpdate = HEvolution,
                         epsTensile = epsilonTensile,
                         nTensile = nTensile)
output("hydro")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.gradhCorrection")
output("hydro.XSPH")
output("hydro.densityUpdate")
output("hydro.HEvolution")
output("hydro.epsilonTensile")
output("hydro.nTensile")

packages = [hydro]

#-------------------------------------------------------------------------------
# Gimme gravity.
#-------------------------------------------------------------------------------
gravity = OctTreeGravity(G = units.G,
                             softeningLength = plummerLength,
                             opening = opening,
                             ftimestep = fdt,
                             timeStepChoice = timeStepChoice)

packages.append(gravity)

#-------------------------------------------------------------------------------
# Construct a time integrator.
#-------------------------------------------------------------------------------
integrator = myIntegrator(db)
for p in packages:
    integrator.appendPhysicsPackage(p)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.rigorousBoundaries = rigorousBoundaries
integrator.verbose = verbosedt
output("integrator")
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
                            restartBaseName = restartBaseName,
                            vizTime = dtSample,
                            vizDir = visitDir,
                            vizBaseName = "wd-pair-3d",
                            restoreCycle = restoreCycle)
output("control")


#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
hstats([nodes])

control.advance(goalTime)







