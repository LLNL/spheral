#ATS:test(SELF, "--CRKSPH=True --nx1=256 --ny1=512 --goalTime=4 --cfl=0.25 --Cl=1.0 --Cq=1.0 --clearDirectories=False --filter=0 --nPerh=2.01  --serialDump=True --compatibleEnergy=True", label="RT crk, 256x512", np=20)
#ATS:test(SELF, "--CRKSPH=True --nx1=128 --ny1=256 --goalTime=4 --cfl=0.25 --Cl=1.0 --Cq=1.0 --clearDirectories=False --filter=0 --nPerh=2.01  --serialDump=True --compatibleEnergy=True", label="RT crk, 128x256", np=20)
#ATS:test(SELF, "--CRKSPH=True --nx1=512 --ny1=1028 --goalTime=4 --cfl=0.25 --Cl=1.0 --Cq=1.0 --clearDirectories=False --filter=0 --nPerh=2.01  --serialDump=True --compatibleEnergy=True", label="RT crk, 512x1028", np=20)

#-------------------------------------------------------------------------------
# This is the basic Rayleigh-Taylor Problem
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

class ExponentialDensity:
    def __init__(self,
                 rho1,
                 rho2,
                 delta):
        self.rho1 = rho1
        self.rho2 = rho2
        self.delta = delta
        return
    def __call__(self, r):
        return self.rho1+(self.rho2-self.rho1)/(1+exp(-(r.y-0.5)/delta))

title("Rayleigh-Taylor test problem in 2D")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(nx1 = 128,
            ny1 = 256,
            rhoT = 2.0,
            rhoB = 1.0,
            x0 = 0.0,
            x1 = 0.5,
            y0 = 0.0,
            y1 = 1.0,
            gval = -0.5,
	    w0  = 0.025,
            delta = 0.025, 
            gamma = 1.4,
            mu = 1.0,
            
            nPerh = 2.01,
            
            SVPH = False,
            CRKSPH = False,
            ASPH = False,
            SPH = True,   # This just chooses the H algorithm -- you can use this with CRKSPH for instance.
            filter = 0.0,   # CRKSPH filtering
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
            Cq = 1.0,
            linearInExpansion = False,
            Qlimiter = False,
            balsaraCorrection = False,
            epsilon2 = 1e-2,
            hmin = 0.0001,
            hmax = 0.5,
            hminratio = 0.1,
            cfl = 0.5,
            useVelocityMagnitudeForDt = False,
            XSPH = False,
            epsilonTensile = 0.0,
            nTensile = 8,
            
            IntegratorConstructor = CheapSynchronousRK2Integrator,
            goalTime = 20.0,
            steps = None,
            vizCycle = None,
            vizTime = 0.01,
            dt = 0.0001,
            dtMin = 1.0e-8,
            dtMax = 0.1,
            dtGrowth = 2.0,
            maxSteps = None,
            statsStep = 10,
            smoothIters = 0,
            HUpdate = IdealH,
            domainIndependent = False,
            rigorousBoundaries = False,
            dtverbose = False,
            
            densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
            compatibleEnergy = True,            # <--- Important!  rigorousBoundaries does not work with the compatibleEnergy algorithm currently.
            gradhCorrection = False,
            
            clearDirectories = False,
            restoreCycle = None,
            restartStep = 100,
            redistributeStep = 500,
            checkRestart = False,
            dataDir = "dumps-Rayleigh-Taylor-2d_hopkins",
            outputFile = "None",
            comparisonFile = "None",
            
            serialDump = False, #whether to dump a serial ascii file at the end for viz
            
            bArtificialConduction = False,
            arCondAlpha = 0.5,
            )

# Decide on our hydro algorithm.
if SVPH:
    if ASPH:
        HydroConstructor = ASVPHFacetedHydro
    else:
        HydroConstructor = SVPHFacetedHydro
elif CRKSPH:
    Qconstructor = CRKSPHMonaghanGingoldViscosity
    if ASPH:
        HydroConstructor = ACRKSPHHydro
    else:
        HydroConstructor = CRKSPHHydro
else:
    if ASPH:
        HydroConstructor = ASPHHydro
    else:
        HydroConstructor = SPHHydro

dataDir = os.path.join(dataDir,
                       "gval=%g" % (gval),
                       str(HydroConstructor).split("'")[1].split(".")[-1],
                       "densityUpdate=%s" % (densityUpdate),
                       "XSPH=%s" % XSPH,
                       "filter=%s" % filter,
                       "compatible=%s" % compatibleEnergy,
                       "%s-Cl=%g-Cq=%g" % (str(Qconstructor).split("'")[1].split(".")[-1], Cl, Cq),
                       "%ix%i" % (nx1, ny1),
                       "nPerh=%g-Qhmult=%g" % (nPerh, Qhmult))
restartDir = os.path.join(dataDir, "restarts")
vizDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, "Rayleigh-Taylor-2d")
vizBaseName = "Rayleigh-Taylor-2d"

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
nodes = makeFluidNodeList("High density gas", eos,
                           hmin = hmin,
                           hmax = hmax,
                           hminratio = hminratio,
                           nPerh = nPerh)
output("nodes.name")
output("nodes.hmin")
output("nodes.hmax")
output("nodes.hminratio")
output("nodes.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    # Add some points above and below the problem to represent the infinite atmosphere.
    nybound = 10
    dy = (y1 - y0)/ny1
    generator = GenerateNodeDistribution2d(nx1, ny1 + 2*nybound,
                                           rho = ExponentialDensity(rhoB,
                                                                    rhoT,
                                                                    delta),
                                           distributionType = "xstaggeredLattice",
                                           xmin = (x0,y0 - nybound*dy),
                                           xmax = (x1,y1 + nybound*dy),
                                           nNodePerh = nPerh,
                                           SPH = SPH)

    if mpi.procs > 1:
        from VoronoiDistributeNodes import distributeNodes2d
    else:
        from DistributeNodes import distributeNodes2d

    distributeNodes2d((nodes, generator))

    #Set IC
    vel = nodes.velocity()
    eps = nodes.specificThermalEnergy()
    pos = nodes.positions()
    rho = nodes.massDensity()
    for i in xrange(nodes.numInternalNodes):
        xi, yi = pos[i]
        P0 = rhoT/gamma
        Pi = P0 + gval*rho[i]*(yi-0.5)
        velx = 0.0
        vely = 0.0
        if yi > 0.3 and yi < 0.7:
          vely = w0*(1+cos(8.0*pi*(xi+0.25)))*(1+cos(5.0*pi*(yi-0.5)))
        vel[i]=Vector(velx,vely)
        eps0 = Pi/((gamma - 1.0)*rho[i])
        eps[i]=eps0

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
db.appendNodeList(nodes)
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the artificial viscosity.
#-------------------------------------------------------------------------------
q = Qconstructor(Cl, Cq, linearInExpansion)
q.epsilon2 = epsilon2
q.limiter = Qlimiter
q.balsaraShearCorrection = balsaraCorrection
output("q")
output("q.Cl")
output("q.Cq")
output("q.epsilon2")
output("q.limiter")
output("q.balsaraShearCorrection")
output("q.linearInExpansion")
output("q.quadraticInExpansion")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if SVPH:
    hydro = HydroConstructor(WT, q,
                             cfl = cfl,
                             useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
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
elif CRKSPH:
    hydro = HydroConstructor(WT, WTPi, q,
                             filter = filter,
                             cfl = cfl,
                             useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                             compatibleEnergyEvolution = compatibleEnergy,
                             XSPH = XSPH,
                             densityUpdate = densityUpdate,
                             HUpdate = HUpdate)
else:
    hydro = HydroConstructor(WT,
                             WTPi,
                             q,
                             cfl = cfl,
                             useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
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
# Construct the Artificial Conduction physics object.
#-------------------------------------------------------------------------------
if bArtificialConduction:
    #q.reducingViscosityCorrection = True
    ArtyCond = ArtificialConduction(WT,arCondAlpha)
    
    packages.append(ArtyCond)

#-------------------------------------------------------------------------------
# Construct the gravitational acceleration object.
#-------------------------------------------------------------------------------
pos = nodes.positions()
nodeIndices = vector_of_int()
for i in xrange(nodes.numInternalNodes):
    if pos[i].y > y0 and pos[i].y < y1:
        nodeIndices.append(i)

gravity = ConstantAcceleration2d(Vector2d(0.0, gval),
                                  nodes,
                                  nodeIndices)

packages.append(gravity)

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xp1 = Plane(Vector(x0, y0), Vector( 1.0,  0.0))
xp2 = Plane(Vector(x1, y0), Vector(-1.0,  0.0))
yp1 = Plane(Vector(x0, y0), Vector( 0.0,  1.0))
yp2 = Plane(Vector(x0, y1), Vector( 0.0, -1.0))
xbc = PeriodicBoundary(xp1, xp2)

# The y boundary will be a snapshot of the state of the points above and below
# the y-cutoffs.
pos = nodes.positions()
ylow, yhigh = vector_of_int(), vector_of_int()
for i in xrange(nodes.numInternalNodes):
    if pos[i].y < y0:
        ylow.append(i)
    elif pos[i].y > y1:
        yhigh.append(i)
ybc1 = ConstantBoundary(nodes, ylow, yp1)
print "### BLAGO ####"

ybc2 = ConstantBoundary(nodes, yhigh, yp2)

bcSet = [ybc1, ybc2, xbc]  # <-- ybc should be first!

for bc in bcSet:
    for p in packages:
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
                            initializeDerivatives = True,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            redistributeStep = None,
                            vizBaseName = vizBaseName,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime,
                            vizDerivs = True,
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

if serialDump:
    procs = mpi.procs
    rank = mpi.rank
    serialData = []
    i,j = 0,0
    for i in xrange(procs):
            if rank == i:
                for j in xrange(nodes.numInternalNodes):
                    serialData.append([nodes.positions()[j],3.0/(nodes.Hfield()[j].Trace()),nodes.mass()[j],nodes.massDensity()[j],nodes.specificThermalEnergy()[j]])
    serialData = mpi.reduce(serialData,mpi.SUM)
    if rank == 0:
        f = open(dataDir + "/serialDump.ascii",'w')
        for i in xrange(len(serialData)):
            f.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(i,serialData[i][0][0],serialData[i][0][1],0.0,serialData[i][1],serialData[i][2],serialData[i][3],serialData[i][4]))
        f.close()
