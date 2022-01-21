#ATS:test(SELF, "--CRKSPH=False --PSPH=True --KernelConstructor=QuinticSplineKernel --order=5 --goalTime=4 --cfl=0.25 --Cl=1.0 --Cq=1.0 --clearDirectories=False --filter=0 --nPerh=4.01  --serialDump=True --compatibleEnergy=False --evolveTotalEnergy=True --boolCullenViscosity=True --HopkinsConductivity=True", label="comp", np=80)



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
from NodeHistory import NodeHistory
from RTMixLength import RTMixLength

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
            
            nPerh = 1.01,
            
            SVPH = False,
            CRKSPH = False,
            PSPH = False,
            SPH = True,   # This just chooses the H algorithm -- you can use this with CRKSPH for instance.
            filter = 0.0,   # CRKSPH filtering
            Qconstructor = MonaghanGingoldViscosity,
            #Qconstructor = TensorMonaghanGingoldViscosity,
            KernelConstructor = BSplineKernel,
            order = 5, 
            linearConsistent = False,
            fcentroidal = 0.0,
            fcellPressure = 0.0,
            boolReduceViscosity = False,
            nh = 5.0,
            aMin = 0.1,
            aMax = 2.0,
            Qhmult = 1.0,
            boolCullenViscosity = False,
            alphMax = 2.0,
            alphMin = 0.02,
            betaC = 0.7,
            betaD = 0.05,
            betaE = 1.0,
            fKern = 1.0/3.0,
            boolHopkinsCorrection = True,
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
            goalTime = 4.0,
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
           
            correctionOrder = LinearOrder, 
            densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
            compatibleEnergy = True,            # <--- Important!  rigorousBoundaries does not work with the compatibleEnergy algorithm currently.
            gradhCorrection = True,
            correctVelocityGradient = True,
            evolveTotalEnergy= False,
            HopkinsConductivity=False,
            
            clearDirectories = False,
            restoreCycle = -1,
            restartStep = 100,
            redistributeStep = 500,
            sampleFreq = 20,
            dataDir = "dumps-Rayleigh-Taylor-2d_hopkins",
            outputFile = "RT_Hopkins.txt",
            comparisonFile = "None",
            
            serialDump = False, #whether to dump a serial ascii file at the end for viz
            useVoronoiOutput = False,
            
            bArtificialConduction = False,
            arCondAlpha = 0.5,
            )

assert not(boolReduceViscosity and boolCullenViscosity)
# Decide on our hydro algorithm.
if SVPH:
    if SPH:
        HydroConstructor = SVPHFacetedHydro
    else:
        HydroConstructor = ASVPHFacetedHydro
elif CRKSPH:
    Qconstructor = CRKSPHMonaghanGingoldViscosity
    if SPH:
        HydroConstructor = CRKSPHHydro
    else:
        HydroConstructor = ACRKSPHHydro
elif PSPH:
    if SPH:
        HydroConstructor = PSPHHydro
    else:
        HydroConstructor = APSPHHydro
else:
    if SPH:
        HydroConstructor = SPHHydro
    else:
        HydroConstructor = ASPHHydro

dataDir = os.path.join(dataDir,
                       "gval=%g" % (gval),
                       "w0=%g" % w0,
                       HydroConstructor.__name__,
                       Qconstructor.__name__,
                       KernelConstructor.__name__,
                       "densityUpdate=%s" % (densityUpdate),
                       "correctionOrder=%s" % (correctionOrder),
                       "XSPH=%s" % XSPH,
                       "filter=%s" % filter,
                       "compatible=%s" % compatibleEnergy,
                       "Cullen=%s" % boolCullenViscosity,
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
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS(gamma, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
if KernelConstructor==NBSplineKernel:
  WT = TableKernel(NBSplineKernel(order), 10000)
  WTPi = TableKernel(NBSplineKernel(order), 10000, Qhmult)
else:
  WT = TableKernel(KernelConstructor(), 10000)
  WTPi = TableKernel(KernelConstructor(), 10000, Qhmult)
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
                          nPerh = nPerh,
                          kernelExtent = kernelExtent)
output("nodes.name")
output("nodes.hmin")
output("nodes.hmax")
output("nodes.hminratio")
output("nodes.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
# Add some points above and below the problem to represent the infinite atmosphere.
nybound = 20
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
    hydro = HydroConstructor(W = WT,
                             Q = q,
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
    hydro = HydroConstructor(W = WT,
                             WPi = WTPi,
                             Q = q,
                             filter = filter,
                             cfl = cfl,
                             useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                             compatibleEnergyEvolution = compatibleEnergy,
                             XSPH = XSPH,
                             correctionOrder = correctionOrder,
                             densityUpdate = densityUpdate,
                             HUpdate = HUpdate)
elif PSPH:
    hydro = HydroConstructor(W = WT,
                             Q = q,
                             filter = filter,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             evolveTotalEnergy = evolveTotalEnergy,
                             HopkinsConductivity = HopkinsConductivity,
                             densityUpdate = densityUpdate,
                             correctVelocityGradient = correctVelocityGradient,
                             HUpdate = HUpdate,
                             XSPH = XSPH)
else:
    hydro = HydroConstructor(W = WT,
                             WPi = WTPi,
                             Q = q,
                             cfl = cfl,
                             useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                             compatibleEnergyEvolution = compatibleEnergy,
                             gradhCorrection = gradhCorrection,
                             correctVelocityGradient = correctVelocityGradient,
                             evolveTotalEnergy = evolveTotalEnergy,
                             XSPH = XSPH,
                             densityUpdate = densityUpdate,
                             HUpdate = HUpdate,
                             epsTensile = epsilonTensile)
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
elif boolCullenViscosity:
    evolveCullenViscosityMultiplier = CullenDehnenViscosity(q,WTPi,alphMax,alphMin,betaC,betaD,betaE,fKern,boolHopkinsCorrection)
    packages.append(evolveCullenViscosityMultiplier)

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
# Build the history object to measure the time-dependent mix length.
#-------------------------------------------------------------------------------
mixlengthhistory = NodeHistory(nodes, [], 
                               sampleMethod = RTMixLength(nodes, 0.5, 95.0),
                               filename = os.path.join(dataDir, "mix_length.gnu"),
                               labels = ("yhigh", "ylow", "L"))

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
if useVoronoiOutput:
    import SpheralVoronoiSiloDump
    vizMethod = SpheralVoronoiSiloDump.dumpPhysicsState
else:
    import SpheralPointmeshSiloDump
    vizMethod = SpheralPointmeshSiloDump.dumpPhysicsState

control = SpheralController(integrator, WT,
                            initializeDerivatives = True,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            redistributeStep = None,
                            vizMethod = vizMethod,
                            vizBaseName = vizBaseName,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime,
                            SPH = SPH)
output("control")

# Add the periodic work.
#control.appendPeriodicWork(mixlengthhistory.sample, sampleFreq)
#mixlengthhistory.flushHistory()  # <-- in case of restart

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
        f = open(outputFile,'w')
        for i in xrange(len(serialData)):
            f.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(i,serialData[i][0][0],serialData[i][0][1],0.0,serialData[i][1],serialData[i][2],serialData[i][3],serialData[i][4]))
        f.close()
