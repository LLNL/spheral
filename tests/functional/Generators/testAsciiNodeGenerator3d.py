import shutil
from math import *
from Spheral import *
from Spheral3d import *
from GenerateNodeDistribution3d import *
from AsciiFileNodeGenerator import *
from VoronoiDistributeNodes import distributeNodes3d as distributeNodes
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from findLastRestart import *
from SpheralVisitDump import *
import mpi

title("White Dwarf pair test from ic readin file")

class WDUnits(PhysicalConstants):
    def __init__(self):
        PhysicalConstants.__init__(self,
                                   6.3678e6,  #Unit Length Earth Radii (meters)
                                   1.989e30,  #Unit Mass Solar Mass (Kg)
                                   1.0      ) #Unit Time (s)
        return

def hashEnergy(generator):
    e_hash = []
    e_heps = []
    for i in xrange(generator.localNumNodes()):
        e_hash.append( int( 1.0e15*abs(generator.x[i]) + 1.0e10*abs(generator.y[i]) + 1.0e5*abs(generator.z[i]) ) )
        e_heps.append( generator.eps[i] )
    #Pass it around
    e_hash = mpi.allreduce(e_hash,mpi.SUM)
    e_heps = mpi.allreduce(e_heps,mpi.SUM)
    return e_hash,e_heps

def energyHash(nodeList,e_hash,e_heps):
    eps = nodeList.specificThermalEnergy()
    pos = nodeList.positions()
    
    for i in xrange(nodeList.numInternalNodes):
        cval = int( 1.0e15*abs(pos[i].x) + 1.0e10*abs(pos[i].y) + 1.0e5*abs(pos[i].z) )
        indx = -100
        for j in xrange(len(e_hash)):
            if(e_hash[j] == cval):
                indx = j
                break
        if(indx == -100):
            sys.exit('Problem in Hash')
        eps[i] = e_heps[indx]
    return




#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(star1 = "AsciiTest.ascii",
            #Sim params
            useGravity      = True,
            plummerLength   = 0.05,    # (cm) Plummer softening scale
            opening         = 1.0,      # (dimensionless, OctTreeGravity) opening parameter for tree walk
            fdt             = 0.5,      # (dimensionless, OctTreeGravity) timestep multiplier
            timeStepType    = DynamicalTime,
            
            nPerh   = 1.51, #2.01,
            CRKSPH  = False,
            filter  = 0.0,
            
            
            gamma   = 5.0/3.0,
            kappa   = 1.382e15,   #8.63e13,      # this is the constant of proportionality for the sun
            #The above kappa is taken from public poly from FXT for a 5M n=5 polytrope
            mue = 13.71,
            muc = 1.3,
            
            #Artificial Viscosity and things
            momentumConserving  = True,  #for CSPH
            Qconstructor        = MonaghanGingoldViscosity,
            Cl      = 1.0,
            Cq      = 1.0,
            Qhmult  = 1.0,
            Qlimiter = False,
            balsaraCorrection = False,
            epsilon2 = 1e-2,
            epsilonTensile = 0.0,
            nTensile = 8,
            negligibleSoundSpeed = 1e-10,
            csMultiplier = 1e-4,
            hmin = 1e-8,
            hmax = 1e5,
            hminratio = 0.05,
            cfl = 0.5,
            XSPH = True,
            
            #Hydro params
            HEvolution                  = IdealH,
            densityUpdate               = RigorousSumDensity,
            compatibleEnergyEvolution   = False,
            gradhCorrection = False,
            

            
            #Time and sim control
            timeStepChoice = AccelerationRatio,
            myIntegrator = CheapSynchronousRK2Integrator3d,
            steps = None,
            goalTime = 100.,
            vizTime = 10e0,
            vizCycle = 100,
            dt = 1e-10,
            dtMin = 1e-15,
            dtMax = 1e5,
            dtGrowth = 2.0,
            dtSample = 1,
            rigorousBoundaries = False,
            verbosedt = False,
            maxSteps = None,
            statsStep = 10,
            redistributeStep = 100,
            restartStep = 50,
            restoreCycle = None,
            
            clearDirectories = False,
            dataDir = "dumps-ascii-test",
            
            serialDump=False,
            )


restartDir = dataDir + "/restarts"
visitDir = dataDir + "/visit"
restartBaseName = restartDir + "/wd"

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if clearDirectories and os.path.exists(dataDir):
        shutil.rmtree(dataDir)
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
    if not os.path.exists(visitDir):
        os.makedirs(visitDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
units = WDUnits()
Pmin = 1e-6
Pmax = 1e35
Tmin = 100.0
#eos = HelmholtzEquationOfState(units,Pmin,Pmax,Tmin)
eos = GammaLawGas(4.0/3.0,1.0,
                  constants = units)

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
nodes1 = makeFluidNodeList("nodes1", eos,
                           hmin = hmin,
                           hmax = hmax,
                           nPerh = nPerh,
                           xmin = Vector.one * -1e20,
                           xmax = Vector.one * 1e20)
output("nodes1")
output("nodes1.hmin")
output("nodes1.hmax")
output("nodes1.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Generate them nodes.
#-------------------------------------------------------------------------------
generator1 = AsciiFileNodeGenerator3D(filename=star1,
                                      materialName = "Default",
                                      nNodePerh = nPerh,
                                      nodes=nodes1)

#msum = mpi.allreduce(sum(generator1.m + [0.0]), mpi.SUM)
#assert msum > 0.0
#print "Found star1 mass = %g." % (msum)

#ehash1,eheps1 = hashEnergy(generator1)

distributeNodes((nodes1, generator1))

eps = nodes1.specificThermalEnergy()
for i in xrange(nodes1.numInternalNodes):
    eps[i] = generator1.eps[i]

#energyHash(nodes1,ehash1,eheps1)

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase3d()
output("db")
output("db.appendNodeList(nodes1)")
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
hydro = SPH(dataBase=db,
            Q=q,
            W=WT,
            WPi=WTPi,
            cfl = cfl,
            compatibleEnergyEvolution = compatibleEnergyEvolution,
            #gradhCorrection = gradhCorrection,
                         XSPH = XSPH,
            densityUpdate = densityUpdate,
            HUpdate = HEvolution,
            #epsTensile = epsilonTensile,
            #nTensile = nTensile
                         )
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

print "building controller"

control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            vizTime = vizTime,
                            vizDir = visitDir,
                            vizStep = vizCycle,
                            vizBaseName = "ascii",
                            restoreCycle = restoreCycle)
output("control")


#control.appendPeriodicWork(restartBurn,restartStep)


#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------

if steps is None:
    control.advance(goalTime)
else:
    control.step(steps)
