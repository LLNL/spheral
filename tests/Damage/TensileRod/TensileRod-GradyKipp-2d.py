#-------------------------------------------------------------------------------
# A rod of stainless steel undergoing tensile strain.  This is intended as a
# test of the cracking/failure models.
#
# This version imposes a constant velocity gradient across the rod, and uses
# a Grady-Kipp flaw distribution.
#
# See Benz & Asphaug (1994), Icarus, 107, 98
#-------------------------------------------------------------------------------
from SolidSpheral import *
from SpheralTestUtilities import *
from findLastRestart import *
from SpheralVisitDump import dumpPhysicsState
from identifyFragments import identifyFragments, fragmentProperties
from math import *
import ParMETISDistributeNodes
import PeanoHilbertDistributeNodes
import MortonOrderDistributeNodes
import mpi

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("2-D Tensile rod strength/damage model test")

#-------------------------------------------------------------------------------
# Stupid little class to override the mass density evolution of the control
# boundary nodes.
#-------------------------------------------------------------------------------
class OverrideNodeProperties(RestartableObject):
    def __init__(self,
                 nodeList,
                 rho0,
                 eps0,
                 controlNodeIDs):
        RestartableObject.__init__(self)
        self.nodeList = nodeList
        self.rho0 = rho0
        self.eps0 = eps0
        self.controlNodeFlags = ScalarField2d("flag nodes",
                                              nodeList,
                                              0.0)
        for i in controlNodeIDs:
            assert i < nodeList.numInternalNodes
            self.controlNodeFlags[i] = 1.0
        return

    def override(self, cycle, t, dt):
        ids = [i for i in xrange(self.nodeList.numInternalNodes)
               if self.controlNodeFlags[i] > 0.5]
        for i in ids:
            self.nodeList.massDensity()[i] = self.rho0
            self.nodeList.specificThermalEnergy()[i] = self.eps0
        return

    def label(self):
        return "OverrideNodeProperties"

    def dumpState(self, file, path):
        file.writeObject(self.rho0, path + "/rho0")
        file.writeObject(self.eps0, path + "/eps0")
        file.write(self.controlNodeFlags, path + "/controlNodeFlags")
        return

    def restoreState(self, file, path):
        self.rho0 = file.readObject(path + "/rho0")
        self.eps0 = file.readObject(path + "/eps0")
        file.read(self.controlNodeFlags, path + "/controlNodeFlags")
        return

#-------------------------------------------------------------------------------
# Generic problem parameters
# All CGS units.
#-------------------------------------------------------------------------------
commandLine(SolidNodeListConstructor = AsphSolidNodeList2d,
            FluidNodeListConstructor = AsphNodeList2d,
            distributor = PeanoHilbertDistributeNodes.distributeNodes2d,

            problemLabel = "uniqueLabel",

            seed = "lattice",

            xlength = 3.0,
            ylength = 1.0,
            nx = 150,
            ny = 50,
            nPerh = 2.01,

            rho0 = 7.9,

            # Material specific bounds on the mass density.
            etamin = 0.8,
            etamax = 1.2,

            # Parameters for the time dependent strain and cracking.
            DamageModelConstructor = GradyKippTensorDamageOwen2d, # StochasticWeibullTensorDamageModel2d, # 
            v0 = 1.0e4,
            kWeibullFactor = 1.0,
            mWeibullFactor = 1.0,
            randomSeed = 548928513,
            strainType = PseudoPlasticStrain,
            effectiveDamage = Copy,
            averageFailure = -1.0,
            useDamageGradient = True,
            cullToWeakestFlaws = False,
            effectiveFlawAlgorithm = FullSpectrumFlaws,

            IntegratorConstructor = CheapSynchronousRK2Integrator2d,
            Qconstructor = MonaghanGingoldViscosity2d,
            Cl = 1.0,
            Cq = 1.0,
            Qlimiter = False,
            balsaraCorrection = False,
            epsilon2 = 1e-2,
            negligibleSoundSpeed = 1e-5,
            csMultiplier = 1e-4,
            hmin = 1e-5,
            hmax = 0.1,
            hminratio = 0.05,
            cfl = 0.5,
            useVelocityMagnitudeForDt = False,
            XSPH = True,
            epsilonTensile = 0.3,
            nTensile = 4,
            hybridMassDensityThreshold = 0.01,

            goalTime = 200.0e-6,
            steps = None,
            dt = 1e-10,
            dtMin = 1e-12,
            dtMax = 1e-5,
            dtGrowth = 2.0,
            dumpFrac = 0.005,
            maxSteps = None,
            statsStep = 10,
            redistributeStep = None,
            smoothIters = 0,
            HEvolution = IdealH,
            sumForMassDensity = RigorousSumDensity, # IntegrateDensity, # HybridDensity # CorrectedSumDensity
            compatibleEnergyEvolution = True,
            gradhCorrection = False,
            dtverbose = False,

            baseDir = ".",

            restoreCycle = None,
            restartStep = 500,

            plotFlaws = False,
            )

kWeibull = 8.8e4 * kWeibullFactor
#kWeibull = 6.52e5 * kWeibullFactor
mWeibull = 2.63  * mWeibullFactor

volumeMultiplier = ylength/ny
numFlawsPerNode = 1 # ny

dataDir = ("%s/dumps-TensileRod-2d/%s/%ix%i/nPerh=%4.2f/%s/%s/XSPH=%s/%s/%s/gradDamage=%s/k=%4.2f_m=%4.2f" %
           (baseDir,
            problemLabel,
            nx, ny,
            nPerh,
            str(DamageModelConstructor).split("'")[1],
            str(SolidNodeListConstructor).split(" ")[1],
            XSPH,
            str(strainType),
            str(effectiveDamage),
            useDamageGradient,
            kWeibull, mWeibull))
restartDir = dataDir + "/restarts/proc-%04i" % mpi.rank
visitDir = dataDir + "/visit"
restartBaseName = restartDir + "/TensileRod-%ix%i" % (nx, ny)

xmin = (-0.5*xlength, -0.5*ylength)
xmax = ( 0.5*xlength,  0.5*ylength)

thickness = ylength
volume = xlength*ylength*thickness
dx = xlength/nx
dy = ylength/ny

origin = Vector2d(-xlength, -ylength)

dtSample = dumpFrac*goalTime

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
    if not os.path.exists(visitDir):
        os.makedirs(visitDir)
if not os.path.exists(restartDir):
    os.makedirs(restartDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# If we're restarting, find the set of most recent restart files.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    restoreCycle = findLastRestart(restartBaseName)

#-------------------------------------------------------------------------------
# Stainless steel material properties.
#-------------------------------------------------------------------------------
eos = GruneisenEquationOfStateCGS2d(rho0,    # reference density  
                                    etamin,  # etamin             
                                    etamax,  # etamax             
                                    0.457e6, # C0                 
                                    1.49,    # S1                 
                                    0.0,     # S2                 
                                    0.0,     # S3                 
                                    1.93,    # gamma0             
                                    0.5,     # b                  
                                    55.350,  # atomic weight
                                    0.0,     # external pressure
                                   -1.0e50,  # minimum pressure
                                    1.0e11)  # maximum pressure
coldFit = NinthOrderPolynomialFit(-1.06797724e10,
                                  -2.06872020e10,
                                   8.24893246e11,
                                  -2.39505843e10,
                                  -2.44522017e10,
                                   5.38030101e10,
                                   0.0,
                                   0.0,
                                   0.0,
                                   0.0)
meltFit = NinthOrderPolynomialFit(7.40464217e10,
                                  2.49802214e11,
                                  1.00445029e12,
                                 -1.36451475e11,
                                  7.72897829e9,
                                  5.06390305e10,
                                  0.0,
                                  0.0,
                                  0.0,
                                  0.0)
strengthModel = SteinbergGuinanStrengthCGS2d(eos,
                                             7.700000e11,        # G0
                                             2.2600e-12,         # A
                                             4.5500e-04,          # B
                                             3.4000e9,           # Y0
                                             2.5e10,             # Ymax
                                             1.0e-3,             # Yp
                                             43.0000,            # beta
                                             0.0,                # gamma0
                                             0.35,               # nhard
                                             coldFit,
                                             meltFit)

# Construct another equation of state for the damaged material.
eosDamaged = GruneisenEquationOfStateCGS2d(rho0,    # reference density  
                                           etamin,  # etamin             
                                           etamax,  # etamax             
                                           0.457e6, # C0                 
                                           1.49,    # S1                 
                                           0.0,     # S2                 
                                           0.0,     # S3                 
                                           1.93,    # gamma0             
                                           0.5,     # b                  
                                           55.350,  # atomic weight
                                           0.0,     # external pressure,
                                           0.0)     # minimum pressure

#-------------------------------------------------------------------------------
# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
#-------------------------------------------------------------------------------
WT = TableKernel2d(BSplineKernel2d(), 1000)
WTPi = TableKernel2d(BSplineKernel2d(), 1000)
output("WT")
output("WTPi")

#-------------------------------------------------------------------------------
# Create the NodeLists.
#-------------------------------------------------------------------------------
nodes = SolidNodeListConstructor("Stainless steel", eos, strengthModel, WT, WTPi)
nodeSet = [nodes]
for n in nodeSet:
    n.nodesPerSmoothingScale = nPerh
    n.epsilonTensile = epsilonTensile
    n.nTensile = nTensile
    n.XSPH = XSPH
    n.hmin = hmin
    n.hmax = hmax
    n.hminratio = hminratio
    n.rhoMin = etamin*rho0
    n.rhoMax = etamax*rho0
    output("n.name")
    output("  n.nodesPerSmoothingScale")
    output("  n.epsilonTensile")
    output("  n.nTensile")
    output("  n.XSPH")
    output("  n.hmin")
    output("  n.hmax")
    output("  n.hminratio")
    output("  n.rhoMin")
    output("  n.rhoMax")
del n

#-------------------------------------------------------------------------------
# Construct the neighbor objects and associate them with the node lists.
#-------------------------------------------------------------------------------
cache = []
for n in nodeSet:
    neighbor = NestedGridNeighbor2d(n,
                                    kernelExtent = WT.kernelExtent)
    n.registerNeighbor(neighbor)
    cache.append(neighbor)
del n

#-------------------------------------------------------------------------------
# Set node properties (positions, masses, H's, etc.)
#-------------------------------------------------------------------------------
eps0 = 0.0
if restoreCycle is None:
    print "Generating node distribution."
    from GenerateNodeDistribution2d import *
    generator = GenerateNodeDistribution2d(nx,
                                           ny,
                                           rho0,
                                           seed,
                                           xmin = xmin,
                                           xmax = xmax,
                                           nNodePerh = nPerh)
    distributor((nodes, generator))
    output("mpi.reduce(nodes.numInternalNodes, mpi.MIN)")
    output("mpi.reduce(nodes.numInternalNodes, mpi.MAX)")
    output("mpi.reduce(nodes.numInternalNodes, mpi.SUM)")

    # Set node specific thermal energies
    eps0 = eos.specificThermalEnergy(rho0, 300.0)
    nodes.specificThermalEnergy(ScalarField2d("tmp", nodes, eps0))

    # Set node velocites.
    for i in xrange(nodes.numInternalNodes):
        nodes.velocity()[i].x = nodes.positions()[i].x/(0.5*xlength)*v0

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase2d()
for n in nodeSet:
    db.appendNodeList(n)
del n
output("db")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct constant velocity boundary conditions to be applied to the rod ends.
#-------------------------------------------------------------------------------
x0Nodes = vector_of_int()
x1Nodes = vector_of_int()
[x0Nodes.append(i) for i in xrange(nodes.numInternalNodes)
 if nodes.positions()[i].x < -0.5*xlength + 5*dx]
[x1Nodes.append(i) for i in xrange(nodes.numInternalNodes)
 if nodes.positions()[i].x >  0.5*xlength - 5*dx]
print "Selected %i constant velocity nodes." % (mpi.allreduce(len(x0Nodes) + len(x1Nodes), mpi.SUM))

# Set the nodes we're going to control to one single velocity at each end.
v0 = mpi.allreduce(min([nodes.velocity()[i].x for i in x0Nodes] + [100.0]), mpi.MIN)
v1 = mpi.allreduce(max([nodes.velocity()[i].x for i in x1Nodes] + [-100.0]), mpi.MAX)
for i in x0Nodes:
    nodes.velocity()[i].x = v0
for i in x1Nodes:
    nodes.velocity()[i].x = v1

xbc0 = ConstantVelocityBoundary2d(nodes, x0Nodes)
xbc1 = ConstantVelocityBoundary2d(nodes, x1Nodes)

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
hydro = Hydro2d(WT,
                WTPi,
                q,
                compatibleEnergyEvolution,
                gradhCorrection,
                sumForMassDensity,
                HEvolution,
                hmin,
                hmax)
hydro.cfl = cfl
hydro.useVelocityMagnitudeForDt = useVelocityMagnitudeForDt
output("hydro")
output("hydro.cfl")
output("hydro.useVelocityMagnitudeForDt")
output("hydro.HEvolution")
output("hydro.sumForMassDensity")
output("hydro.compatibleEnergyEvolution")
output("hydro.gradhCorrection")
output("hydro.hmin")
output("hydro.hmax")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.valid()")

#-------------------------------------------------------------------------------
# Construct a strength physics object.
#-------------------------------------------------------------------------------
strength = Strength2d()
output("strength")

#-------------------------------------------------------------------------------
# Construct a damage model.
#-------------------------------------------------------------------------------
if DamageModelConstructor is GradyKippTensorDamage2d:
    nflaw = int((nx*ny*ny)*log(nx*ny*ny))
    damageModel = DamageModelConstructor(nodes,
                                         kWeibull,
                                         mWeibull,
                                         volume,
                                         1.0,
                                         WT,
                                         randomSeed,
                                         strainType,
                                         effectiveDamage,
                                         useDamageGradient,
                                         0.4,
                                         effectiveFlawAlgorithm,
                                         1,
                                         nflaw,
                                         averageFailure)

elif DamageModelConstructor is GradyKippTensorDamageOwen2d:
    damageModel = DamageModelConstructor(nodes,
                                         kWeibull,
                                         mWeibull,
                                         WT,
                                         randomSeed,
                                         volumeMultiplier,
                                         strainType,
                                         effectiveDamage,
                                         useDamageGradient,
                                         0.4,
                                         effectiveFlawAlgorithm,
                                         numFlawsPerNode)

output("damageModel")
output("damageModel.useDamageGradient")
output("damageModel.effectiveDamageAlgorithm")
output("damageModel.effectiveFlawAlgorithm")

if cullToWeakestFlaws:
    damageModel.cullToWeakestFlaws()

# damageModel.excludeNodes = xNodes
output("damageModel")

#-------------------------------------------------------------------------------
# Construct a time integrator.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
integrator.appendPhysicsPackage(hydro)
integrator.appendPhysicsPackage(strength)
integrator.appendPhysicsPackage(damageModel)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.verbose = dtverbose
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.havePhysicsPackage(strength)")
output("integrator.havePhysicsPackage(damageModel)")
output("integrator.valid()")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.verbose")

#-------------------------------------------------------------------------------
# Add the boundary conditions.
#-------------------------------------------------------------------------------
for package in integrator.physicsPackages():
    for bc in [xbc0, xbc1]:
        package.appendBoundary(bc)

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            redistributeStep = redistributeStep,
                            initializeMassDensity = False)
output("control")

#-------------------------------------------------------------------------------
# Since we are creating a forced velocity gradient on the control nodes, we have
# to override their mass density evolution or they falsely wander away from the
# reference density.
#-------------------------------------------------------------------------------
## override = OverrideNodeProperties(nodes,
##                                   rho0,
##                                   eps0,
##                                   xNodes)
## control.appendPeriodicWork(override.override, 1)

#-------------------------------------------------------------------------------
# Monitor the evolution of the mass averaged strain.
#-------------------------------------------------------------------------------
sys.path.append("../Utilities")
from AverageStrain import AverageStrain
strainHistory = AverageStrain(damageModel,
                              dataDir + "/strainhistory.txt")
control.appendPeriodicWork(strainHistory.sample, 5)

#-------------------------------------------------------------------------------
# Helper method to dump viz files.
#-------------------------------------------------------------------------------
def viz(nodes, damageModel):
    vx = ScalarField2d("x velocity", nodes)
    vy = ScalarField2d("y velocity", nodes)
    for i in xrange(nodes.numInternalNodes):
        vx[i] = nodes.velocity()[i].x
        vy[i] = nodes.velocity()[i].y

    tdamage = nodes.damage()
    etdamage = nodes.effectiveDamage()
    tstrain = damageModel.effectiveStrain()
    sdamage = ScalarField2d("damage magnitude", nodes)
    sdamagemin = ScalarField2d("damage magnitude min", nodes)
    sdamagemax = ScalarField2d("damage magnitude max", nodes)
    edamage = ScalarField2d("effective damage magnitude", nodes)
    edamagemin = ScalarField2d("effective damage magnitude min", nodes)
    edamagemax = ScalarField2d("effective damage magnitude max", nodes)
    sstrain = ScalarField2d("strain average", nodes)
    sstrainmin = ScalarField2d("strain min", nodes)
    sstrainmax = ScalarField2d("strain max", nodes)
    for i in xrange(nodes.numInternalNodes):
        sdamage[i] = tdamage[i].Trace()
        sdamagemin[i] = tdamage[i].eigenValues().minElement()
        sdamagemax[i] = tdamage[i].eigenValues().maxElement()
        edamage[i] = etdamage[i].Trace()
        edamagemin[i] = etdamage[i].eigenValues().minElement()
        edamagemax[i] = etdamage[i].eigenValues().maxElement()
        sstrain[i] = tstrain[i].Trace()/2.0
        sstrainmin[i] = tstrain[i].eigenValues().minElement()
        sstrainmax[i] = tstrain[i].eigenValues().maxElement()
    ufragIndex = identifyFragments(nodes, 2.0, tdamage, 1.0, True)
    dumpPhysicsState(integrator,
                     "TensileRod-2d",
                     visitDir,
                     fields = [vx, vy,
                               damageModel.sumActivationEnergiesPerNode(),
                               damageModel.numFlawsPerNode(),
                               sdamage,
                               sdamagemin,
                               sdamagemax,
                               edamage,
                               edamagemin,
                               edamagemax,
                               sstrain,
                               sstrainmin,
                               sstrainmax,
                               ufragIndex]
                     )

#-------------------------------------------------------------------------------
# Hack!
#-------------------------------------------------------------------------------
## for i in xrange(nodes.numInternalNodes):
##     ri = nodes.positions()[i]
##     if ri.x >= -0.02 and ri.x <= 0.02:
##         nodes.damage()[i] = SymTensor2d.one

#-------------------------------------------------------------------------------
# Smooth the initial conditions/restore state.
#-------------------------------------------------------------------------------
if restoreCycle is not None:
    control.loadRestartFile(restoreCycle)
    strainHistory.flushHistory()

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)
else:
    while control.time() < goalTime:
        nextGoalTime = min(control.time() + dtSample, goalTime)
        control.advance(nextGoalTime, maxSteps)
        control.dropRestartFile()
        viz(nodes, damageModel)

#-------------------------------------------------------------------------------
# Plot the final flaw distribution, if requested.
#-------------------------------------------------------------------------------
if plotFlaws and DamageModelConstructor is StochasticWeibullTensorDamageModel2d:
    from SpheralGnuPlotUtilities import generateNewGnuPlot
    import Gnuplot
    flaws = [x for x in damageModel.sortedFlaws()]
    f = [i/volume for i in xrange(len(flaws))]
    fans = [kWeibull * x**mWeibull for x in flaws]
    d = Gnuplot.Data(flaws, f,
                     title = "Simulation",
                     inline = True)
    dans = Gnuplot.Data(flaws, fans,
                        with_ = "lines",
                        title = "Analytic",
                        inline = True)
    p = generateNewGnuPlot()
    p("set logscale xy")
    p("set key top left")
    p.plot(d)
    p.replot(dans)

#-------------------------------------------------------------------------------
# Compute the fragment properties and spit them out to a file.
#-------------------------------------------------------------------------------
## import pickle
## if DamageModelConstructor == GradyKippScalarDamage2d:
##     ufragIndex, dfragIndex = identifyFragments(nodes,
##                                                1.0,
##                                                damageModel,
##                                                nodesDamaged)
##     fragProps = fragmentProperties(nodes, ufragIndex,
##                                    nodesDamaged, dfragIndex)
## else:
##     ufragIndex = identifyFragments(nodes, 1.0)
##     fragProps = fragmentProperties(nodes, ufragIndex)

## if procID == 0:
##     f = open(dataDir + "/fragment-properties-t=%3.1fmicrosec" % (1.0e6*control.time()), "w")
##     pickle.dump(fragProps, f)
##     f.close()
