#-------------------------------------------------------------------------------
# A rod of stainless steel undergoing tensile strain.  This is intended as a
# test of the cracking/failure models.
#
# This version imposes a constant velocity gradient across the rod, and uses
# a Grady-Kipp flaw distribution.
#
# See Benz & Asphaug (1994), Icarus, 107, 98
#-------------------------------------------------------------------------------
from SolidSpheral3d import *
from SpheralTestUtilities import *
from findLastRestart import *
from SpheralVisitDump import dumpPhysicsState
from identifyFragments import identifyFragments, fragmentProperties
from math import *
import mpi

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("3-D Tensile rod strength/damage model test")

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
        self.controlNodeFlags = ScalarField("flag nodes",
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
commandLine(# How much of the 2 Pi geometry are we doing?
            phiFactor = 2.0,

            problemLabel = "uniqueLabel",

            length = 3.0,
            radius = 0.5,
            nx = 150,
            ny = 25,
            nPerh = 1.51,

            rho0 = 7.9,

            # Material specific bounds on the mass density.
            etamin = 0.5,
            etamax = 1.5,

            # Parameters for the time dependent strain and cracking.
            DamageModelConstructor = GradyKippTensorDamageOwen,
            v0 = 1.0e4,
            kWeibullFactor = 1.0,
            mWeibullFactor = 1.0,
            randomSeed = 548928513,
            strainType = PseudoPlasticStrain,
            damageMethod = Copy,
            useDamageGradient = True,
            cullToWeakestFlaws = False,
            effectiveFlawAlgorithm = SampledFlaws,

            IntegratorConstructor = CheapSynchronousRK2Integrator,
            Qconstructor = MonaghanGingoldViscosity,
            HydroConstructor = SolidASPHHydro,
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
            vizCycle = 100,
            dt = 1e-10,
            dtMin = 1e-12,
            dtMax = 1e-5,
            dtGrowth = 2.0,
            dumpFrac = 0.005,
            maxSteps = None,
            statsStep = 10,
            smoothIters = 0,
            redistributeStep = None,
            HEvolution = IdealH,
            densityUpdate = IntegrateDensity,
            compatibleEnergyEvolution = True,
            gradhCorrection = False,
            domainIndependent = False,
            dtverbose = False,


            restoreCycle = None,
            restartStep = 1000,

            plotFlaws = False,

            dataDirBase = "dumps-TensileRod-3d",
            )

phi = pi * phiFactor

#kWeibull = 8.8e4 * kWeibullFactor
kWeibull = 6.52e5 * kWeibullFactor
mWeibull = 2.63  * mWeibullFactor

dataDir = os.path.join(dataDirBase,
                       problemLabel,
                       "%ix%i" % (nx, 2*ny),
                       "nPerh=%4.2f" % nPerh,
                       str(DamageModelConstructor).split("'")[1],
                       "XSPH=%s" % XSPH,
                       str(strainType),
                       str(damageMethod),
                       "gradDamage=%s" % useDamageGradient,
                       "k=%4.2f_m=%4.2f" % (kWeibull, mWeibull))
restartDir = os.path.join(dataDir,
                          "restarts",
                          "proc-%04i" % mpi.rank)
vizDir = os.path.join(dataDir, "viz")
restartBaseName = os.path.join(restartDir, "TensileRod-%ix%i" % (nx, ny))

volume = pi*radius**2 * length
dx = length/nx
dy = radius/ny

dtSample = dumpFrac*goalTime

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
    if not os.path.exists(vizDir):
        os.makedirs(vizDir)
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
eos = GruneisenEquationOfStateCGS(rho0,    # reference density  
                                  etamin,  # etamin             
                                  etamax,  # etamax             
                                  0.457e6, # C0                 
                                  1.49,    # S1                 
                                  0.0,     # S2                 
                                  0.0,     # S3                 
                                  1.93,    # gamma0             
                                  0.5,     # b                  
                                  55.350)  # atomic weight
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
strengthModel = SteinbergGuinanStrengthCGS(eos,
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
eosDamaged = GruneisenEquationOfStateCGS(rho0,    # reference density  
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
WT = TableKernel(BSplineKernel(), 1000)
WTPi = TableKernel(BSplineKernel(), 1000)
output("WT")
output("WTPi")

#-------------------------------------------------------------------------------
# Create the NodeLists.
#-------------------------------------------------------------------------------
nodes = makeSolidNodeList("Stainless steel", eos, strengthModel,
                          nPerh = nPerh,
                          hmin = hmin,
                          hmax = hmax,
                          hminratio = hminratio,
                          rhoMin = etamin*rho0,
                          rhoMax = etamax*rho0)
nodeSet = [nodes]

#-------------------------------------------------------------------------------
# Set node properties (positions, masses, H's, etc.)
#-------------------------------------------------------------------------------
eps0 = 0.0
if restoreCycle is None:
    print "Generating node distribution."
    from GenerateNodeDistribution3d import *
    #from ParMETISDistributeNodes import distributeNodes3d
    #from MortonOrderDistributeNodes import distributeNodes3d
    from VoronoiDistributeNodes import distributeNodes3d
    #from PeanoHilbertDistributeNodes import distributeNodes3d
    generator = GenerateNodeDistribution3d(ny,
                                           nx,
                                           0,
                                           rho0,
                                           "cylindrical",
                                           rmin = 0.0,
                                           rmax = radius,
                                           thetamin = 0.0,
                                           thetamax = phi,
                                           zmin = -0.5*length,
                                           zmax = 0.5*length,
                                           nNodePerh = nPerh,
                                           SPH = (HydroConstructor is SolidSPHHydro),
                                           )

    # For consistency with our 2-D case, we spin the coordinates so
    # that the rod is aligned along the x axis rather than z.
    R = Tensor(0, 0, 1,
               0, 1, 0,
              -1, 0, 0)
    for i in xrange(generator.localNumNodes()):
        vi = generator.localPosition(i)
        Hi = generator.localHtensor(i)
        v1 = R.dot(vi)
        Hi.rotationalTransform(R)
        generator.x[i] = v1.x
        generator.y[i] = v1.y
        generator.z[i] = v1.z
        generator.H[i] = Hi
    
    distributeNodes3d((nodes, generator))
    output("mpi.reduce(nodes.numInternalNodes, mpi.MIN)")
    output("mpi.reduce(nodes.numInternalNodes, mpi.MAX)")
    output("mpi.reduce(nodes.numInternalNodes, mpi.SUM)")

    # # Set node specific thermal energies
    # eps0 = eos.specificThermalEnergy(rho0, 300.0)
    # nodes.specificThermalEnergy(ScalarField("tmp", nodes, eps0))

    # Set node velocites.
    for i in xrange(nodes.numInternalNodes):
        nodes.velocity()[i].x = nodes.positions()[i].x/(0.5*length)*v0

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
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
 if nodes.positions()[i].x < -0.5*length + 5*dx]
[x1Nodes.append(i) for i in xrange(nodes.numInternalNodes)
 if nodes.positions()[i].x >  0.5*length - 5*dx]
print "Selected %i constant velocity nodes." % (mpi.allreduce(len(x0Nodes) + len(x1Nodes), mpi.SUM))

# Set the nodes we're going to control to one single velocity at each end.
v0 = mpi.allreduce(min([nodes.velocity()[i].x for i in x0Nodes] + [100.0]), mpi.MIN)
v1 = mpi.allreduce(max([nodes.velocity()[i].x for i in x1Nodes] + [-100.0]), mpi.MAX)
for i in x0Nodes:
    nodes.velocity()[i].x = v0
for i in x1Nodes:
    nodes.velocity()[i].x = v1

xbc0 = ConstantVelocityBoundary(nodes, x0Nodes)
xbc1 = ConstantVelocityBoundary(nodes, x1Nodes)

bcs = [xbc0, xbc1]

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
                         useVelocityMagnitudeForDt = True,
                         compatibleEnergyEvolution = compatibleEnergyEvolution,
                         gradhCorrection = gradhCorrection,
                         densityUpdate = densityUpdate,
                         HUpdate = HEvolution,
                         XSPH = XSPH,
                         epsTensile = epsilonTensile,
                         nTensile = nTensile)
output("hydro")
output("hydro.cfl")
output("hydro.useVelocityMagnitudeForDt")
output("hydro.HEvolution")
output("hydro.sumForMassDensity")
output("hydro.compatibleEnergyEvolution")
output("hydro.gradhCorrection")
output("hydro.kernel()")
output("hydro.PiKernel()")

#-------------------------------------------------------------------------------
# Construct a damage model.
#-------------------------------------------------------------------------------
if DamageModelConstructor is GradyKippTensorDamage:
    nflaw = int((nx*ny*ny)*log(nx*ny*ny))
    damageModel = DamageModelConstructor(nodes,
                                         kWeibull,
                                         mWeibull,
                                         volume,
                                         1.0,
                                         WT,
                                         randomSeed,
                                         strainType,
                                         damageMethod,
                                         useDamageGradient,
                                         0.4,
                                         effectiveFlawAlgorithm,
                                         1,
                                         nflaw,
                                         averageFailure)

elif DamageModelConstructor is GradyKippTensorDamageOwen:
    numFlawsPerNode = 1
    volumeMultiplier = 2.0*pi/phi
    damageModel = DamageModelConstructor(nodes,
                                         kWeibull,
                                         mWeibull,
                                         WT,
                                         randomSeed,
                                         volumeMultiplier,
                                         strainType,
                                         damageMethod,
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
integrator.appendPhysicsPackage(damageModel)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.domainDecompositionIndependent = domainIndependent
integrator.verbose = dtverbose
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.havePhysicsPackage(damageModel)")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.domainDecompositionIndependent")
output("integrator.verbose")

#-------------------------------------------------------------------------------
# Add the boundary conditions.
#-------------------------------------------------------------------------------
for package in integrator.physicsPackages():
    for bc in bcs:
        package.appendBoundary(bc)

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            redistributeStep = redistributeStep,
                            vizBaseName = "TensileRod-3d",
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = dtSample)
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
    vx = ScalarField("x velocity", nodes)
    vy = ScalarField("y velocity", nodes)
    vz = ScalarField("z velocity", nodes)
    for i in xrange(nodes.numInternalNodes):
        vx[i] = nodes.velocity()[i].x
        vy[i] = nodes.velocity()[i].y
        vz[i] = nodes.velocity()[i].z

    tdamage = nodes.damage()
    etdamage = nodes.effectiveDamage()
    tstrain = damageModel.effectiveStrain()
    sdamage = ScalarField("damage magnitude", nodes)
    sdamagemin = ScalarField("damage magnitude min", nodes)
    sdamagemax = ScalarField("damage magnitude max", nodes)
    edamage = ScalarField("effective damage magnitude", nodes)
    edamagemin = ScalarField("effective damage magnitude min", nodes)
    edamagemax = ScalarField("effective damage magnitude max", nodes)
    sstrain = ScalarField("strain average", nodes)
    sstrainmin = ScalarField("strain min", nodes)
    sstrainmax = ScalarField("strain max", nodes)
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
    #ufragIndex = identifyFragments(nodes, 1.0, tdamage, 1.0, False)
    dumpPhysicsState(integrator,
                     "TensileRod-3d",
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
                               sstrainmax]
                               #ufragIndex]
                     )

#-------------------------------------------------------------------------------
# Smooth the initial conditions/restore state.
#-------------------------------------------------------------------------------
if restoreCycle is not None:
    control.loadRestartFile(restoreCycle)
    strainHistory.flushHistory()
hstats(nodeSet)

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
        #viz(nodes, damageModel)

#-------------------------------------------------------------------------------
# Plot the final flaw distribution, if requested.
#-------------------------------------------------------------------------------
if plotFlaws and DamageModelConstructor is StochasticWeibullTensorDamageModel:
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
## if DamageModelConstructor == GradyKippScalarDamage:
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
