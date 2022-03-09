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
from SpheralVisitDump import dumpPhysicsState
from identifyFragments import identifyFragments, fragmentProperties
from math import *
import os, shutil
import mpi

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("3-D Tensile rod strength/damage model test")

#-------------------------------------------------------------------------------
# Stupid little class to override the mass density evolution of the control
# boundary nodes.
#-------------------------------------------------------------------------------
class OverrideNodeProperties:
    def __init__(self,
                 nodeList,
                 rho0,
                 eps0,
                 controlNodeIDs):
        self.restart = RestartableObject(self)
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
# Build our units.
#-------------------------------------------------------------------------------
units = PhysicalConstants(0.01,  # Unit length in m
                          0.001, # Unit mass in kg
                          1e-6)  # Unit length in sec

#-------------------------------------------------------------------------------
# Generic problem parameters
# All (cm, gm, usec) units.
#-------------------------------------------------------------------------------
commandLine(# How much of the 2 Pi geometry are we doing?
            phiFactor = 2.0,

            length = 3.0,
            radius = 0.5,
            nx = 150,
            ny = 25,

            rho0 = 7.9,

            # Material specific bounds on the mass density.
            etamin = 0.5,
            etamax = 1.5,

            # Parameters for the time dependent strain and cracking.
            DamageModelConstructor = GradyKippTensorDamageOwen,
            numFlawsPerNode = 1,
            v0 = 1e-2,
            kWeibullFactor = 1.0,
            mWeibullFactor = 1.0,
            randomSeed = 548928513,
            strainType = PseudoPlasticStrain,
            damageCoupling = ThreePointDamage,
            cullToWeakestFlaws = False,
            damageInCompression = False,
            negativePressureInDamage = False,

            # Johnson-Cook choices
            D1 = 0.0,
            D2 = 2.0,
            D3 = -1.5,
            D4 = 0.0,
            D5 = 0.0,
            epsilondot0 = 0.01,
            Tcrit = -2.0,
            sigmamax = -3.0,
            efailmin = 0.01,
            sigmaD1 = 0.0,  # Gaussian
            aD1 = 0.0,      # Weibull
            bD1 = 0.0,      # Weibull
            eps0D1 = 0.0,   # Weibull
            sigmaD2 = 0.05, # Gaussian
            aD2 = 0.065,    # Weibull
            bD2 = 2.0,      # Weibull
            eps0D2 = 0.165, # Weibull

            # Optionally we can initialize a break near the origin.
            initialBreakRadius = 0.0,
            
            crksph = False,
            ASPH = True,
            hminratio = 0.05,
            cfl = 0.25,
            useVelocityMagnitudeForDt = False,
            XSPH = False,
            epsilonTensile = 0.0,
            nTensile = 4,
            hybridMassDensityThreshold = 0.01,
            filter = 0.0,
            volumeType = RKSumVolume,

            IntegratorConstructor = CheapSynchronousRK2Integrator,
            goalTime = 200.0,
            steps = None,
            dt = 1e-10,
            dtMin = 1e-6,
            dtMax = 10.0,
            dtGrowth = 2.0,
            dumpFrac = 0.005,
            maxSteps = None,
            statsStep = 10,
            redistributeStep = None,
            vizCycle = None,
            vizTime = 1.0,
            smoothIters = 0,
            HUpdate = IdealH,
            densityUpdate = IntegrateDensity,
            compatibleEnergy = True,
            gradhCorrection = True,
            correctVelocityGradient = True,
            domainIndependent = True,
            dtverbose = False,

            restoreCycle = -1,
            restartStep = 500,

            plotFlaws = False,

            clearDirectories = False,
            dataDirBase = "dumps-TensileRod-3d",
            timerFile = "time.table",
            )

phi = pi * phiFactor

if crksph:
    hydroname = os.path.join("CRKSPH", str(volumeType))
    nPerh = 1.51
    order = 5
else:
    hydroname = "SPH"
    nPerh = 1.51
    order = 5
if ASPH:
    hydroname = "A" + hydroname

#kWeibull = 8.8e4 * kWeibullFactor
#kWeibull = 6.52e3 * kWeibullFactor
kWeibull = 6.52e5 * kWeibullFactor
mWeibull = 2.63   * mWeibullFactor

dataDir = os.path.join(dataDirBase,
                       hydroname,
                       DamageModelConstructor.__name__,
                       "nx=%i" % nx,
                       "k=%4.2f_m=%4.2f" % (kWeibull, mWeibull))
restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "TensileRod-%i" % nx)
vizDir = os.path.join(dataDir, "visit")
vizBaseName = "TensileRod-%i" % nx

volume = pi*radius**2 * length
dx = length/nx
dy = radius/ny

dtSample = dumpFrac*goalTime

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
# Stainless steel material properties.
#-------------------------------------------------------------------------------
eos = GruneisenEquationOfState(rho0,    # reference density  
                               etamin,  # etamin             
                               etamax,  # etamax             
                               0.457, # C0                 
                               1.49,    # S1                 
                               0.0,     # S2                 
                               0.0,     # S3                 
                               1.93,    # gamma0             
                               0.5,     # b                  
                               55.350,  # atomic weight
                               units)
coldFit = NinthOrderPolynomialFit(-1.06797724e-2,
                                  -2.06872020e-2,
                                   8.24893246e-1,
                                  -2.39505843e-2,
                                  -2.44522017e-2,
                                   5.38030101e-2,
                                   0.0,
                                   0.0,
                                   0.0,
                                   0.0)
meltFit = NinthOrderPolynomialFit(7.40464217e-2,
                                  2.49802214e-2,
                                  1.00445029e-2,
                                 -1.36451475e-1,
                                  7.72897829e-3,
                                  5.06390305e-2,
                                  0.0,
                                  0.0,
                                  0.0,
                                  0.0)
strengthModel = SteinbergGuinanStrength(eos,
                                        7.700000e-1,        # G0
                                        2.2600,             # A
                                        4.5500e-04,         # B
                                        3.4000e-3,          # Y0
                                        2.5e-2,             # Ymax
                                        1.0e-3,             # Yp
                                        43.0000,            # beta
                                        0.0,                # gamma0
                                        0.35,               # nhard
                                        coldFit,
                                        meltFit)

#-------------------------------------------------------------------------------
# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
#-------------------------------------------------------------------------------
WT = TableKernel(NBSplineKernel(order), 1000)
output("WT")

#-------------------------------------------------------------------------------
# Create the NodeLists.
#-------------------------------------------------------------------------------
hmin = 0.1*nPerh*min(length/nx, radius/ny)
hmax = 2.0*nPerh*max(length/nx, radius/ny)
nodes = makeSolidNodeList("Stainless steel", eos, strengthModel,
                          nPerh = nPerh,
                          hmin = hmin,
                          hmax = hmax,
                          hminratio = hminratio,
                          rhoMin = etamin*rho0,
                          rhoMax = etamax*rho0,
                          xmin = -100.0*Vector.one,
                          xmax =  100.0*Vector.one)
nodeSet = [nodes]

#-------------------------------------------------------------------------------
# Set node properties (positions, masses, H's, etc.)
#-------------------------------------------------------------------------------
eps0 = 0.0
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
                                       SPH = not ASPH,
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

# Set an initial damage if requested.
if initialBreakRadius > 0.0:
    pos = nodes.positions()
    D = nodes.damage()
    for i in xrange(nodes.numInternalNodes):
        if abs(pos[i].x) < initialBreakRadius:
            D[i] = SymTensor.one

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
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if crksph:
    hydro = CRKSPH(dataBase = db,
                   W = WT,
                   filter = filter,
                   cfl = cfl,
                   compatibleEnergyEvolution = compatibleEnergy,
                   XSPH = XSPH,
                   densityUpdate = densityUpdate,
                   volumeType = volumeType,
                   HUpdate = HUpdate,
                   ASPH = ASPH)
else:
    hydro = SPH(dataBase = db,
                W = WT,
                filter = filter,
                cfl = cfl,
                compatibleEnergyEvolution = compatibleEnergy,
                gradhCorrection = gradhCorrection,
                correctVelocityGradient = correctVelocityGradient,
                densityUpdate = densityUpdate,
                HUpdate = HUpdate,
                XSPH = XSPH,
                epsTensile = epsilonTensile,
                nTensile = nTensile,
                ASPH = ASPH,
                negativePressureInDamage = negativePressureInDamage)
output("hydro")
output("hydro.cfl")
output("hydro.useVelocityMagnitudeForDt")
output("hydro.HEvolution")
output("hydro.densityUpdate")
output("hydro.compatibleEnergyEvolution")
output("hydro.kernel")
output("hydro.PiKernel")
output("hydro.negativePressureInDamage")

#-------------------------------------------------------------------------------
# Construct a damage model.
#-------------------------------------------------------------------------------
if DamageModelConstructor is GradyKippTensorDamage:
    damageModel = DamageModelConstructor(nodeList = nodes,
                                         kWeibull = kWeibull,
                                         mWeibull = mWeibull,
                                         volume = volume,
                                         volumeStretchFactor = 1.0,
                                         kernel = WT,
                                         seed = randomSeed,
                                         strainAlgorithm = strainType,
                                         damageCouplingAlgorithm = damageCoupling,
                                         damageInCompression = damageInCompression)

elif DamageModelConstructor is GradyKippTensorDamageOwen:
    volumeMultiplier = 2.0*pi/phi
    damageModel = DamageModelConstructor(nodeList = nodes,
                                         kWeibull = kWeibull,
                                         mWeibull = mWeibull,
                                         kernel = WT,
                                         seed = randomSeed,
                                         volumeMultiplier = volumeMultiplier,
                                         damageCouplingAlgorithm = damageCoupling,
                                         strainAlgorithm = strainType,
                                         minFlawsPerNode = numFlawsPerNode,
                                         damageInCompression = damageInCompression)

elif DamageModelConstructor is JohnsonCookDamageWeibull:
    damageModel = DamageModelConstructor(nodes,
                                         D1 = D1,
                                         D2 = D2,
                                         D3 = D3,
                                         D4 = D4,
                                         D5 = D5,
                                         aD1 = aD1,
                                         bD1 = bD1,
                                         eps0D1 = eps0D1,
                                         aD2 = aD2,
                                         bD2 = bD2,
                                         eps0D2 = eps0D2,
                                         epsilondot0 = epsilondot0,
                                         Tcrit = Tcrit,
                                         sigmamax = sigmamax,
                                         efailmin = efailmin,
                                         seed = randomSeed,
                                         domainIndependent = domainIndependent)

output("damageModel")

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
                            restoreCycle = restoreCycle,
                            vizBaseName = vizBaseName,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime,
                            timerName = timerFile,
                            SPH = not ASPH)
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
from AverageStrain import AverageStrain
strainHistory = AverageStrain(damageModel,
                              dataDir + "/strainhistory.txt")
control.appendPeriodicWork(strainHistory.sample, 1)

#-------------------------------------------------------------------------------
# Smooth the initial conditions/restore state.
#-------------------------------------------------------------------------------
if restoreCycle is not None:
    strainHistory.flushHistory()
hstats(nodeSet)

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)
else:
    control.advance(goalTime)
    control.dropRestartFile()

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
