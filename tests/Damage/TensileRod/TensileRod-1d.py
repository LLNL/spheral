#ATS:t10 = test(SELF, "--graphics False --clearDirectories True --domainIndependent True --outputFile 'TensileRod-GradyKipp-1d-1proc-reproducing.txt'", np=1, label="Tensile rod domain independence test SERIAL RUN")
#ATS:t11 = testif(t10, SELF, "--graphics False --clearDirectories False --domainIndependent True --outputFile 'TensileRod-GradyKipp-1d-4proc-reproducing.txt' --comparisonFile 'TensileRod-GradyKipp-1d-1proc-reproducing.txt'", np=4, label="Tensile rod domain independence test 4 DOMAIN RUN")
#ATS:t12 = testif(t10, SELF, "--graphics False --clearDirectories False --domainIndependent True --outputFile 'TensileRod-GradyKipp-1d-1proc-reproducing-restart.txt' --comparisonFile 'TensileRod-GradyKipp-1d-1proc-reproducing.txt' --restoreCycle 500", np=1, label="Tensile rod domain independence test SERIAL RESTART RUN")
#ATS:t13 = testif(t11, SELF, "--graphics False --clearDirectories False --domainIndependent True --outputFile 'TensileRod-GradyKipp-1d-4proc-reproducing-restart.txt' --comparisonFile 'TensileRod-GradyKipp-1d-1proc-reproducing.txt' --restoreCycle 500", np=4, label="Tensile rod domain independence test 4 DOMAIN RESTART RUN")
#-------------------------------------------------------------------------------
# A rod of stainless steel undergoing tensile strain.  This is intended as a
# test of the cracking/failure models.
#
# See Benz & Asphaug (1994), Icarus, 107, 98
#-------------------------------------------------------------------------------
from SolidSpheral1d import *
from SpheralTestUtilities import *
from identifyFragments import identifyFragments, fragmentProperties
from math import *
import os, shutil
import mpi

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("1-D Tensile rod strength/damage model test")

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
commandLine(length = 3.0,
            radius = 0.5,
            nx = 100,

            rho0 = 7.9,

            # Material specific bounds on the mass density.
            etamin = 0.5,
            etamax = 1.5,

            # Parameters for the time dependent strain and cracking.
            DamageModelConstructor = GradyKippTensorDamageOwen,
            volumeMultiplier = (3.0/100.0)**2,
            numFlawsPerNode = 1,
            v0 = 1e-2,
            kWeibullFactor = 1.0,
            mWeibullFactor = 1.0,
            randomSeed = 548928513,
            strainType = PseudoPlasticStrain,
            damageMethod = CopyDamage,
            useDamageGradient = True,
            cullToWeakestFlaws = False,
            effectiveFlawAlgorithm = SampledFlaws,
            damageInCompression = False,

            # Johnson-Cook choices
            D1 = 0.0,
            D2 = 2.0,
            D3 = -1.5,
            D4 = 0.0,
            D5 = 0.0,
            aD1 = 0.0,
            bD1 = 0.0,
            eps0D1 = 0.0,
            aD2 = 0.065,
            bD2 = 2.0,
            eps0D2 = 0.165,
            epsilondot0 = 0.01,
            Tcrit = -2.0,
            sigmamax = -3.0,
            efailmin = 0.1,

            # Optionally we can initialize a break near the origin.
            initialBreakRadius = 0.0,
            
            crksph = False,
            hmin = 1e-5,
            hmax = 1.0,
            cfl = 0.5,
            useVelocityMagnitudeForDt = False,
            XSPH = False,
            epsilonTensile = 0.3,
            nTensile = 4,
            hybridMassDensityThreshold = 0.01,
            filter = 0.0,
            volumeType = CRKSumVolume,

            IntegratorConstructor = CheapSynchronousRK2Integrator,
            goalTime = 50.0,
            steps = None,
            dt = 1e-10,
            dtMin = 1e-6,
            dtMax = 10.0,
            dtGrowth = 2.0,
            dumpFrac = 0.005,
            maxSteps = None,
            statsStep = 10,
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

            graphics = True,

            testtol = 1.0e-4,
            clearDirectories = False,
            referenceFile = "Reference/TensileRod-GradyKipp-1d-1proc-reproducing-20180530.txt",
            dataDirBase = "dumps-TensileRod-1d",
            outputFile = "None",
            comparisonFile = "None",
            )

if crksph:
    hydroname = "CRKSPH"
    nPerh = 1.51
    order = 5
else:
    hydroname = "SPH"
    nPerh = 1.51
    order = 5

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

xmin = -0.5*length
xmax =  0.5*length

volume = pi*radius**2 * length
dx = length/nx

dtSample = dumpFrac*goalTime

#-------------------------------------------------------------------------------
# Sampling function to measure the average strain in the volume of the rod.
#-------------------------------------------------------------------------------
class AverageStrain:
    def __init__(self, damageModel, filename):
        self.restart = RestartableObject(self)
        self.damageModel = damageModel
        self.filename = filename
        self.timeHistory = []
        self.strainHistory = []

        self.file = None
        if mpi.rank == 0:
            self.file = open(self.filename, "w")
            assert self.file is not None

        return

    def sample(self, cycle, atime, dt):
        nodes = self.damageModel.nodeList()
        mass = nodes.mass()
        strain = self.damageModel.strain()

        n = nodes.numInternalNodes
        result = (mpi.allreduce(sum([mass[i]*(strain[i].Trace()) for i in xrange(n)]), mpi.SUM)/
                  mpi.allreduce(sum(mass.internalValues()), mpi.SUM))

        self.timeHistory.append(atime)
        self.strainHistory.append(result)

        if mpi.rank == 0:
            self.file.write("%g   %g\n" % (atime, result))
            self.file.flush()

##         # Check for out of control nodes!
##         vel = [v.x for v in nodes.velocity().internalValues()]
##         for i in xrange(nodes.numInternalNodes):
##             if abs(vel[i]) > 1e6:
##                 print "Bad node: ", i, " vel = ", vel[i]

        return

    def flushHistory(self):
        if mpi.rank == 0:
            assert self.file is not None
            n = len(self.timeHistory)
            assert len(self.strainHistory) == n
            if mpi.rank == 0:
                for i in xrange(n):
                    self.file.write("%g   %g\n" % (self.timeHistory[i],
                                                   self.strainHistory[i]))
            self.file.flush()

    def label(self):
        return "AverageStrain"

    def dumpState(self, file, path):
        file.writeObject(self.timeHistory, path + "/timeHistory")
        file.writeObject(self.strainHistory, path + "/strainHistory")
        return

    def restoreState(self, file, path):
        self.timeHistory = file.readObject(path + "/timeHistory")
        self.strainHistory = file.readObject(path + "/strainHistory")
        return

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if clearDirectories and os.path.exists(dataDir):
        shutil.rmtree(dataDir)
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
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
nodes = makeSolidNodeList("Stainless steel", eos, strengthModel,
                          nPerh = nPerh,
                          hmin = hmin,
                          hmax = hmax,
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
from DistributeNodes import distributeNodesInRange1d
distributeNodesInRange1d([(nodes, nx, rho0, (xmin, xmax))])
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
    fragIDs = nodes.fragmentIDs()
    for i in xrange(nodes.numInternalNodes):
        if abs(pos[i].x) < initialBreakRadius:
            D[i] = SymTensor.one
        if pos[i].x < 0.0:
            fragIDs[i] = 1
        else:
            fragIDs[i] = 2

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
dummy = [x0Nodes.append(i) for i in xrange(nodes.numInternalNodes)
         if nodes.positions()[i].x < -0.5*length + 5*dx]
dummy = [x1Nodes.append(i) for i in xrange(nodes.numInternalNodes)
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
                   HUpdate = HUpdate)
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
                nTensile = nTensile)
output("hydro")
output("hydro.cfl")
output("hydro.useVelocityMagnitudeForDt")
output("hydro.HEvolution")
output("hydro.densityUpdate")
output("hydro.compatibleEnergyEvolution")
output("hydro.kernel()")
output("hydro.PiKernel()")

#-------------------------------------------------------------------------------
# Construct a damage model.
#-------------------------------------------------------------------------------
if DamageModelConstructor is GradyKippTensorDamage:
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
                                         flawAlgorithm = effectiveFlawAlgorithm,
                                         damageInCompression = damageInCompression)

elif DamageModelConstructor is GradyKippTensorDamageOwen:
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
                                         numFlawsPerNode,
                                         damageInCompression = damageInCompression)

elif DamageModelConstructor is JohnsonCookDamage:
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
if DamageModelConstructor in (GradyKippTensorDamageBenzAsphaug, GradyKippTensorDamageOwen):
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
                            restartBaseName = restartBaseName)
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
if DamageModelConstructor in (GradyKippTensorDamageBenzAsphaug, GradyKippTensorDamageOwen):
    strainHistory = AverageStrain(damageModel,
                                  os.path.join(dataDir, "strainhistory.txt"))
    control.appendPeriodicWork(strainHistory.sample, 1)

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)
else:
    control.advance(goalTime)
    control.dropRestartFile()

#-------------------------------------------------------------------------------
# Plot the state.
#-------------------------------------------------------------------------------
if graphics:
    from SpheralGnuPlotUtilities import *
    state = State(db, integrator.physicsPackages())
    H = state.symTensorFields("H")
    h = db.newFluidScalarFieldList(0.0, "h")
    for i in xrange(nodes.numInternalNodes):
        h[0][i] = 1.0/H[0][i].xx
    rhoPlot = plotFieldList(state.scalarFields("mass density"),
                            plotStyle="linespoints",
                            winTitle="rho @ %g %i" % (control.time(), mpi.procs))
    velPlot = plotFieldList(state.vectorFields("velocity"),
                            yFunction = "%s.x",
                            plotStyle="linespoints",
                            winTitle="vel @ %g %i" % (control.time(), mpi.procs))
    mPlot = plotFieldList(state.scalarFields("mass"),
                          plotStyle="linespoints",
                          winTitle="mass @ %g %i" % (control.time(), mpi.procs))
    PPlot = plotFieldList(state.scalarFields("pressure"),
                          plotStyle="linespoints",
                          winTitle="pressure @ %g %i" % (control.time(), mpi.procs))
    hPlot = plotFieldList(h,
                          plotStyle="linespoints",
                          winTitle="h @ %g %i" % (control.time(), mpi.procs))

    d = state.symTensorFields("tensor damage")
    dPlot = plotFieldList(d,
                          yFunction = "%s.xx",
                          winTitle="damage @ %g %i" % (control.time(), mpi.procs),
                          plotStyle="linespoints")
    ed = state.symTensorFields("effective tensor damage")
    edPlot = plotFieldList(ed,
                           yFunction = "%s.xx",
                           winTitle="Effective damage @ %g %i" % (control.time(), mpi.procs),
                           plotStyle="linespoints")

    if DamageModelConstructor in (GradyKippTensorDamageBenzAsphaug, GradyKippTensorDamageOwen):
        ts = damageModel.strain()
        s = ScalarField("strain", nodes)
        for i in xrange(nodes.numInternalNodes):
            s[i] = ts[i].xx
        sl = ScalarFieldList()
        sl.appendField(s)
        sPlot = plotFieldList(sl, winTitle="strain @ %g %i" % (control.time(), mpi.procs),
                              plotStyle="linespoints")
        eps = damageModel.sumActivationEnergiesPerNode()
        nflaws = damageModel.numFlawsPerNode()
        for i in xrange(nodes.numInternalNodes):
            assert nflaws[i] > 0
            eps[i] /= nflaws[i]
        epsl = ScalarFieldList()
        epsl.appendField(eps)
        epsPlot = plotFieldList(epsl, winTitle="Flaw activation strains",
                                plotStyle="linespoints")
      
        eflawsPlot = plotFieldList(state.scalarFields("effective flaws"),
                                   plotStyle = "linespoints",
                                   winTitle = "Effective Flaws @ %g %i" % (control.time(), mpi.procs))

    elif DamageModelConstructor is JohnsonCookDamage:
        eps = damageModel.failureStrain()
        epsl = ScalarFieldList()
        epsl.appendField(eps)
        epsPlot = plotFieldList(epsl, winTitle="JC failure strains",
                                plotStyle="linespoints")
        D1 = damageModel.D1()
        D1l = ScalarFieldList()
        D1l.appendField(D1)
        D1Plot = plotFieldList(D1l, winTitle="JC D1",
                                plotStyle="linespoints")
        D2 = damageModel.D2()
        D2l = ScalarFieldList()
        D2l.appendField(D2)
        D2Plot = plotFieldList(D2l, winTitle="JC D2",
                                plotStyle="linespoints")

    fragPlot = plotFieldList(state.intFields(SolidFieldNames.fragmentIDs),
                             plotStyle = "linespoints",
                             winTitle = "Fragments @  %g %i" % (control.time(), mpi.procs))

#-------------------------------------------------------------------------------
# If requested, write out the state in a global ordering to a file.
#-------------------------------------------------------------------------------
if outputFile != "None":
    from SpheralTestUtilities import multiSort
    state = State(db, integrator.physicsPackages())
    outputFile = os.path.join(dataDir, outputFile)
    pos = state.vectorFields(HydroFieldNames.position)
    rho = state.scalarFields(HydroFieldNames.massDensity)
    P = state.scalarFields(HydroFieldNames.pressure)
    vel = state.vectorFields(HydroFieldNames.velocity)
    eps = state.scalarFields(HydroFieldNames.specificThermalEnergy)
    Hfield = state.symTensorFields(HydroFieldNames.H)
    S = state.symTensorFields(SolidFieldNames.deviatoricStress)
    D = state.symTensorFields(SolidFieldNames.effectiveTensorDamage)
    xprof = mpi.reduce([x.x for x in internalValues(pos)], mpi.SUM)
    rhoprof = mpi.reduce(internalValues(rho), mpi.SUM)
    Pprof = mpi.reduce(internalValues(P), mpi.SUM)
    vprof = mpi.reduce([v.x for v in internalValues(vel)], mpi.SUM)
    epsprof = mpi.reduce(internalValues(eps), mpi.SUM)
    hprof = mpi.reduce([1.0/sqrt(H.Determinant()) for H in internalValues(Hfield)], mpi.SUM)
    sprof = mpi.reduce([x.xx for x in internalValues(S)], mpi.SUM)
    dprof = mpi.reduce([x.xx for x in internalValues(D)], mpi.SUM)
    mof = mortonOrderIndices(db)
    mo = mpi.reduce(internalValues(mof), mpi.SUM)
    if mpi.rank == 0:
        multiSort(mo, xprof, rhoprof, Pprof, vprof, epsprof, hprof, sprof, dprof)
        f = open(outputFile, "w")
        f.write(("#" + 8*" %16s" + "\n") % ("x", "rho", "P", "v", "eps", "h", "S", "D"))
        for (xi, rhoi, Pi, vi, epsi, hi, si, di) in zip(xprof, rhoprof, Pprof, vprof, epsprof, hprof, sprof, dprof):
            f.write((8*"%16.12e " + "\n") %
                    (xi, rhoi, Pi, vi, epsi, hi, si, di))
        f.close()

        #---------------------------------------------------------------------------
        # Check the floating values for the state against reference data.
        #---------------------------------------------------------------------------
        import filearraycmp as fcomp
        assert fcomp.filearraycmp(outputFile, referenceFile, testtol, testtol)
        print "Floating point comparison test passed."

        #---------------------------------------------------------------------------
        # Also we can optionally compare the current results with another file for
        # bit level consistency.
        #---------------------------------------------------------------------------
        if comparisonFile != "None":
            comparisonFile = os.path.join(dataDir, comparisonFile)
            import filecmp
            assert filecmp.cmp(outputFile, comparisonFile)

