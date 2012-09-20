#ATS:if SYS_TYPE.startswith('darwin'):
#ATS:    t10 = test(SELF, "--graphics False --clearDirectories True --domainIndependent True --outputFile 'TensileRod-GradyKipp-1d-1proc-reproducing.txt' --referenceFile 'Reference/TensileRod-GradyKipp-1d-1proc-reproducing-darwin-20110307.txt'", np=1, label="Tensile rod domain independence test SERIAL RUN")
#ATS:    t11 = testif(t10, SELF, "--graphics False --clearDirectories False --domainIndependent True --outputFile 'TensileRod-GradyKipp-1d-4proc-reproducing.txt' --comparisonFile 'TensileRod-GradyKipp-1d-1proc-reproducing.txt' --referenceFile 'Reference/TensileRod-GradyKipp-1d-1proc-reproducing-darwin-20110307.txt'", np=4, label="Tensile rod domain independence test 4 DOMAIN RUN")
#ATS:else:
#ATS:    t10 = test(SELF, "--graphics False --clearDirectories True --domainIndependent True --outputFile 'TensileRod-GradyKipp-1d-1proc-reproducing.txt'", np=1, label="Tensile rod domain independence test SERIAL RUN")
#ATS:    t11 = testif(t10, SELF, "--graphics False --clearDirectories False --domainIndependent True --outputFile 'TensileRod-GradyKipp-1d-4proc-reproducing.txt' --comparisonFile 'TensileRod-GradyKipp-1d-1proc-reproducing.txt'", np=4, label="Tensile rod domain independence test 4 DOMAIN RUN")
#-------------------------------------------------------------------------------
# A rod of stainless steel undergoing tensile strain.  This is intended as a
# test of the cracking/failure models.
#
# See Benz & Asphaug (1994), Icarus, 107, 98
#-------------------------------------------------------------------------------
from SolidSpheral1d import *
from SpheralTestUtilities import *
from findLastRestart import *
from SpheralVisitDump import dumpPhysicsState
from identifyFragments import identifyFragments, fragmentProperties
from math import *
import shutil
import mpi

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("1-D Tensile rod strength/damage model test")

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
commandLine(length = 3.0,
            radius = 0.5,
            nx = 100,
            nPerh = 1.25,

            rho0 = 7.9,

            # Material specific bounds on the mass density.
            etamin = 0.5,
            etamax = 1.5,

            # Parameters for the time dependent strain and cracking.
            DamageModelConstructor = GradyKippTensorDamageOwen, # GradyKippTensorDamage, # GradyKippScalarDamage # GradyKippVectorDamage # WeibullTensorDamage, # 
            volumeMultiplier = (3.0/100.0)**2,
            numFlawsPerNode = 1,
            v0 = 1e4,
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
            HydroConstructor = SolidSPHHydro,
            Cl = 1.0,
            Cq = 1.0,
            Qlimiter = False,
            balsaraCorrection = False,
            epsilon2 = 1e-2,
            negligibleSoundSpeed = 1e-5,
            csMultiplier = 1e-4,
            hmin = 1e-5,
            hmax = 20.0,
            cfl = 0.5,
            useVelocityMagnitudeForDt = False,
            XSPH = True,
            epsilonTensile = 0.3,
            nTensile = 4,
            hybridMassDensityThreshold = 0.01,

            goalTime = 50.0e-6,
            steps = None,
            dt = 1e-10,
            dtMin = 1e-12,
            dtMax = 1e-5,
            dtGrowth = 2.0,
            dumpFrac = 0.005,
            maxSteps = None,
            statsStep = 10,
            smoothIters = 0,
            HEvolution = IdealH,
            densityUpdate = IntegrateDensity,
            compatibleEnergyEvolution = True,
            gradhCorrection = False,
            domainIndependent = False,

            restoreCycle = None,
            restartStep = 1000,

            graphics = True,

            testtol = 1.0e-3,
            clearDirectories = False,
            referenceFile = "Reference/TensileRod-GradyKipp-1d-1proc-reproducing-20120920.txt",
            dataDirBase = "dumps-TensileRod-1d",
            outputFile = "None",
            comparisonFile = "None",
            )

#kWeibull = 8.8e4 * kWeibullFactor
#kWeibull = 6.52e3 * kWeibullFactor
kWeibull = 6.52e5 * kWeibullFactor
mWeibull = 2.63   * mWeibullFactor

dataDir = os.path.join(dataDirBase,
                       "nx=%i" % nx,
                       str(DamageModelConstructor).split("'")[1],
                       "k=%4.2f_m=%4.2f" % (kWeibull, mWeibull))
restartDir = os.path.join(dataDir, "restarts")
visitDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, "TensileRod-%i" % nx)

xmin = -0.5*length
xmax =  0.5*length

volume = pi*radius**2 * length
dx = length/nx

dtSample = dumpFrac*goalTime

#-------------------------------------------------------------------------------
# Sampling function to measure the average strain in the volume of the rod.
#-------------------------------------------------------------------------------
class AverageStrain(RestartableObject):
    def __init__(self, damageModel, filename):
        RestartableObject.__init__(self)
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
    if not os.path.exists(visitDir):
        os.makedirs(visitDir)
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
                          rhoMin = etamin*rho0,
                          rhoMax = etamax*rho0)
nodeSet = [nodes]

#-------------------------------------------------------------------------------
# Set node properties (positions, masses, H's, etc.)
#-------------------------------------------------------------------------------
eps0 = 0.0
if restoreCycle is None:
    print "Generating node distribution."
    from DistributeNodes import distributeNodesInRange1d
    distributeNodesInRange1d([(nodes, nx, rho0, (xmin, xmax))])
    output("mpi.reduce(nodes.numInternalNodes, mpi.MIN)")
    output("mpi.reduce(nodes.numInternalNodes, mpi.MAX)")
    output("mpi.reduce(nodes.numInternalNodes, mpi.SUM)")

    # Set node specific thermal energies
    eps0 = eos.specificThermalEnergy(rho0, 300.0)
    nodes.specificThermalEnergy(ScalarField("tmp", nodes, eps0))

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
                                         flawAlgorithm = effectiveFlawAlgorithm)

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
strainHistory = AverageStrain(damageModel,
                              dataDir + "/strainhistory.txt")
control.appendPeriodicWork(strainHistory.sample, 1)

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

#-------------------------------------------------------------------------------
# Plot the state.
#-------------------------------------------------------------------------------
if graphics:
    from SpheralGnuPlotUtilities import *
    state = State(db, integrator.physicsPackages())
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

#-------------------------------------------------------------------------------
# If requested, write out the state in a global ordering to a file.
#-------------------------------------------------------------------------------
if outputFile != "None":
    from SpheralGnuPlotUtilities import multiSort
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
    mof = mortonOrderIndicies(db)
    mo = mpi.reduce(internalValues(mof), mpi.SUM)
    if mpi.rank == 0:
        multiSort(mo, xprof, rhoprof, Pprof, vprof, epsprof, hprof, sprof, dprof)
        f = open(outputFile, "w")
        for (xi, rhoi, Pi, vi, epsi, hi, si, di, mi) in zip(xprof, rhoprof, Pprof, vprof, epsprof, hprof, sprof, dprof, mo):
            f.write((8*"%16.12e " + 9*"%i " + "\n") %
                    (xi, rhoi, Pi, vi, epsi, hi, si, di, mi,
                     unpackElementUL(packElementDouble(xi)),
                     unpackElementUL(packElementDouble(rhoi)),
                     unpackElementUL(packElementDouble(Pi)),
                     unpackElementUL(packElementDouble(vi)),
                     unpackElementUL(packElementDouble(epsi)),
                     unpackElementUL(packElementDouble(hi)),
                     unpackElementUL(packElementDouble(si)),
                     unpackElementUL(packElementDouble(di))))
        f.close()

        #---------------------------------------------------------------------------
        # Check the floating values for the state against reference data.
        #---------------------------------------------------------------------------
        f = open(referenceFile, "r")
        reflines = f.readlines()
        f.close()
        xref =   [float(line.split()[0]) for line in reflines]
        rhoref = [float(line.split()[1]) for line in reflines]
        Pref =   [float(line.split()[2]) for line in reflines]
        vref =   [float(line.split()[3]) for line in reflines]
        epsref = [float(line.split()[4]) for line in reflines]
        href =   [float(line.split()[5]) for line in reflines]
        sref =   [float(line.split()[6]) for line in reflines]
        dref =   [float(line.split()[7]) for line in reflines]

        for f, fref, name, tt in ((xprof, xref, "position", testtol),
                                  (rhoprof, rhoref, "density", testtol),
                                  (Pprof, Pref, "pressure", testtol),
                                  (vprof, vref, "velocity", testtol),
                                  (epsprof, epsref, "specific thermal energy", testtol),
                                  (hprof, href, "h", testtol),
                                  (sprof, sref, "strain", testtol),
                                  (dprof, dref, "damage", testtol)):
            assert len(f) == len(fref)
            for fi, frefi in zip(f, fref):
                if not fuzzyEqual(fi, frefi, tt):
                    raise ValueError, "Comparison to reference %s failed : %g != %g : %g" % (name, fi, frefi, tt)
        print "Floating point comparison test passed."

        #---------------------------------------------------------------------------
        # Also we can optionally compare the current results with another file for
        # bit level consistency.
        #---------------------------------------------------------------------------
        if comparisonFile != "None":
            comparisonFile = os.path.join(dataDir, comparisonFile)
            import filecmp
            assert filecmp.cmp(outputFile, comparisonFile)
