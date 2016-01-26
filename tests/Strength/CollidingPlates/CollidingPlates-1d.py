# Non-compatible tests
#ATS:t10 = test(SELF, "--graphics False --restartStep 50 --clearDirectories True --domainIndependent True --compatibleEnergy False --outputFile 'CollidingPlates-1d-1proc-reproducing.txt' --referenceFile 'Reference/CollidingPlates-1d-reference-noncompatible-20160125.txt'", np=1, label="Colliding plates domain independence test SERIAL Non-compatible RUN")
#ATS:t11 = testif(t10, SELF, "--graphics False --clearDirectories False --domainIndependent True --compatibleEnergy False --outputFile 'CollidingPlates-1d-1proc-reproducing.txt' --referenceFile 'Reference/CollidingPlates-1d-reference-noncompatible-20160125.txt' --restoreCycle 50", np=1, label="Colliding plates domain independence test SERIAL RESTART Non-compatible RUN")
#ATS:t12 = testif(t10, SELF, "--graphics False --clearDirectories False --domainIndependent True --compatibleEnergy False --outputFile 'CollidingPlates-1d-4proc-reproducing.txt' --referenceFile 'Reference/CollidingPlates-1d-reference-noncompatible-20160125.txt' --comparisonFile 'dumps-CollidingPlates-1d/SolidSPHHydro1d/MonaghanGingoldViscosity1d/100/CollidingPlates-1d-1proc-reproducing.txt'", np=4, label="Colliding plates domain independence test 4 DOMAIN Non-compatible RUN")
#
# Compatible tests
#ATS:t50 = test(SELF, "--graphics False --restartStep 50 --dataDirBase 'dumps-CollidingPlates-compatible-1d' --clearDirectories True --domainIndependent True --compatibleEnergy True --outputFile 'CollidingPlates-1d-1proc-reproducing.txt' --referenceFile 'Reference/CollidingPlates-1d-reference-compatible-20160125.txt'", np=1, label="Colliding plates domain independence test SERIAL Compatible RUN")
#ATS:t51 = testif(t50, SELF, "--graphics False --dataDirBase 'dumps-CollidingPlates-compatible-1d' --clearDirectories False --domainIndependent True --compatibleEnergy True --outputFile 'CollidingPlates-1d-1proc-reproducing.txt' --referenceFile 'Reference/CollidingPlates-1d-reference-compatible-20160125.txt' --restoreCycle 50", np=1, label="Colliding plates domain independence test SERIAL RESTART Compatible RUN")
#ATS:t52 = testif(t50, SELF, "--graphics False --dataDirBase 'dumps-CollidingPlates-compatible-1d' --clearDirectories False --domainIndependent True --compatibleEnergy True --outputFile 'CollidingPlates-1d-4proc-reproducing.txt' --referenceFile 'Reference/CollidingPlates-1d-reference-compatible-20160125.txt' --comparisonFile 'dumps-CollidingPlates-compatible-1d/SolidSPHHydro1d/MonaghanGingoldViscosity1d/100/CollidingPlates-1d-1proc-reproducing.txt'", np=4, label="Colliding plates domain independence test 4 DOMAIN Compatible RUN")

#-------------------------------------------------------------------------------
# A pair of steel plates colliding at the origin.  This is a useful test of
# how our reflecting boundary condition handles problems with strength.
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
title("1-D colliding plates strength test")

#-------------------------------------------------------------------------------
# Generic problem parameters
# All CGS units.
#-------------------------------------------------------------------------------
commandLine(# Geometry
            length = 3.0,
            radius = 0.5,
            nx = 100,
            reflect = False,
            v0 = 1.0e4,        # Initial collision velocity

            # Material specific bounds on the mass density.
            rho0 = 7.9,
            etamin = 0.5,
            etamax = 1.5,

            # Hydro
            CRKSPH = False,
            Qconstructor = MonaghanGingoldViscosity,
            Cl = 1.0,
            Cq = 1.0,
            Qlimiter = False,
            linearInExpansion = False,
            balsaraCorrection = False,
            epsilon2 = 1e-2,
            negligibleSoundSpeed = 1e-5,
            csMultiplier = 1e-4,
            cfl = 0.5,
            useVelocityMagnitudeForDt = False,
            XSPH = True,
            epsilonTensile = 0.3,
            nTensile = 4,
            nPerh = 1.25,
            hmin = 1e-5,
            hmax = 1.0,
            HUpdate = IdealH,
            densityUpdate = IntegrateDensity,
            compatibleEnergy = True,
            evolveTotalEnergy = False,
            gradhCorrection = False,
            correctVelocityGradient = True,
            filter = 0.0,

            # Time advancement.
            IntegratorConstructor = CheapSynchronousRK2Integrator,
            goalTime = 2.0e-6,
            steps = None,
            dt = 1e-10,
            dtMin = 1e-12,
            dtMax = 1e-5,
            dtGrowth = 2.0,
            dumpFrac = 1.0,
            maxSteps = None,
            statsStep = 10,
            smoothIters = 0,
            domainIndependent = False,

            restoreCycle = None,
            restartStep = 500,

            # Diagnostics
            graphics = True,
            testtol = 1.0e-3,
            clearDirectories = False,
            referenceFile = "Reference/CollidingPlates-1d-reference-compatible-20160125.txt",
            dataDirBase = "dumps-CollidingPlates-1d",
            outputFile = "None",
            comparisonFile = "None",
            )

if CRKSPH:
    HydroConstructor = SolidCRKSPHHydro
    Qconstructor = CRKSPHMonaghanGingoldViscosity
else:
    HydroConstructor = SolidSPHHydro

dataDir = os.path.join(dataDirBase,
                       str(HydroConstructor).split("'")[1].split(".")[-1],
                       str(Qconstructor).split("'")[1].split(".")[-1],
                       str(nx))
restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "CollidingPlates-%i" % nx)

if reflect:
    xmin = 0.0
else:
    xmin = -0.5*length
xmax =  0.5*length

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
    T = ScalarField("temperature", nodes, 300.0)
    rho = nodes.massDensity()
    eps = nodes.specificThermalEnergy()
    eos.setSpecificThermalEnergy(eps, rho, T)

    # Set node velocites.
    pos = nodes.positions()
    vel = nodes.velocity()
    for i in xrange(nodes.numInternalNodes):
        if pos[i].x < 0.0:
            vel[i].x = v0
        else:
            vel[i].x = -v0

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
output("db.numSolidNodeLists")

#-------------------------------------------------------------------------------
# Construct the artificial viscosities for the problem.
#-------------------------------------------------------------------------------
q = Qconstructor(Cl, Cq, linearInExpansion)
q.limiter = Qlimiter
q.balsaraShearCorrection = balsaraCorrection
q.epsilon2 = epsilon2
q.negligibleSoundSpeed = negligibleSoundSpeed
q.csMultiplier = csMultiplier
output("q")
output("q.Cl")
output("q.Cq")
output("q.linearInExpansion")
output("q.limiter")
output("q.epsilon2")
output("q.negligibleSoundSpeed")
output("q.csMultiplier")
output("q.balsaraShearCorrection")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if CRKSPH:
    hydro = HydroConstructor(W = WT,
                             Q = q,
                             filter = filter,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             XSPH = XSPH,
                             densityUpdate = densityUpdate,
                             HUpdate = HUpdate)
else:
    hydro = HydroConstructor(W = WT,
                             Q = q,
                             filter = filter,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             evolveTotalEnergy = evolveTotalEnergy,
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
# Boundary conditions.
#-------------------------------------------------------------------------------
if reflect:
    plane0 = Plane(Vector(0.0), Vector(1.0))
    xbc0 = ReflectingBoundary(plane0)
    hydro.appendBoundary(xbc0)

#-------------------------------------------------------------------------------
# Construct a time integrator.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
integrator.appendPhysicsPackage(hydro)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.domainDecompositionIndependent = domainIndependent
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.domainDecompositionIndependent")

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle)
output("control")

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if steps is not None:
    control.step(steps)

    pos = nodes.positions()
    rho = nodes.massDensity()
    vel = nodes.velocity()
    eps = nodes.specificThermalEnergy()
    P = ScalarField("pressure", nodes)
    nodes.pressure(P)
    S = nodes.deviatoricStress()
    if reflect:
        points = [50, 0]
    else:
        points = [49, 50]
    for p in points:
        print pos[p], rho[p], vel[p], eps[p], P[p], S[p]

else:
    while control.time() < goalTime:
        nextGoalTime = min(control.time() + dtSample, goalTime)
        control.advance(nextGoalTime, maxSteps)
        control.dropRestartFile()

Eerror = (control.conserve.EHistory[-1] - control.conserve.EHistory[0])/control.conserve.EHistory[0]
print "Total energy error: %g" % Eerror
if compatibleEnergy and abs(Eerror) > 1e-13:
    raise ValueError, "Energy error outside allowed bounds."

#-------------------------------------------------------------------------------
# Plot the state.
#-------------------------------------------------------------------------------
if graphics:
    from SpheralGnuPlotUtilities import *
    state = State(db, integrator.physicsPackages())
    rhoPlot = plotFieldList(state.scalarFields("mass density"),
                            plotStyle = "linespoints",
                            winTitle = "rho @ %g %i" % (control.time(), mpi.procs))
    velPlot = plotFieldList(state.vectorFields("velocity"),
                            yFunction = "%s.x",
                            plotStyle = "linespoints",
                            winTitle = "vel @ %g %i" % (control.time(), mpi.procs))
    mPlot = plotFieldList(state.scalarFields("mass"),
                          plotStyle = "linespoints",
                          winTitle = "mass @ %g %i" % (control.time(), mpi.procs))
    PPlot = plotFieldList(state.scalarFields("pressure"),
                          plotStyle = "linespoints",
                          winTitle = "pressure @ %g %i" % (control.time(), mpi.procs))
    SPlot = plotFieldList(state.symTensorFields(SolidFieldNames.deviatoricStress),
                                                yFunction = "%s.xx",
                                                plotStyle = "linespoints",
                                                winTitle = "Deviatoric stress @ %g %i" % (control.time(), mpi.procs))

#-------------------------------------------------------------------------------
# If requested, write out the state in a global ordering to a file.
#-------------------------------------------------------------------------------
if outputFile != "None":
    state = State(db, integrator.physicsPackages())
    outputFile = os.path.join(dataDir, outputFile)
    pos = state.vectorFields(HydroFieldNames.position)
    rho = state.scalarFields(HydroFieldNames.massDensity)
    P = state.scalarFields(HydroFieldNames.pressure)
    vel = state.vectorFields(HydroFieldNames.velocity)
    eps = state.scalarFields(HydroFieldNames.specificThermalEnergy)
    Hfield = state.symTensorFields(HydroFieldNames.H)
    S = state.symTensorFields(SolidFieldNames.deviatoricStress)
    xprof = mpi.reduce([x.x for x in internalValues(pos)], mpi.SUM)
    rhoprof = mpi.reduce(internalValues(rho), mpi.SUM)
    Pprof = mpi.reduce(internalValues(P), mpi.SUM)
    vprof = mpi.reduce([v.x for v in internalValues(vel)], mpi.SUM)
    epsprof = mpi.reduce(internalValues(eps), mpi.SUM)
    hprof = mpi.reduce([1.0/sqrt(H.Determinant()) for H in internalValues(Hfield)], mpi.SUM)
    sprof = mpi.reduce([x.xx for x in internalValues(S)], mpi.SUM)
    mof = mortonOrderIndices(db)
    mo = mpi.reduce(internalValues(mof), mpi.SUM)
    if mpi.rank == 0:
        thpt = zip(mo, xprof, rhoprof, Pprof, vprof, epsprof, hprof, sprof)
        thpt = [thpt[i] for i in xrange(len(thpt)) if thpt[i][1] > 0.0]
        thpt.sort()
        mo, xprof, rhoprof, Pprof, vprof, epsprof, hprof, sprof = zip(*thpt)
        f = open(outputFile, "w")
        f.write(("#" + 7*" %16s" + "\n") % ("x", "rho", "P", "v", "eps", "h", "S"))
        for (xi, rhoi, Pi, vi, epsi, hi, si, mi) in zip(xprof, rhoprof, Pprof, vprof, epsprof, hprof, sprof, mo):
            f.write((7*"%16.12e " + "\n") %
                    (xi, rhoi, Pi, vi, epsi, hi, si))
        f.close()

        #---------------------------------------------------------------------------
        # Check the floating values for the state against reference data.
        #---------------------------------------------------------------------------
        if referenceFile != "None":
            import filearraycmp as fcomp
            assert fcomp.filearraycmp(outputFile, referenceFile, testtol, testtol)
            print "Floating point comparison test passed."

        #---------------------------------------------------------------------------
        # Also we can optionally compare the current results with another file for
        # bit level consistency.
        #---------------------------------------------------------------------------
        if comparisonFile != "None":
            import filecmp
            print "Compare files : %s     <--->     %s" % (outputFile, comparisonFile)
            assert filecmp.cmp(outputFile, comparisonFile)
