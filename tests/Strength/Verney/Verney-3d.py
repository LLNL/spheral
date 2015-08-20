#-------------------------------------------------------------------------------
# An idealized strength test test where an imploding shell is ultimately stopped
# at a known radius due to plastic work dissipation.
#
# See LA-14379
#-------------------------------------------------------------------------------
from math import *
import shutil
import mpi

from SolidSpheral3d import *
from SpheralTestUtilities import *
from findLastRestart import *
from NodeHistory import NodeHistory
from GenerateNodeDistribution3d import *
from VoronoiDistributeNodes import distributeNodes3d

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("3-D Verney imploding shell test.")

#-------------------------------------------------------------------------------
# Use (cm, gm, usec) units as in paper.
#-------------------------------------------------------------------------------
units = PhysicalConstants(0.01,  # Unit length in m
                          0.001, # Unit mass in kg
                          1e-6)  # Unit length in sec

#-------------------------------------------------------------------------------
# Generic problem parameters
# All (cm, gm, usec) units.
#-------------------------------------------------------------------------------
commandLine(nr = 10,              # Radial resolution of the shell in points
            seed = "lattice",     # "lattice" or "icosahedral"
            nPerh = 2.01,

            # Material specific bounds on the mass density.
            etamin = 1e-3,
            etamax = 1e3,

            # Should we run with strength?
            useStrength = True,

            # Hydro parameters.
            CRKSPH = False,
            Qconstructor = MonaghanGingoldViscosity,
            Cl = 1.0,
            Cq = 1.0,
            linearInExpansion = False,
            Qlimiter = False,
            balsaraCorrection = False,
            epsilon2 = 1e-2,
            negligibleSoundSpeed = 1e-5,
            csMultiplier = 1e-4,
            hmin = 1e-5,
            hmax = 0.1,
            cfl = 0.25,
            useVelocityMagnitudeForDt = False,
            XSPH = False,
            epsilonTensile = 0.0,
            nTensile = 4,
            filter = 0.0,
            HUpdate = IdealH,
            densityUpdate = IntegrateDensity,
            compatibleEnergy = True,
            gradhCorrection = False,

            # Time integration
            IntegratorConstructor = CheapSynchronousRK2Integrator,
            goalTime = 100.0,
            steps = None,
            dt = 1e-6,
            dtMin = 1e-6,
            dtMax = 10.0,
            dtGrowth = 2.0,
            maxSteps = None,
            statsStep = 10,
            domainIndependent = False,
            dtverbose = False,
            restoreCycle = -1,
            restartStep = 500,
            sampleFreq = 10,

            graphics = True,

            clearDirectories = False,
            dataDirBase = "dumps-Verney-Cu-3d",
            outputFile = "Verney-Cu-3d.gnu",
        )

# Material parameters for this test problem.
rho0Cu = 9.30333
R0, R1 = 5.0, 10.0        # Inner, outer initial radius
G0, Y0 = 0.48, 0.0012     # Shear modulus and yield strength
r0 = 2.75                 # Expected final inner radius
delta = R1 - R0
lamb = r0/R0
alpha = delta/R0
beta = lamb**3 + alpha**3 + 3.0*alpha*(alpha* + 1.0)
F = (1.0 + alpha)**3*log(1.0 + alpha) + lamb**3*log(lamb) + beta*log(beta)/3.0
U0 = sqrt(4.0*Y0/(3.0*rho0Cu)*R1/delta*F)

# Hydro constructor.
if CRKSPH:
    HydroConstructor = SolidCRKSPHHydro
    Qconstructor = CRKSPHMonaghanGingoldViscosity
else:
    HydroConstructor = SolidSPHHydro

# Directories.
dataDir = os.path.join(dataDirBase,
                       "strength=%s" % useStrength,
                       str(HydroConstructor).split("'")[1].split(".")[-1],
                       str(Qconstructor).split("'")[1].split(".")[-1],
                       "densityUpdate=%s" % densityUpdate,
                       "nr=%i" % nr)
restartDir = os.path.join(dataDir, "restarts")
visitDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, "Verney-%i" % nr)

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
eosCu = GruneisenEquationOfState(rho0Cu,  # reference density  
                                 etamin,  # etamin             
                                 etamax,  # etamax             
                                 0.3940,  # C0                 
                                 1.489,   # S1                 
                                 0.0,     # S2                 
                                 0.0,     # S3                 
                                 2.02,    # gamma0             
                                 0.47,    # b                  
                                 63.57,   # atomic weight
                                 units)
strengthModelCu = ConstantStrength(0.48,        # G0
                                   0.0012)      # Y0

#-------------------------------------------------------------------------------
# If we're not using strength, override the strength models.
#-------------------------------------------------------------------------------
if not useStrength:
    strengthModelCu = NullStrength()

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
nodesCu = makeSolidNodeList("Copper", eosCu, strengthModelCu,
                            nPerh = nPerh,
                            hmin = hmin,
                            hmax = hmax,
                            rhoMin = etamin*rho0Cu,
                            rhoMax = etamax*rho0Cu,
                            xmin = -100.0*Vector.one,
                            xmax =  100.0*Vector.one)
nodeSet = [nodesCu]

#-------------------------------------------------------------------------------
# Set node properties (positions, masses, H's, etc.)
#-------------------------------------------------------------------------------
print "Generating node distribution."
from DistributeNodes import distributeNodesInRange1d
if seed == "lattice":
    nx = int(R1/(R1 - R0) * nr + 0.5)
    gen = GenerateNodeDistribution3d(2*nx, 2*nx, 2*nx, rho0Cu,
                                     distributionType = seed,
                                     xmin = (-R1, -R1, -R1),
                                     xmax = ( R1,  R1,  R1),
                                     rmin = R0, 
                                     rmax = R1)
elif seed == "icosohedral":
    gen = GenerateIcosahedronMatchingProfile(nr, rho0Cu,
                                             rmin = R0,
                                             rmax = R1,
                                             nNodePerh = nPerh)

distributeNodes3d((nodesCu, gen))
output("mpi.reduce(nodesCu.numInternalNodes, mpi.MIN)")
output("mpi.reduce(nodesCu.numInternalNodes, mpi.MAX)")
output("mpi.reduce(nodesCu.numInternalNodes, mpi.SUM)")

# Set node velocites.
pos = nodesCu.positions()
vel = nodesCu.velocity()
for i in xrange(nodesCu.numInternalNodes):
    ri = pos[i].magnitude()
    rhat = pos[i].unitVector()
    vel[i] = U0 * (R0/ri)**2 * rhat

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
    hydro = HydroConstructor(WT, WTPi, q,
                             filter = filter,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             XSPH = XSPH,
                             densityUpdate = densityUpdate,
                             HUpdate = HUpdate)
else:
    hydro = HydroConstructor(WT, WTPi, q,
                             filter = filter,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             gradhCorrection = gradhCorrection,
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
integrator.verbose = dtverbose
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
if not steps is None:
    control.step(steps)
else:
    control.advance(goalTime)
    control.dropRestartFile()

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
    rprof = mpi.reduce([x.magnitude() for x in internalValues(pos)], mpi.SUM)
    rhoprof = mpi.reduce(internalValues(rho), mpi.SUM)
    Pprof = mpi.reduce(internalValues(P), mpi.SUM)
    vprof = mpi.reduce([v.x for v in internalValues(vel)], mpi.SUM)
    epsprof = mpi.reduce(internalValues(eps), mpi.SUM)
    hprof = mpi.reduce([1.0/sqrt(H.Determinant()) for H in internalValues(Hfield)], mpi.SUM)
    sprof = mpi.reduce([x.xx for x in internalValues(S)], mpi.SUM)
    mof = mortonOrderIndices(db)
    mo = mpi.reduce(internalValues(mof), mpi.SUM)
    if mpi.rank == 0:
        multiSort(mo, rprof, rhoprof, Pprof, vprof, epsprof, hprof, sprof)
        f = open(outputFile, "w")
        f.write(("#" + 15*" %16s" + "\n") % ("r", "rho", "P", "v", "eps", "h", "S", "m", 
                                             "int(r)", "int(rho)", "int(P)", "int(v)", "int(eps)", "int(h)", "int(S)"))
        for (ri, rhoi, Pi, vi, epsi, hi, si, mi) in zip(rprof, rhoprof, Pprof, vprof, epsprof, hprof, sprof, mo):
            f.write((7*"%16.12e " + 8*"%i " + "\n") %
                    (ri, rhoi, Pi, vi, epsi, hi, si, mi,
                     unpackElementUL(packElementDouble(ri)),
                     unpackElementUL(packElementDouble(rhoi)),
                     unpackElementUL(packElementDouble(Pi)),
                     unpackElementUL(packElementDouble(vi)),
                     unpackElementUL(packElementDouble(epsi)),
                     unpackElementUL(packElementDouble(hi)),
                     unpackElementUL(packElementDouble(si))))
        f.close()

#-------------------------------------------------------------------------------
# Plot the state.
#-------------------------------------------------------------------------------
if graphics:
    from SpheralGnuPlotUtilities import *
    state = State(db, integrator.physicsPackages())
    rhoPlot = plotFieldList(state.scalarFields("mass density"),
                            xFunction = "%s.magnitude()",
                            plotStyle="points",
                            winTitle="rho @ %g" % (control.time()))
    velPlot = plotFieldList(state.vectorFields("velocity"),
                            xFunction = "%s.magnitude()",
                            yFunction = "%s.magnitude()",
                            plotStyle="points",
                            winTitle="vel @ %g" % (control.time()))
    mPlot = plotFieldList(state.scalarFields("mass"),
                          xFunction = "%s.magnitude()",
                          plotStyle="points",
                          winTitle="mass @ %g" % (control.time()))
    PPlot = plotFieldList(state.scalarFields("pressure"),
                          xFunction = "%s.magnitude()",
                          plotStyle="points",
                          winTitle="pressure @ %g" % (control.time()))
    hPlot = plotFieldList(state.symTensorFields("H"),
                          xFunction = "%s.magnitude()",
                          yFunction = "1.0/%s.xx",
                          plotStyle="points",
                          winTitle="h @ %g" % (control.time()))
