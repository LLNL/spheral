#-------------------------------------------------------------------------------
# An idealized strength test test mocking a flyer plate of Aluminum impacting 
# a fixed plate of Tantalum.
#
# See Weseloh & Najjar, LLNL-TR-484921
#-------------------------------------------------------------------------------
from SolidSpheral1d import *
from SpheralTestUtilities import *
from findLastRestart import *
from NodeHistory import NodeHistory
from math import *
import shutil
import mpi

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("1-D TP106 flyer plate test.")

#-------------------------------------------------------------------------------
# Use (cm, gm, usec) units as in paper.
#-------------------------------------------------------------------------------
units = PhysicalConstants(0.01,  # Unit length in m
                          0.001, # Unit mass in kg
                          1e-6)  # Unit length in sec

# #-------------------------------------------------------------------------------
# # The strength model used in this test.
# #-------------------------------------------------------------------------------
# class TP106StrengthModel(StrengthModel):

#     def __init__(self, 
#                  rho0,    # reference density
#                  G0,      # initial shear modulus
#                  Gmax,    # maximum shear modulus
#                  Y0,      # initial yield
#                  Ymax,    # maximum yield
#                  alpha,   # strain hardening coeff
#                  ps0,     # initial plastic strain
#                  beta,    # strain hardening exponent, thermal softening coeff
#                  Em,      # melt specific energy
#                  gamma,   # pressure hardening coeff for yield
#                  gammap): # pressure hardening coeff for shear modulus
#         StrengthModel.__init__(self)
#         self.rho0 = rho0
#         self.Gmax = Gmax
#         self.Y0 = Y0
#         self.Ymax = Ymax
#         self.alpha = alpha
#         self.ps0 = ps0
#         self.beta = beta
#         self.Em = Em
#         self.gamma = gamma
#         self.gammap = gammap
#         return

#     def shearModulus(shearModulus,
#                      density,
#                      specificEnergy,
#                      pressure):
#         n = shearModulus.numInternalElements
#         for i in xrange(n):
#             G[i] = min(self.Gmax,
#                        self.G0*(1.0 + self.gamma*pressure[i]*(density[i]/self.rho0)**(-1.0/3.0))*exp(-

#         return

#-------------------------------------------------------------------------------
# Generic problem parameters
# All (cm, gm, usec) units.
#-------------------------------------------------------------------------------
commandLine(nxAl = 100,
            nPerh = 2.01,

            # Material specific bounds on the mass density.
            etamin = 0.2,
            etamax = 1.8,

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
            XSPH = True,
            epsilonTensile = 0.0,
            nTensile = 4,
            filter = 0.0,
            HUpdate = IdealH,
            densityUpdate = IntegrateDensity,
            compatibleEnergy = True,
            gradhCorrection = False,

            # Time integration
            IntegratorConstructor = CheapSynchronousRK2Integrator,
            goalTime = 3.0,
            steps = None,
            dt = 1e-10,
            dtMin = 1e-6,
            dtMax = 10.0,
            dtGrowth = 2.0,
            dumpFrac = 0.005,
            maxSteps = None,
            statsStep = 10,
            domainIndependent = False,
            dtverbose = False,
            restoreCycle = -1,
            restartStep = 500,
            sampleFreq = 10,

            graphics = True,

            clearDirectories = False,
            dataDirBase = "dumps-TP106-1d",
            )

# Figure out the mass matched resolution.
rho0Al = 2.785
rho0Ta = 16.654
nxTa = int(rho0Ta/rho0Al * nxAl)
print "Selected %i Ta points to mass match with %i Al points." % (nxTa, nxAl)

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
                       "nxAl=%i_nxTa=%i" % (nxAl, nxTa))
restartDir = os.path.join(dataDir, "restarts")
visitDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, "TP106-%i" % (nxAl + nxTa))
tracerOutputName = os.path.join(dataDir, "TP106-tracer%i.dat")

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
# Al
eosAl = GruneisenEquationOfState(rho0Al,  # reference density  
                                 etamin,  # etamin             
                                 etamax,  # etamax             
                                 0.5328,  # C0                 
                                 1.338,   # S1                 
                                 0.0,     # S2                 
                                 0.0,     # S3                 
                                 2.0,     # gamma0             
                                 0.0,     # b                  
                                 26.98,  # atomic weight
                                 units)
coldFit = NinthOrderPolynomialFit(0.0,
                                  0.0,
                                  0.0,
                                  0.0,
                                  0.0,
                                  0.0,
                                  0.0,
                                  0.0,
                                  0.0,
                                  0.0)
meltFitAl = NinthOrderPolynomialFit(0.0796,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0)
strengthModelAl = SteinbergGuinanStrength(eosAl,
                                          0.2860,        # G0
                                          6.52,          # A
                                          0.0,           # B
                                          0.0026,        # Y0
                                          0.0076,        # Ymax
                                          1.0e-3,        # Yp (delta)
                                          310.0,         # beta
                                          0.0,           # gamma0
                                          0.185,         # nhard
                                          coldFit,
                                          meltFitAl)

# Ta
eosTa = GruneisenEquationOfState(rho0Ta,  # reference density  
                                 etamin,  # etamin             
                                 etamax,  # etamax             
                                 0.3414,  # C0                 
                                 1.201,   # S1                 
                                 0.0,     # S2                 
                                 0.0,     # S3                 
                                 2.0,     # gamma0             
                                 0.0,     # b                  
                                 180.95,  # atomic weight
                                 units)
meltFitTa = NinthOrderPolynomialFit(0.0546,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0)
strengthModelTa = SteinbergGuinanStrength(eosTa,
                                          0.6900,        # G0
                                          1.45,          # A
                                          0.0,           # B
                                          0.0077,        # Y0
                                          0.0110,        # Ymax
                                          1.0e-3,        # Yp (delta)
                                          10.0,          # beta
                                          0.0,           # gamma0
                                          0.1,           # nhard
                                          coldFit,
                                          meltFitTa)

#-------------------------------------------------------------------------------
# If we're not using strength, override the strength models.
#-------------------------------------------------------------------------------
if not useStrength:
    strengthModelAl = NullStrength()
    strengthModelTa = NullStrength()

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
nodesAl = makeSolidNodeList("Aluminum", eosAl, strengthModelAl,
                            nPerh = nPerh,
                            hmin = hmin,
                            hmax = hmax,
                            rhoMin = etamin*rho0Al,
                            rhoMax = etamax*rho0Al,
                            xmin = -100.0*Vector.one,
                            xmax =  100.0*Vector.one)
nodesTa = makeSolidNodeList("Tantalum", eosTa, strengthModelTa,
                            nPerh = nPerh,
                            hmin = hmin,
                            hmax = hmax,
                            rhoMin = etamin*rho0Ta,
                            rhoMax = etamax*rho0Ta,
                            xmin = -100.0*Vector.one,
                            xmax =  100.0*Vector.one)
nodeSet = [nodesAl, nodesTa]

#-------------------------------------------------------------------------------
# Set node properties (positions, masses, H's, etc.)
#-------------------------------------------------------------------------------
print "Generating node distribution."
from DistributeNodes import distributeNodesInRange1d
distributeNodesInRange1d([(nodesAl, nxAl, rho0Al, (-1.025, -0.025)),
                          (nodesTa, nxTa, rho0Ta, ( 0.025,  1.025))])

# Set node velocites.
vel = nodesAl.velocity()
for i in xrange(nodesAl.numInternalNodes):
    vel[i].x = 0.18

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
# Select points to follow histories.
#-------------------------------------------------------------------------------
def tp106tracersample(nodes, indices):
    assert mpi.allreduce(len(indices), mpi.SUM) == 1
    if len(indices) == 1:
        i = indices[0]
        assert i < nodes.numInternalNodes
        pos = nodes.positions()
        vel = nodes.velocity()
        rho = nodes.massDensity()
        eps = nodes.specificThermalEnergy()
        P = ScalarField("pressure", nodes)
        nodes.pressure(P)
        result = [pos[i].x, vel[i].x, rho[i], eps[i], P[i]]
    else:
        result = []
    result = mpi.allreduce(result, mpi.SUM)
    return tuple(result)

# Al sample points.
histories = []
tracerNumber = 1
for nodes, samplePositions in ((nodesAl, (-0.0375, -0.0625, -0.9875)),
                               (nodesTa, ( 0.0375,  0.0625,  0.9875))):
    pos = nodes.positions()
    for x0 in samplePositions:
        thpt = [(abs(pos[i].x - x0), i) for i in xrange(nodes.numInternalNodes)] + [(1e100, -1)]
        thpt.sort()
        dxmin = mpi.allreduce(thpt[0][0], mpi.MIN)
        if thpt[0][0] == dxmin:
            i = thpt[0][1]
            indices = [i]
            sys.stderr.write("Tracer %i is node %i @ %s.\n" % (tracerNumber, i, pos[i]))
        else:
            indices = []
        histories.append(NodeHistory(nodes, indices, tp106tracersample, tracerOutputName % tracerNumber, 
                                     labels = ("pos", "vel", "rho", "eps", "P")))
        tracerNumber += 1

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
# Add the diagnostics to the controller.
#-------------------------------------------------------------------------------
for hist in histories:
    control.appendPeriodicWork(hist.sample, sampleFreq)
    hist.flushHistory()

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
    hPlot = plotFieldList(state.symTensorFields("H"),
                          yFunction = "1.0/%s.xx",
                          plotStyle="linespoints",
                          winTitle="h @ %g %i" % (control.time(), mpi.procs))
