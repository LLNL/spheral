#-------------------------------------------------------------------------------
# Use forced terms for the inputs to the P-alpha model (pressure, derivatives,
# hydro time derivatives, etc.) to force the P-alpha distension evolution to
# follow a known analytic solution.
#-------------------------------------------------------------------------------
from SolidSpheral1d import *
from SpheralTestUtilities import *
from math import *
import os, shutil
import mpi

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("Analytic evolution test of P-alpha porosity model")

#-------------------------------------------------------------------------------
# Fake hydro class to force known time derviative values
#-------------------------------------------------------------------------------
class FakeHydro(Physics):
    def __init__(self,
                 db,
                 DrhoDt0,
                 DuDt0):
        Physics.__init__(self)
        self.DrhoDt0 = DrhoDt0
        self.DuDt0 = DuDt0
        self.pressure = db.newFluidScalarFieldList(0.0, HydroFieldNames.pressure)
        self.DrhoDt = db.newFluidScalarFieldList(0.0, "delta " + HydroFieldNames.massDensity)
        self.DuDt = db.newFluidScalarFieldList(0.0, "delta " + HydroFieldNames.specificThermalEnergy)
        return

    def evaluateDerivatives(self, t, dt, db, state, derivs):
        DrhoDt = derivs.scalarFields("delta " + HydroFieldNames.massDensity)
        DuDt = derivs.scalarFields("delta " + HydroFieldNames.specificThermalEnergy)
        assert (DrhoDt.numFields == 1 and DuDt.numFields == 1 and
                DrhoDt[0].numElements == 1 and DuDt[0].numElements == 1)
        DrhoDt[0][0] = self.DrhoDt0
        DuDt[0][0] = self.DuDt0
        return

    def dt(self, db, state, derivs, t):
        print("FakeHydro::dt ", state.keys())
        rho = state.scalarFields(HydroFieldNames.massDensity)
        print("              ", rho)
        print("              ", rho[0].internalValues())
        return pair_double_string(1.0, "Fake timestep")

    def registerState(self, db, state):
        rho = db.fluidMassDensity
        state.enroll(rho, ScalarIncrementBoundedFieldList(0.0))
        state.enroll(db.fluidSpecificThermalEnergy, ScalarIncrementFieldList())
        state.enroll(self.pressure, PressurePolicy())
        print("All keys registered in state following FakeHydro::registerState: ", state.keys())
        return

    def registerDerivatives(self, db, derivs):
        derivs.enroll(self.DrhoDt)
        derivs.enroll(self.DuDt)
        print("All keys registered in derivs following FakeHydro::registerDerivatives: ", derivs.keys())
        return

    def label(self):
        return "FakeHydro"

#-------------------------------------------------------------------------------
# Build our units (cm, gm, microsec)
#-------------------------------------------------------------------------------
units = CGuS()

#-------------------------------------------------------------------------------
# Generic problem parameters
# All (cm, gm, usec) units.
#-------------------------------------------------------------------------------
commandLine(

    # Porosity
    alpha0 = 1.275,                    # Initial distention
    alphae = 1.2,                      # Elastic limiting distention
    alphat = 1.1,                      # Plastic limiting distention
    c0frac = 0.5,                      # Fraction of initial solid sound speed for initial porous sound speed
    Pe = 0.1,                          # Elastic limiting pressure
    Pt = 0.5,                          # Transition pressure in plastic regime
    n1 = 1.0,
    n2 = 1.0,

    # Fake hydro
    dPdrho0 = 0.5,
    dPdu0 = 0.75,
    DrhoDt0 = 1.5,
    DuDt0 = 0.1,

    # Initial conditions (initial density set by alpha0 and rhoS0)
    eps0 = 0.0,

    # Linear-polynomial EOS
    rhoS0 = 2.0,
    a0 = 0.0,
    a1 = 1.5,
    a2 = 0.0,
    a3 = 0.0,
    b0 = 0.75,
    b1 = 0.0,
    b2 = 0.0,
    mu = 2.0,

    IntegratorConstructor = SynchronousRK1Integrator,
    goalTime = 1.0,
    steps = None,
    dt = 1e-2,
    dtMin = 1e-2,
    dtMax = 1e-2,
    dtverbose = False,
    graphics = True,
    clearDirectories = False,
    dataDir = "dumps-PalphaAnalyticEvolution",
)

assert 1.0 <= alpha0
assert 1.0 <= alphat <= alphae <= alpha0
assert c0frac <= 1.0
assert Pe <= Pt

phi0 = 1.0 - 1.0/alpha0

# restartDir = os.path.join(dataDir, "restarts")
# restartBaseName = os.path.join(restartDir, "PalphaAnalyticEvolution")

# #-------------------------------------------------------------------------------
# # Check if the necessary output directories exist.  If not, create them.
# #-------------------------------------------------------------------------------
# import os, sys
# if mpi.rank == 0:
#     if clearDirectories and os.path.exists(dataDir):
#         shutil.rmtree(dataDir)
#     if not os.path.exists(restartDir):
#         os.makedirs(restartDir)
# mpi.barrier()

#-------------------------------------------------------------------------------
# Material models
#-------------------------------------------------------------------------------
eosS = LinearPolynomialEquationOfState(referenceDensity = rhoS0,
                                       etamin = 0.0,
                                       etamax = 1e5,
                                       a0 = a0,
                                       a1 = a1,
                                       a2 = a2,
                                       a3 = a3,
                                       b0 = b0,
                                       b1 = b1,
                                       b2 = b2,
                                       atomicWeight = mu,
                                       constants = units)

strengthModelS = NullStrength()
eos = PorousEquationOfState(eosS)
strengthModel = PorousStrengthModel(strengthModelS)
eps0 = 0.0
Ps0 = eosS.pressure(rhoS0, eps0)
cS0 = eosS.soundSpeed(rhoS0, eps0)
print("Material initial properties:")
output(" rhoS0")
output("  eps0")
output("   Ps0")
output("   cS0")

#-------------------------------------------------------------------------------
# Create our interpolation kernels
#-------------------------------------------------------------------------------
WT = TableKernel(WendlandC4Kernel(), 200)
output("WT")

#-------------------------------------------------------------------------------
# Create the NodeLists.
#-------------------------------------------------------------------------------
nodes = makeSolidNodeList("aluminium", eos, strengthModel,
                          numInternal = 1,
                          nPerh = 4.01,
                          xmin = -10.0*Vector.one,
                          xmax =  10.0*Vector.one)
nodeSet = [nodes]

#-------------------------------------------------------------------------------
# Set node properties (positions, masses, H's, etc.)
#-------------------------------------------------------------------------------
pos = nodes.positions()
vel = nodes.velocity()
rho = nodes.massDensity()
eps = nodes.specificThermalEnergy()
H = nodes.Hfield()
rho[0] = rhoS0 * (1.0 - phi0)
eps[0] = eps0
H[0].xx = 1.0

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
# Construct the fake hydro physics object.  We force known derivatives for
# the hydro state in order to drive and test the porosity evolution.
#-------------------------------------------------------------------------------
hydro = FakeHydro(db, DrhoDt0, DuDt0)
output("hydro")
output("hydro.DrhoDt0")
output("hydro.DuDt0")

packages = [hydro]

#-------------------------------------------------------------------------------
# Construct a porosity model
#-------------------------------------------------------------------------------
phi0 = 1.0 - 1.0/alpha0
c0 = c0frac * cS0

print("P-alpha porosity model parameters:")
output("  alpha0")
output("    phi0")
porosityAl = PalphaPorosity(eos, strengthModel, nodes,
                            phi0 = phi0,
                            Pe = Pe,
                            Pt = Pt,
                            alphae = alphae,
                            alphat = alphat,
                            n1 = n1,
                            n2 = n2,
                            cS0 = cS0,
                            c0 = c0)

output("porosityAl")
output("  porosityAl.Pe")
output("  porosityAl.Pt")
output("  porosityAl.alphae")
output("  porosityAl.alphat")
output("  porosityAl.n1")
output("  porosityAl.n2")
output("  porosityAl.cS0")
output("  porosityAl.phi().internalValues()")
packages.append(porosityAl)

#-------------------------------------------------------------------------------
# Construct a time integrator.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
for package in packages:
    integrator.appendPhysicsPackage(package)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.verbose = dtverbose
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.havePhysicsPackage(porosityAl)")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.domainDecompositionIndependent")

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            restartStep = None,
                            iterateInitialH = False)
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
# Plot the state.
#-------------------------------------------------------------------------------
# if graphics:
#     from SpheralMatplotlib import *
#     state = State(db, integrator.physicsPackages())
#     H = state.symTensorFields("H")
#     h = db.newFluidScalarFieldList(0.0, "h")
#     for i in range(nodes.numInternalNodes):
#         h[0][i] = 1.0/H[0][i].xx
#     rhoPlot = plotFieldList(state.scalarFields("mass density"),
#                             plotStyle = "o-",
#                             winTitle = "rho @ %g %i" % (control.time(), mpi.procs))
#     velPlot = plotFieldList(state.vectorFields("velocity"),
#                             yFunction = "%s.x",
#                             plotStyle = "o-",
#                             winTitle = "vel @ %g %i" % (control.time(), mpi.procs))
#     PPlot = plotFieldList(state.scalarFields("pressure"),
#                           plotStyle = "o-",
#                           winTitle = "pressure @ %g %i" % (control.time(), mpi.procs))
#     uPlot = plotFieldList(state.scalarFields("specific thermal energy"),
#                           plotStyle = "o-",
#                           winTitle = "specific thermal energy @ %g %i" % (control.time(), mpi.procs))
#     SPlot = plotFieldList(state.symTensorFields("deviatoric stress"),
#                           yFunction = "%s.xx",
#                           plotStyle = "o-",
#                           winTitle = "$S_{xx}$ @ %g %i" % (control.time(), mpi.procs))
#     hPlot = plotFieldList(h,
#                           plotStyle = "o-",
#                           winTitle = "h @ %g %i" % (control.time(), mpi.procs))
#     csPlot = plotFieldList(state.scalarFields(HydroFieldNames.soundSpeed),
#                            plotStyle = "o-",
#                            xlabel = r"$x$",
#                            ylabel = r"$c_S$",
#                            winTitle = "Sound speed @ %g %i" %  (control.time(), mpi.procs))
#     alpha = porosityAl.alpha
#     alphaPlot = plotField(alpha,
#                           plotStyle = "o-",
#                           winTitle = r"$\alpha$ @ %g %i" %  (control.time(), mpi.procs))
#     DalphaDt = porosityAl.DalphaDt
#     DalphaDtPlot = plotField(DalphaDt,
#                              plotStyle = "o-",
#                              xlabel = r"$x$",
#                              ylabel = r"$D\alpha/Dt$",
#                              winTitle = r"$D\alpha/Dt$ @ %g %i" %  (control.time(), mpi.procs))
#     phi = porosityAl.phi()
#     phiPlot = plotField(phi,
#                         plotStyle = "o-",
#                         winTitle = r"$\phi$ @ %g %i" %  (control.time(), mpi.procs))

#     dPdRplot = plotField(porosityAl.partialPpartialRho,
#                          plotStyle = "o-",
#                          winTitle = r"$\partial P/\partial \rho$ @ %g %i" %  (control.time(), mpi.procs))
#     dPdUplot = plotField(porosityAl.partialPpartialEps,
#                          plotStyle = "o-",
#                          winTitle = r"$\partial P/\partial \varepsilon$ @ %g %i" %  (control.time(), mpi.procs))

#     plots = [(rhoPlot, "rho.png"),
#              (velPlot, "vel.png"),
#              (PPlot, "pressure.png"),
#              (SPlot, "devstress.png"),
#              (uPlot, "u.png"),
#              (hPlot, "h.png"),
#              (alphaPlot, "alpha.png"),
#              (phiPlot, "phi.png")]

#     # Save the figures.
#     for p, fname in plots:
#         savefig(p, os.path.join(dataDir, fname))

# if graphics:
#     plt.show()
