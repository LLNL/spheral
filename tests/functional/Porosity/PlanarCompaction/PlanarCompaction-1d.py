#-------------------------------------------------------------------------------
# A rod of aluminum undergoing compaction by a piston.  This test is based on the
# first example test from
#
# Jutzi, M., Benz, W., & Michel, P. (2008). Numerical simulations of impacts
# involving porous bodies.I. Implementing sub-resolution porosity in a 3D SPH
# hydrocode. Icarus, 198(1), 242–255.
#-------------------------------------------------------------------------------
from SolidSpheral1d import *
from SpheralTestUtilities import *
from math import *
import os, shutil
import mpi

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("1-D planar compaction of a poroous Al rod")

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
        ids = [i for i in range(self.nodeList.numInternalNodes)
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
# Build our units (cm, gm, microsec)
#-------------------------------------------------------------------------------
units = CGuS()

#-------------------------------------------------------------------------------
# Generic problem parameters
# All (cm, gm, usec) units.
#-------------------------------------------------------------------------------
commandLine(nx = 500,                          # Number of internal free points
            vright = -45.8e-3,                 # Fixed velocity on right-end
            vleft = 0.0,                       # Fixed velocity on left-end
            nbound = 10,                       # Number of points in each fixed velocity boundaries

            PorousModel = PalphaPorosity,

            crksph = False,
            cfl = 0.5,
            useVelocityMagnitudeForDt = False,
            XSPH = False,
            epsilonTensile = 0.3,
            nTensile = 4,
            volumeType = RKSumVolume,

            IntegratorConstructor = VerletIntegrator,
            goalTime = 3.5,
            steps = None,
            dt = 1e-10,
            dtMin = 1e-6,
            dtMax = 10.0,
            dtGrowth = 2.0,
            maxSteps = None,
            statsStep = 10,
            densityUpdate = IntegrateDensity,
            compatibleEnergy = True,
            domainIndependent = True,
            dtverbose = False,

            restoreCycle = -1,
            restartStep = 500,

            graphics = True,

            clearDirectories = False,
            dataDirBase = "dumps-PlanarCompaction-1d",
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

dataDir = os.path.join(dataDirBase,
                       PorousModel.__name__,
                       hydroname,
                       "nx=%i" % nx)

restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "PlanarCompaction-%i" % nx)

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
# Al material parameters
# Note, in the test as presented in the reference there is no strength model,
# just pressure response.
#-------------------------------------------------------------------------------
eosS = TillotsonEquationOfState("aluminum",
                                etamin = 0.1,
                                etamax = 10.0,
                                units = units)
strengthModelS = NullStrength()
eos = PorousEquationOfState(eosS)
strengthModel = PorousStrengthModel(strengthModelS)
rho0 = eosS.referenceDensity
eps0 = 0.0
Ps0 = eosS.pressure(rho0, eps0)
cS0 = eosS.soundSpeed(rho0, eps0)
print("Aluminum equation of state initial properties:")
output("  rho0")
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
                          nPerh = nPerh,
                          xmin = -10.0*Vector.one,
                          xmax =  10.0*Vector.one)
nodeSet = [nodes]

#-------------------------------------------------------------------------------
# Set node properties (positions, masses, H's, etc.)
#-------------------------------------------------------------------------------
print("Generating node distribution.")
from DistributeNodes import distributeNodesInRange1d
dx = 2.0/nx
xmax = 1.0 + (nbound + 0.5)*dx
distributeNodesInRange1d([(nodes, nx, rho0, (-xmax, xmax))])
output("mpi.reduce(nodes.numInternalNodes, mpi.MIN)")
output("mpi.reduce(nodes.numInternalNodes, mpi.MAX)")
output("mpi.reduce(nodes.numInternalNodes, mpi.SUM)")

# Set node velocites and find the boundary nodes
pos = nodes.positions()
vel = nodes.velocity()
xmin_nodes = [i for i in range(nodes.numInternalNodes) if pos[i].x < -1.0]
xmax_nodes = [i for i in range(nodes.numInternalNodes) if pos[i].x >  1.0]
for i in xmin_nodes:
    vel[i].x = vleft
for i in xmax_nodes:
    vel[i].x = vright

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
# Construct constant velocity boundary conditions to be applied to the rod ends.
#-------------------------------------------------------------------------------
print("Selected %i constant velocity nodes." % (mpi.allreduce(len(xmin_nodes) + len(xmax_nodes), mpi.SUM)))
xbc0 = ConstantVelocityBoundary(nodes, xmin_nodes)
xbc1 = ConstantVelocityBoundary(nodes, xmax_nodes)

bcs = [xbc0, xbc1]

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if crksph:
    hydro = CRKSPH(dataBase = db,
                   cfl = cfl,
                   compatibleEnergyEvolution = compatibleEnergy,
                   XSPH = XSPH,
                   densityUpdate = densityUpdate)
else:
    hydro = SPH(dataBase = db,
                W = WT,
                cfl = cfl,
                compatibleEnergyEvolution = compatibleEnergy,
                densityUpdate = densityUpdate,
                XSPH = XSPH,
                epsTensile = epsilonTensile,
                nTensile = nTensile)
output("hydro")
output("hydro.cfl")
output("hydro.useVelocityMagnitudeForDt")
output("hydro.HEvolution")
output("hydro.densityUpdate")
output("hydro.compatibleEnergyEvolution")

packages = [hydro]

#-------------------------------------------------------------------------------
# Construct a porosity model
#-------------------------------------------------------------------------------
# First we need some conversion factors from CGS to our units
cgs = CGS()
mCGSconv = cgs.unitMassKg/units.unitMassKg
lCGSconv = cgs.unitLengthMeters/units.unitLengthMeters
tCGSconv = cgs.unitTimeSec/units.unitTimeSec
PCGSconv = mCGSconv/(tCGSconv*tCGSconv*lCGSconv)
vCGSconv = lCGSconv/tCGSconv

if PorousModel is PalphaPorosity:
    alpha0 = 1.275
    Pe = 8e8      # dynes/cm^2
    Ps = 7e9      # dynes/cm^2
    #cS0 = 5.35e5  # cm/sec
    ce = 4.11e5   # cm/sec
    alphae = (alpha0 - 1.0)*((Ps - Pe)/(Ps - Ps0))**2 + 1.0
    print("P-alpha porosity model parameters:")
    output("  alpha0")
    output("      Pe")
    output("      Ps")
    output("     cS0")
    output("      ce")
    output("  alphae")
    print("Computing cS0 from solid EOS yields ", eosS.soundSpeed(rho0, 0.0))
    porosityAl = PalphaPorosity(eos, strengthModel, nodes,
                                phi0 = 1.0 - 1.0/alpha0,
                                Pe = Pe * PCGSconv,
                                Pt = Pe * PCGSconv,
                                alphae = alphae,
                                alphat = alphae,
                                n1 = 0.0,
                                n2 = 2.0,
                                cS0 = cS0 * vCGSconv,
                                c0 = ce * vCGSconv)

output("porosityAl")
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
integrator.dtGrowth = dtGrowth
integrator.domainDecompositionIndependent = domainIndependent
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
# Add the boundary conditions.
#-------------------------------------------------------------------------------
for package in integrator.physicsPackages():
    for bc in bcs:
        package.appendBoundary(bc)

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            volumeType = volumeType,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle)
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
    from SpheralMatplotlib import *
    state = State(db, integrator.physicsPackages())
    H = state.symTensorFields("H")
    h = db.newFluidScalarFieldList(0.0, "h")
    for i in range(nodes.numInternalNodes):
        h[0][i] = 1.0/H[0][i].xx
    rhoPlot = plotFieldList(state.scalarFields("mass density"),
                            plotStyle = "o-",
                            winTitle = "rho @ %g %i" % (control.time(), mpi.procs))
    velPlot = plotFieldList(state.vectorFields("velocity"),
                            yFunction = "%s.x",
                            plotStyle = "o-",
                            winTitle = "vel @ %g %i" % (control.time(), mpi.procs))
    PPlot = plotFieldList(state.scalarFields("pressure"),
                          plotStyle = "o-",
                          winTitle = "pressure @ %g %i" % (control.time(), mpi.procs))
    uPlot = plotFieldList(state.scalarFields("specific thermal energy"),
                          plotStyle = "o-",
                          winTitle = "specific thermal energy @ %g %i" % (control.time(), mpi.procs))
    SPlot = plotFieldList(state.symTensorFields("deviatoric stress"),
                          yFunction = "%s.xx",
                          plotStyle = "o-",
                          winTitle = "$S_{xx}$ @ %g %i" % (control.time(), mpi.procs))
    hPlot = plotFieldList(h,
                          plotStyle = "o-",
                          winTitle = "h @ %g %i" % (control.time(), mpi.procs))
    csPlot = plotFieldList(state.scalarFields(HydroFieldNames.soundSpeed),
                           plotStyle = "o-",
                           xlabel = r"$x$",
                           ylabel = r"$c_S$",
                           winTitle = "Sound speed @ %g %i" %  (control.time(), mpi.procs))
    alpha = porosityAl.alpha
    alphaPlot = plotField(alpha,
                          plotStyle = "o-",
                          winTitle = r"$\alpha$ @ %g %i" %  (control.time(), mpi.procs))
    DalphaDt = porosityAl.DalphaDt
    DalphaDtPlot = plotField(DalphaDt,
                             plotStyle = "o-",
                             xlabel = r"$x$",
                             ylabel = r"$D\alpha/Dt$",
                             winTitle = r"$D\alpha/Dt$ @ %g %i" %  (control.time(), mpi.procs))
    phi = porosityAl.phi()
    phiPlot = plotField(phi,
                        plotStyle = "o-",
                        winTitle = r"$\phi$ @ %g %i" %  (control.time(), mpi.procs))

    dPdRplot = plotField(porosityAl.partialPpartialRho,
                         plotStyle = "o-",
                         winTitle = r"$\partial P/\partial \rho$ @ %g %i" %  (control.time(), mpi.procs))
    dPdUplot = plotField(porosityAl.partialPpartialEps,
                         plotStyle = "o-",
                         winTitle = r"$\partial P/\partial \varepsilon$ @ %g %i" %  (control.time(), mpi.procs))

    plots = [(rhoPlot, "rho.png"),
             (velPlot, "vel.png"),
             (PPlot, "pressure.png"),
             (SPlot, "devstress.png"),
             (uPlot, "u.png"),
             (hPlot, "h.png"),
             (alphaPlot, "alpha.png"),
             (phiPlot, "phi.png")]

    # Save the figures.
    for p, fname in plots:
        savefig(p, os.path.join(dataDir, fname))

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
    xprof = mpi.reduce([x.x for x in internalValues(pos)], mpi.SUM)
    rhoprof = mpi.reduce(internalValues(rho), mpi.SUM)
    Pprof = mpi.reduce(internalValues(P), mpi.SUM)
    vprof = mpi.reduce([v.x for v in internalValues(vel)], mpi.SUM)
    epsprof = mpi.reduce(internalValues(eps), mpi.SUM)
    hprof = mpi.reduce([1.0/sqrt(H.Determinant()) for H in internalValues(Hfield)], mpi.SUM)
    sprof = mpi.reduce([x.xx for x in internalValues(S)], mpi.SUM)
    phiprof = mpi.reduce(phi.internalValues(), mpi.SUM)
    mof = mortonOrderIndices(db)
    mo = mpi.reduce(internalValues(mof), mpi.SUM)
    if mpi.rank == 0:
        multiSort(mo, xprof, rhoprof, Pprof, vprof, epsprof, hprof, sprof, phiprof)
        with open(outputFile, "w") as f:
            f.write(("#" + 8*" %16s" + "\n") % ("x", "rho", "P", "v", "eps", "h", "S", "phi"))
            for (xi, rhoi, Pi, vi, epsi, hi, si, phii) in zip(xprof, rhoprof, Pprof, vprof, epsprof, hprof, sprof, phiprof):
                f.write((8*"%16.12e " + "\n") %
                        (xi, rhoi, Pi, vi, epsi, hi, si, phii))

        #---------------------------------------------------------------------------
        # Check the floating values for the state against reference data.
        #---------------------------------------------------------------------------
        import filearraycmp as fcomp
        assert fcomp.filearraycmp(outputFile, referenceFile, testtol, testtol)
        print("Floating point comparison test passed.")

        #---------------------------------------------------------------------------
        # Also we can optionally compare the current results with another file for
        # bit level consistency.
        #---------------------------------------------------------------------------
        if comparisonFile != "None" and BuildData.cxx_compiler_id != "IntelLLVM":
            comparisonFile = os.path.join(dataDir, comparisonFile)
            import filecmp
            assert filecmp.cmp(outputFile, comparisonFile)

if graphics:
    plt.show()
