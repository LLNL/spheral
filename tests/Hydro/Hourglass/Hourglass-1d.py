#ATS:test(SELF, "--graphics False", label="Planar Hourglass test problem -- 1-D (serial)")
#-------------------------------------------------------------------------------
# A made up 1-D problem to test the anti-hourglassing algorithms.
#-------------------------------------------------------------------------------
from Spheral1d import *
from SpheralTestUtilities import *

title("1-D planar hourglassing test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(nx1 = 100,
            rho1 = 1.0,
            eps1 = 1.0,
            x0 = 0.0,
            x1 = 1.0,
            nPerh = 2.01,
            
            wavelength = 0.05,
            amplitude = 0.25,

            a0 = Vector(1.0),

            SVPH = False,
            CRKSPH = False,
            filter = 0.0,
            gamma = 5.0/3.0,
            mu = 1.0,
            Qconstructor = MonaghanGingoldViscosity,
            Cl = 1.0,
            Cq = 1.0,
            Qlimiter = False,
            epsilon2 = 1e-2,
            negligibleSoundSpeed = 1e-5,
            csMultiplier = 1e-4,
            energyMultiplier = 0.1,
            hmin = 0.0001, 
            hmax = 100.0,
            cfl = 0.5,
            XSPH = False,
            epsilonTensile = 0.0,
            nTensile = 8,

            goalTime = 1.0,
            steps = None,
            dt = 0.0001,
            dtMin = 1.0e-5, 
            dtMax = None,
            dtGrowth = 2.0,
            maxSteps = None,
            statsStep = 1,
            smoothIters = 0,
            HUpdate = IdealH,
            densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
            compatibleEnergy = True,
            gradhCorrection = False,
            domainIndependent = False,

            restoreCycle = None,
            restartStep = 10000,
            restartBaseName = "Hourglass-1d",

            graphics = True,
            )

#-------------------------------------------------------------------------------
# CRKSPH Switches to ensure consistency
#-------------------------------------------------------------------------------
if CRKSPH:
    Qconstructor = CRKSPHMonaghanGingoldViscosity

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS(gamma, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(NBSplineKernel(5), 1000)
Wfbase = NBSplineKernel(7)
WTf = TableKernel(Wfbase, 1000, hmult = 1.0/(nPerh*Wfbase.kernelExtent))
kernelExtent = WT.kernelExtent
output("WT")

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
nodes1 = makeFluidNodeList("nodes1", eos, 
                           hmin = hmin,
                           hmax = hmax,
                           nPerh = nPerh)
output("nodes1")
output("nodes1.hmin")
output("nodes1.hmax")
output("nodes1.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
from DistributeNodes import distributeNodesInRange1d
distributeNodesInRange1d([(nodes1, nx1, rho1, (x0, x1))])
nNodesThisDomain1 = nodes1.numInternalNodes
output("nodes1.numNodes")

# Set node specific thermal energies
nodes1.specificThermalEnergy(ScalarField("tmp", nodes1, eps1))

# Displace the nodes in a pattern that looks like the tensile instability clumping.
dx = (x1 - x0)/nx1
for i in xrange(nodes1.numInternalNodes):
    delta = amplitude*((-1.0)**(i % 2))*dx # amplitude*sin(2.0*pi*nodes1.positions()[i].x/wavelength)
    nodes1.positions()[i].x += delta

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
output("db.appendNodeList(nodes1)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the artificial viscosity.
#-------------------------------------------------------------------------------
q = Qconstructor(Cl, Cq)
q.limiter = Qlimiter
q.epsilon2 = epsilon2
q.negligibleSoundSpeed = negligibleSoundSpeed
q.csMultiplier = csMultiplier
q.energyMultiplier = energyMultiplier
output("q")
output("q.Cl")
output("q.Cq")
output("q.limiter")
output("q.epsilon2")
output("q.negligibleSoundSpeed")
output("q.csMultiplier")
output("q.energyMultiplier")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if SVPH:
    hydro = SVPHFacetedHydro(W = WT, 
                             Q = q,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             densityUpdate = densityUpdate,
                             XSVPH = XSPH,
                             linearConsistent = linearConsistent,
                             generateVoid = False,
                             HUpdate = HUpdate,
                             fcentroidal = fcentroidal,
                             fcellPressure = fcellPressure,
                             xmin = Vector(-100.0),
                             xmax = Vector( 100.0))
elif CRKSPH:
    hydro = CRKSPHHydro(W = WT, 
                        Q = q, 
                        filter = filter,
                        cfl = cfl,
                        compatibleEnergyEvolution = compatibleEnergy,
                        XSPH = XSPH,
                        densityUpdate = densityUpdate,
                        HUpdate = HUpdate)
else:
    hydro = SPHHydro(W = WT, 
                     Q = q,
                     cfl = cfl,
                     compatibleEnergyEvolution = compatibleEnergy,
                     gradhCorrection = gradhCorrection,
                     densityUpdate = densityUpdate,
                     HUpdate = HUpdate,
                     XSPH = XSPH,
                     epsTensile = epsilonTensile,
                     nTensile = nTensile)
output("hydro")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.densityUpdate")
output("hydro.HEvolution")

packages = [hydro]

# #-------------------------------------------------------------------------------
# # Construct a constant acceleration package.
# #-------------------------------------------------------------------------------
# indices = vector_of_int()
# for i in xrange(nodes1.numInternalNodes):
#     indices.append(i)
# accel = ConstantAcceleration(a0, nodes1, indices)
# packages.append(accel)

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane0 = Plane(Vector(x0), Vector( 1.0))
xPlane1 = Plane(Vector(x1), Vector(-1.0))
xbc0 = PeriodicBoundary(xPlane0, xPlane1)
for p in packages:
    p.appendBoundary(xbc0)

#-------------------------------------------------------------------------------
# Construct a predictor corrector integrator, and add the one physics package.
#-------------------------------------------------------------------------------
integrator = CheapSynchronousRK2Integrator(db)
output("integrator")
for p in packages:
    integrator.appendPhysicsPackage(p)
integrator.lastDt = dt
output("integrator.lastDt")
if dtMin:
    integrator.dtMin = dtMin
    output("integrator.dtMin")
if dtMax:
    integrator.dtMax = dtMax
    output("integrator.dtMax")
integrator.dtGrowth = dtGrowth
output("integrator.dtGrowth")

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName)
output("control")

# Smooth the initial conditions.
if restoreCycle is not None:
    control.loadRestartFile(restoreCycle)
else:
    control.iterateIdealH()
    control.smoothState(smoothIters)

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if steps is None:
    if control.time() < goalTime:
        control.step(5)
        control.advance(goalTime, maxSteps)
else:
    control.step(steps)

#-------------------------------------------------------------------------------
# Plot the final state.
#-------------------------------------------------------------------------------
import Gnuplot
from SpheralGnuPlotUtilities import *
rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db)
## Eplot = plotEHistory(control.conserve)
## xplot = plotFieldList(db.fluidPosition,
##                       yFunction = "%s.x",
##                       plotStyle = "points",
##                       winTitle = "Position (x)")

a = hydro.DvDt()
aplot = plotFieldList(a,
                      yFunction = "%s.x",
                      plotStyle = "linespoints",
                      winTitle = "Acceleration")

