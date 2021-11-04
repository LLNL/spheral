#ATS:test(SELF, "--graphics False", label="Hourglass test problem -- 2-D (serial)")
#-------------------------------------------------------------------------------
# A made up 2-D problem to test the anti-hourglassing algorithms.
#-------------------------------------------------------------------------------
from Spheral import *
from SpheralTestUtilities import *

title("2-D hourglassing test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(NodeListConstructor = SphNodeList2d,

            nx = 20,
            ny = 20,
            rho1 = 1.0,
            eps1 = 1.0,
            xmin = (0.0, 0.0),
            xmax = (1.0, 1.0),
            nPerh = 2.01,
            
            wavelength = 0.2,
            amplitude = 0.025,

            a0 = Vector2d(1.0, 0.0),

            gamma = 5.0/3.0,
            mu = 1.0,
            Qconstructor = MonaghanGingoldViscosity2d,
            #Qconstructor = TensorMonaghanGingoldViscosity2d,
            Cl = 0.5,
            Cq = 0.5,
            Qlimiter = False,
            epsilon2 = 1e-4,
            negligibleSoundSpeed = 1e-5,
            csMultiplier = 1e-4,
            energyMultiplier = 0.1,
            hmin = 0.0001, 
            hmax = 0.5,
            HsmoothFraction = 0.0,
            cfl = 0.5,
            XSPH = True,
            epsilonTensile = 0.0,
            nTensile = 8,

            hourglassMultiplier = 0.1,
            hourglassAccelerationFactor = 0.01,

            neighborSearchType = Neighbor2d.NeighborSearchType.GatherScatter,
            numGridLevels = 20,
            topGridCellSize = 2.0,
            origin = Vector2d(0.0, 0.0),

            goalTime = 1.0,
            dt = 0.0001,
            dtMin = 1.0e-5, 
            dtMax = None,
            dtGrowth = 2.0,
            maxSteps = None,
            statsStep = 1,
            smoothIters = 0,
            HEvolution = Hydro2d.HEvolutionType.IdealH,
            sumForMassDensity = Hydro2d.MassDensityType.RigorousSumDensity, # VolumeScaledDensity,

            restoreCycle = None,
            restartStep = 10000,
            restartBaseName = "Hourglass-2d",

            graphics = "gnu",
            )

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS2d(gamma, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel2d(BSplineKernel2d(), 100)
kernelExtent = WT.kernelExtent()
WTPi = TableKernel2d(BSplineKernel2d(), 100)
#WTPi = TableKernel2d(HatKernel2d(kernelExtent, kernelExtent), 100)
output("WT")
output("WTPi")

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
nodes1 = NodeListConstructor("nodes1", eos, WT, WTPi)
nodes1.HsmoothFraction = HsmoothFraction
nodes1.nodesPerSmoothingScale = nPerh
nodes1.hmin = hmin
nodes1.hmax = hmax
output("nodes1.HsmoothFraction")
output("nodes1.nodesPerSmoothingScale")
output("nodes1.hmin")
output("nodes1.hmax")

#-------------------------------------------------------------------------------
# Construct the neighbor object.
#-------------------------------------------------------------------------------
neighbor1 = NestedGridNeighbor2d(nodes1,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes1.registerNeighbor(neighbor1)

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
from DistributeNodes import distributeNodes2d
from GenerateNodeDistribution2d import GenerateNodeDistribution2d
generator = GenerateNodeDistribution2d(nx, ny, rho1, "lattice",
                                       xmin = xmin,
                                       xmax = xmax,
                                       nNodePerh = nPerh)
distributeNodes2d((nodes1, generator))
output("nodes1.numNodes")

# Set node specific thermal energies
nodes1.specificThermalEnergy(ScalarField2d("tmp", nodes1, eps1))

# Displace the nodes in a pattern that looks like hourglassing.
for i in xrange(nodes1.numInternalNodes):
    dx = amplitude*sin(2.0*pi*nodes1.positions()[i].x/wavelength)
    nodes1.positions()[i].x += dx

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase2d()
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
# Set the XSPH and tensile corrections for the NodeList
#-------------------------------------------------------------------------------
nodes1.XSPH = XSPH
nodes1.epsilonTensile = epsilonTensile
nodes1.nTensile = nTensile
output("nodes1.XSPH")
output("nodes1.epsilonTensile")
output("nodes1.nTensile")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
hydro = Hydro2d(WT, WTPi, q)
hydro.cfl = cfl
hydro.HEvolution = HEvolution
hydro.sumForMassDensity = sumForMassDensity
hydro.HsmoothMin = hmin
hydro.HsmoothMax = hmax
output("hydro")
output("hydro.cfl")
output("hydro.HEvolution")
output("hydro.sumForMassDensity")
output("hydro.HsmoothMin")
output("hydro.HsmoothMax")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.valid()")

#-------------------------------------------------------------------------------
# Construct a constant acceleration package.
#-------------------------------------------------------------------------------
indicies = vector_of_int()
indicies.extend(range(nodes1.numInternalNodes))
accel = ConstantAcceleration2d(a0, nodes1, indicies)

#-------------------------------------------------------------------------------
# Construct an hour glass control object.
#-------------------------------------------------------------------------------
hourglass = SecondMomentHourglassControl2d(hourglassMultiplier,
                                           hourglassAccelerationFactor)
output("hourglass")
output("hourglass.multiplier")
output("hourglass.maxAccelerationFactor")

packages = [hydro, accel, hourglass]

#-------------------------------------------------------------------------------
# Boundary conditions.
#-------------------------------------------------------------------------------
xbc0 = ReflectingBoundary2d(Plane2d(Vector2d(*xmin), Vector2d(1.0, 0.0)))
xbc1 = ReflectingBoundary2d(Plane2d(Vector2d(*xmax), Vector2d(-1.0, 0.0)))
ybc0 = ReflectingBoundary2d(Plane2d(Vector2d(*xmin), Vector2d(0.0, 1.0)))
ybc1 = ReflectingBoundary2d(Plane2d(Vector2d(*xmax), Vector2d(0.0, -1.0)))

## for bc in [xbc0, xbc1, ybc0, ybc1]:
##     for p in packages:
##         p.appendBoundary(bc)

#-------------------------------------------------------------------------------
# Construct a predictor corrector integrator, and add the one physics package.
#-------------------------------------------------------------------------------
integrator = PredictorCorrectorIntegrator2d(db)
output("integrator")
for p in packages:
    integrator.appendPhysicsPackage(p)
output("integrator.valid()")
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
                            restartBaseName = restartBaseName,
                            initializeMassDensity = True)
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
if control.time() < goalTime:
    control.step(5)
    control.advance(goalTime, maxSteps)

#-------------------------------------------------------------------------------
# Plot the final state.
#-------------------------------------------------------------------------------
import Gnuplot
from SpheralGnuPlotUtilities import *
rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db)
Eplot = plotEHistory(control.conserve)
xplot = plotNodePositions2d(db,
                            title = "Positions")
