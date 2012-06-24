#ATS:test(SELF, "--graphics False", label="Planar Hourglass test problem -- 1-D (serial)")
#-------------------------------------------------------------------------------
# A made up 1-D problem to test the anti-hourglassing algorithms.
#-------------------------------------------------------------------------------
from Spheral import *
from SpheralTestUtilities import *

title("1-D planar hourglassing test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(NodeListConstructor = SphNodeList1d,

            nx1 = 100,
            rho1 = 1.0,
            eps1 = 0.0,
            x0 = 0.0,
            x1 = 1.0,
            nPerh = 2.01,
            
            wavelength = 0.05,
            amplitude = 0.25,

            a0 = Vector1d(1.0),

            gamma = 5.0/3.0,
            mu = 1.0,
            Qconstructor = MonaghanGingoldViscosity1d,
            #Qconstructor = TensorMonaghanGingoldViscosity1d,
            Cl = 1.0,
            Cq = 1.0,
            Qlimiter = False,
            epsilon2 = 1e-2,
            negligibleSoundSpeed = 1e-5,
            csMultiplier = 1e-4,
            energyMultiplier = 0.1,
            HsmoothMin = 0.0001, 
            HsmoothMax = 100.0,
            HsmoothFraction = 0.0,
            cfl = 0.5,
            XSPH = True,
            epsilonTensile = 0.0,
            nTensile = 8,

            hourglassMultiplier = 0.001,

            neighborSearchType = Neighbor1d.NeighborSearchType.GatherScatter,
            numGridLevels = 20,
            topGridCellSize = 200.0,
            origin = Vector1d(0.0),

            goalTime = 1.0,
            steps = None,
            dt = 0.0001,
            dtMin = 1.0e-5, 
            dtMax = None,
            dtGrowth = 2.0,
            maxSteps = None,
            statsStep = 1,
            smoothIters = 0,
            HEvolution = Hydro1d.HEvolutionType.IdealH,
            sumForMassDensity = Hydro1d.MassDensityType.RigorousSumDensity, # VolumeScaledDensity,

            restoreCycle = None,
            restartStep = 10000,
            restartBaseName = "Hourglass-1d",

            graphics = "gnu",
            )

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS1d(gamma, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel1d(BSplineKernel1d(), 100)
kernelExtent = WT.kernelExtent()
WTPi = TableKernel1d(BSplineKernel1d(), 100)
#WTPi = TableKernel1d(HatKernel1d(kernelExtent, kernelExtent), 100)
output("WT")
output("WTPi")

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
nodes1 = NodeListConstructor("nodes1", eos, WT, WTPi)
nodes1.HsmoothFraction = HsmoothFraction
nodes1.nodesPerSmoothingScale = nPerh
output("nodes1.HsmoothFraction")
output("nodes1.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Construct the neighbor object.
#-------------------------------------------------------------------------------
neighbor1 = NestedGridNeighbor1d(nodes1,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes1.registerNeighbor(neighbor1)

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
from DistributeNodes import distributeNodesInRange1d
distributeNodesInRange1d([(nodes1, nx1, rho1, (x0, x1))])
nNodesThisDomain1 = nodes1.numInternalNodes
output("nodes1.numNodes")

# Set node specific thermal energies
nodes1.specificThermalEnergy(ScalarField1d("tmp", nodes1, eps1))

# Displace the nodes in a pattern that looks like the tensile instability clumping.
dx = (x1 - x0)/nx1
for i in xrange(nodes1.numInternalNodes):
    delta = amplitude*((-1.0)**(i % 2))*dx # amplitude*sin(2.0*pi*nodes1.positions()[i].x/wavelength)
    nodes1.positions()[i].x += delta

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase1d()
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
hydro = Hydro1d(WT, WTPi, q)
hydro.cfl = cfl
hydro.HEvolution = HEvolution
hydro.sumForMassDensity = sumForMassDensity
hydro.HsmoothMin = HsmoothMin
hydro.HsmoothMax = HsmoothMax
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
accel = ConstantAcceleration1d(a0, nodes1, indicies)

#-------------------------------------------------------------------------------
# Construct an hour glass control object.
#-------------------------------------------------------------------------------
hourglass = ThirdMomentHourglassControl1d(db, hourglassMultiplier)
output("hourglass")
output("hourglass.multiplier")

packages = [hydro, accel, hourglass]

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane0 = Plane1d(Vector1d(x0), Vector1d( 1.0))
xPlane1 = Plane1d(Vector1d(x1), Vector1d(-1.0))
xbc0 = PeriodicBoundary1d(xPlane0, xPlane1)
for p in packages:
    p.appendBoundary(xbc0)

#-------------------------------------------------------------------------------
# Construct a predictor corrector integrator, and add the one physics package.
#-------------------------------------------------------------------------------
integrator = CheapSynchronousRK2Integrator1d(db)
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
                            initializeMassDensity = False)
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
# Compute the third moment of the given nodes.
# Since this is 1-D we just construct a scalar fields with the (0,0,0) component
# of the tensor.
#-------------------------------------------------------------------------------
def sgn(x):
    if x < 0:
        return -1
    else:
        return 1
    
def thirdMomentField(nodes):
    db.updateConnectivityMap()
    cm = db.connectivityMap()
    result = ScalarField1d("third moment", nodes)
    pos = nodes.positions().internalValues()
    H = nodes.Hfield().internalValues()
    for i in xrange(nodes.numInternalNodes):
        xi = pos[i].x
        Hi = H[i].xx
        neighbors = cm.connectivityForNode(nodes, i)[0]
        for j in neighbors:
            xij = xi - pos[j].x
            etai = abs(Hi*xij)
            Wi = WT.kernelValue(etai, 1.0)
            #Wi = abs(WT.gradValue(etai, 1.0))
            result[i] += sgn(xij) * Wi**3
    return result

def ggThirdMomentField(nodes, thirdMoment):
    db.updateConnectivityMap()
    cm = db.connectivityMap()
    result = ScalarField1d("grad grad third moment", nodes)
    mass = nodes.mass().internalValues()
    rho = nodes.massDensity().internalValues()
    pos = nodes.positions().internalValues()
    H = nodes.Hfield().internalValues()
    T = thirdMoment.internalValues()
    for i in xrange(nodes.numInternalNodes):
        xi = pos[i]
        Hi = H[i].xx
        hi2 = 1.0/Hi**2
        Ti = T[i]
        neighbors = cm.connectivityForNode(nodes, i)[0]
        for j in neighbors:
            Tji = T[j] - Ti
            xij = xi - pos[j]
            etai = Hi*xij
            Hetai = Hi*etai.unitVector()
            gradWi = (Hetai * WT.gradValue(abs(etai.x), Hi)).x
            xij = xij.x
            result[i] += mass[j] * Tji*xij/(xij*xij + 0.01*hi2) * gradWi
        result[i] /= rho[i]
    return result
    

#-------------------------------------------------------------------------------
# Plot the final state.
#-------------------------------------------------------------------------------
import Gnuplot
from SpheralGnuPlotUtilities import *
## rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db)
## Eplot = plotEHistory(control.conserve)
## xplot = plotFieldList(db.fluidPosition,
##                       yFunction = "%s.x",
##                       plotStyle = "points",
##                       winTitle = "Position (x)")
## aplot = plotFieldList(hourglass.acceleration(),
##                       yFunction = "%s.x",
##                       plotStyle = "points",
##                       winTitle = "Hourglass accleration")

#tm = thirdMomentField(nodes1)
#tml = ScalarFieldList1d()
#tml.appendField(tm)
tml = hourglass.thirdMoment()
tmplot = plotFieldList(tml,
                       yFunction = "%s(0,0,0)",
                       plotStyle = "linespoints",
                       winTitle = "Third Moment")

a = db.fluidDvelocityDt
aplot = plotFieldList(a,
                      yFunction = "%s.x",
                      plotStyle = "linespoints",
                      winTitle = "Acceleration")

## ggtm = ggThirdMomentField(nodes1, tm)
## ggtml = ScalarFieldList1d()
## ggtml.appendField(ggtm)
## ggtmplot = plotFieldList(ggtml,
##                          plotStyle = "linespoints",
##                          winTitle = "Grad^2 Third Moment")
