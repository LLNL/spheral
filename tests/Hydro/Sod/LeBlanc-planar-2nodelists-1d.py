from Numeric import *
from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from SodAnalyticSolution import *

import loadmpi
mpi, rank, numprocs = loadmpi.loadmpi()

title("1-D integrated hydro test -- planar LeBlanc shock tube problem")

#-------------------------------------------------------------------------------
# A function to generate our geometric progression of node properties, designed
# to match the requested constant density.
#-------------------------------------------------------------------------------
def geometricNodeDistribution(nx, x0, x1, rho, hmultiplier, m0, m1):
    f = (m1/m0)**(1.0/nx) - 1
    integral = 0.0
    for i in xrange(nx):
        integral = integral + (1.0 + f)**i
    dx0 = (x1 - x0)/integral
    print "x0, x1", x0, x0 + integral*dx0

    x = [x0]
    m = [m0]
    h = [1.0/(hmultiplier*dx0)]
    for i in xrange(1, nx):
        facprev = (1.0 + f)**(i - 1)
        dxprev = dx0*facprev
        fac = (1.0 + f)**i
        dx = dx0*fac
        x.append(x[i - 1] + 0.5*(dxprev + dx))
        m.append(fac*m0)
        h.append(1.0/(hmultiplier*dx))

    return x, m, h

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
nx1, nx2 = 100, 100
rho1, rho2 = 1.0, 0.001
eps1, eps2 = 0.1, 1e-7
x0, x1, x2 = 0.0, 3.0, 10.0

nPerh1, nPerh2 = 2.01, 2.01

m1 = (x1 - x0)*rho1/nx1
m2 = (x2 - x1)*rho2/nx2

gamma = 5.0/3.0
mu = 1.0
Cl = 0.75
Cq = 2.0
epsilon2 = 1e-8
HsmoothMin, HsmoothMax = 0.001, 1.0
cfl = 0.5

P1 = (gamma - 1.0)*rho1*eps1
P2 = (gamma - 1.0)*rho2*eps2

neighborSearchType = Neighbor1d.NeighborSearchType.GatherScatter
numGridLevels = 15
topGridCellSize = 2.5
origin = Vector1d(0.0)

goalTime = 6.0
dt = 1e-4
dtMin, dtMax = 1.0e-5, 0.1
dtGrowth = 2.0
maxSteps = None
statsStep = 10
smoothIters = 0
HEvolution = Hydro1d.HEvolutionType.IdealH # IntegrateH
sumForMassDensity = Hydro1d.MassDensityType.RigorousSumDensity

restoreCycle = None
restartStep = 1000
restartBaseName = "LeBlanc-planar-1d-%i" % (nx1 + nx2)

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS1d(gamma, mu)

#-------------------------------------------------------------------------------
# Make the NodeLists.
#-------------------------------------------------------------------------------
nodes1 = SphNodeList1d("heavy nodes", eos)
nodes1.nodesPerSmoothingScale = nPerh1
output("nodes1.nodesPerSmoothingScale")

nodes2 = SphNodeList1d("light nodes", eos)
nodes2.nodesPerSmoothingScale = nPerh2
output("nodes2.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
# Set node positions
from DistributeNodes import distributeNodesInRange1d
distributeNodesInRange1d([(nodes1, nx1, rho1, (x0, x1)),
                          (nodes2, nx2, rho2, (x1, x2))])
output("nodes1.numNodes")
output("nodes2.numNodes")

# Set node specific thermal energies
nodes1.setSpecificThermalEnergy(ScalarField1d("tmp", nodes1, eps1))
nodes2.setSpecificThermalEnergy(ScalarField1d("tmp", nodes2, eps2))

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel1d(BSplineKernel1d(), 100)
WTPi = WT
output("WT")
output("WTPi")
kernelExtent = WT.kernelExtent()

#-------------------------------------------------------------------------------
# Construct the neighbor objects.
#-------------------------------------------------------------------------------
neighbor1 = NestedGridNeighbor1d(nodes1,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
neighbor2 = NestedGridNeighbor1d(nodes2,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes1.registerNeighbor(neighbor1)
nodes2.registerNeighbor(neighbor2)

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase1d()
output("db")
output("db.appendNodeList(nodes1)")
output("db.appendNodeList(nodes2)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct a standard Monaghan-Gingold artificial viscosity.
#-------------------------------------------------------------------------------
q = MonaghanGingoldViscosity1d(Cl, Cq)
q.epsilon2 = epsilon2
output("q")
output("q.epsilon2")

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
output("hydro.kernel")
output("hydro.PiKernel")
output("hydro.cfl")
output("hydro.HEvolution")
output("hydro.sumForMassDensity")
output("hydro.HsmoothMin")
output("hydro.HsmoothMax")
output("hydro.valid()")

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane0 = Plane1d(Vector1d(x0), Vector1d( 1.0))
xPlane1 = Plane1d(Vector1d(x2), Vector1d(-1.0))
xbc0 = ReflectingBoundary1d(xPlane0)
xbc1 = ReflectingBoundary1d(xPlane1)
hydro.appendBoundary(xbc0)
hydro.appendBoundary(xbc1)

#-------------------------------------------------------------------------------
# Construct a predictor corrector integrator, and add the one physics package.
#-------------------------------------------------------------------------------
integrator = PredictorCorrectorIntegrator1d(db)
integrator.appendPhysicsPackage(hydro)
integrator.lastDt = dt
integrator.dtGrowth = dtGrowth
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.lastDt")
output("integrator.dtMax")
output("integrator.dtMin")
output("integrator.valid()")

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName)
output("control")

if restoreCycle is not None:
    control.loadRestartFile(restoreCycle)
else:
    control.smoothState(smoothIters)

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if control.time() < goalTime:
    control.step(5)
    control.advance(goalTime, maxSteps)
    control.dropRestartFile()

#-------------------------------------------------------------------------------
# Plot the final state.
#-------------------------------------------------------------------------------
rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db)

# Reset the ranges for the density and pressure fields, and make their plots
# log plots.
for plot, fieldList in [(rhoPlot, db.fluidMassDensity),
                        (epsPlot, db.fluidSpecificThermalEnergy),
                        (PPlot, db.fluidPressure)]:
    minY = 1e30
    maxY = -1e30
    for field in fieldList.fields():
        if len(field.internalValues()) > 0:
            minY = min(minY, min(field.internalValues()))
            maxY = max(maxY, max(field.internalValues()))
        if mpi:
            globalMinY = mpi.reduce(minY, mpi.MIN)
            globalMaxY = mpi.reduce(maxY, mpi.MAX)
            if mpi.rank == 0:
                minY = globalMinY
                maxY = globalMaxY
    minY = minY/2.0
    maxY = maxY*2.0
    if plot:
        plot("set yrange[%g:%g]" % (minY, maxY))
        plot("set logscale y")
        plot.replot()

# Now overplot the analytic solution.
dx1 = (x1 - x0)/nx1
dx2 = (x2 - x1)/nx2
h1 = nPerh1*dx1
h2 = nPerh2*dx2
answer = SodSolution(nPoints=nx1 + nx2,
                     gamma = gamma,
                     rho1 = rho1,
                     P1 = P1,
                     rho2 = rho2,
                     P2 = P2,
                     x0 = x0,
                     x1 = x1,
                     x2 = x2,
                     h1 = 1.0/h1,
                     h2 = 1.0/h2)

plotAnswer(answer, control.time(),
           rhoPlot, velPlot, epsPlot, PPlot, HPlot)
