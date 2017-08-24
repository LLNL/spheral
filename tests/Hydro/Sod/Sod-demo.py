import Spheral
import SpheralGnuPlotUtilities
import DistributeNodes

from SodAnalyticSolution import *

#-------------------------------------------------------------------------------
# Create a gamma law gas equation of state.
#-------------------------------------------------------------------------------
gamma, mu = 1.4, 1.0
eos = Spheral.GammaLawGasMKS1d(gamma, mu)

#-------------------------------------------------------------------------------
# Choose an interpolation kernel.
#-------------------------------------------------------------------------------
WT = Spheral.TableKernel1d(Spheral.BSplineKernel1d(), 100)
WTPi = WT

#-------------------------------------------------------------------------------
# Create individual SPH node lists to represent the high and low pressure
# regions, each with 100 nodes.
#-------------------------------------------------------------------------------
hmin = 0.0001
hmax = 0.1
nodes1 = Spheral.SphNodeList1d("heavy material", eos, WT, WTPi)
nodes2 = Spheral.SphNodeList1d("light material", eos, WT, WTPi)
for nodes in [nodes1, nodes2]:
    nodes.hmin = hmin
    nodes.hmax = hmax
del nodes

#-------------------------------------------------------------------------------
# Construct the neighbor objects and associate them with the node lists.
#-------------------------------------------------------------------------------
neighborSearchType = Spheral.Neighbor1d.NeighborSearchType.GatherScatter
numGridLevels = 10
topGridCellSize = 0.25
origin = Spheral.Vector1d(0.0)
neighbor1 = Spheral.NestedGridNeighbor1d(nodes1,
                                         neighborSearchType,
                                         numGridLevels,
                                         topGridCellSize,
                                         origin,
                                         WT.kernelExtent())
nodes1.registerNeighbor(neighbor1)
neighbor2 = Spheral.NestedGridNeighbor1d(nodes2,
                                         neighborSearchType,
                                         numGridLevels,
                                         topGridCellSize,
                                         origin,
                                         WT.kernelExtent())
nodes2.registerNeighbor(neighbor2)

#-------------------------------------------------------------------------------
# Initialize the node positions
#-------------------------------------------------------------------------------
nx1, nx2 = 400, 100
rho1, rho2 = 1.0, 0.25
x0, x1, x2 = -0.5, 0.0, 0.5
DistributeNodes.distributeNodesInRange1d([(nodes1, nx1, rho1, (x0, x1)),
                                          (nodes2, nx2, rho2, (x1, x2))])

# Set node specific thermal energies
P1, P2 = 1.0, 0.1795
eps1 = P1/((gamma - 1.0)*rho1)
eps2 = P2/((gamma - 1.0)*rho2)
nodes1.specificThermalEnergy(Spheral.ScalarField1d("tmp", nodes1, eps1))
nodes2.specificThermalEnergy(Spheral.ScalarField1d("tmp", nodes2, eps2))

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = Spheral.DataBase1d()
db.appendNodeList(nodes1)
db.appendNodeList(nodes2)

#-------------------------------------------------------------------------------
# Construct a standard Monaghan-Gingold artificial viscosity.
#-------------------------------------------------------------------------------
qMG = Spheral.MonaghanGingoldViscosity1d(1.0, 1.5)

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
hydro = Spheral.Hydro1d(WT, WTPi, qMG)
hydro.cfl = 0.5
hydro.HsmoothMin = hmin
hydro.HsmoothMax = hmax
hydro.sumForMassDensity = Spheral.Hydro1d.MassDensityType.RigorousSumDensity

#-------------------------------------------------------------------------------
# Create boundary conditions, and give them to the hydro.
#-------------------------------------------------------------------------------
xPlane0 = Spheral.Plane1d(Spheral.Vector1d(x0), Spheral.Vector1d(1.0))
xPlane1 = Spheral.Plane1d(Spheral.Vector1d(x2), Spheral.Vector1d(-1.0))
xbc0 = Spheral.ReflectingBoundary1d(xPlane0)
xbc1 = Spheral.ReflectingBoundary1d(xPlane1)
hydro.appendBoundary(xbc0)
hydro.appendBoundary(xbc1)

#-------------------------------------------------------------------------------
# Construct a predictor-corrector integrator, and add the one physics package.
#-------------------------------------------------------------------------------
integrator = Spheral.PredictorCorrectorIntegrator1d(db)
integrator.appendPhysicsPackage(hydro)
integrator.lastDt = 1e-4
integrator.dtMin = 1e-5
integrator.dtMax = 0.1
integrator.dtGrowth = 2.0

#-------------------------------------------------------------------------------
# Construct a controller to manage the simulation for us.
#-------------------------------------------------------------------------------
restartStep = 1000
restartBaseName = "Sod-planar-1d-%i" % (nx1 + nx2)
control = Spheral.SpheralController(integrator, WT,
                                    restartStep = restartStep,
                                    restartBaseName = restartBaseName)
control.iterateIdealH()
db.updateFluidMassDensity(WT, db.fluidMassDensity, db.fluidOmegaGradh)

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
goalTime = 0.15
control.advance(goalTime)

#-------------------------------------------------------------------------------
# Plot the final state.
#-------------------------------------------------------------------------------
rhoPlot, velPlot, epsPlot, PPlot, HPlot = SpheralGnuPlotUtilities.plotState(db)

#-------------------------------------------------------------------------------
# Now overplot the analytic solution.
#-------------------------------------------------------------------------------
dx1 = (x1 - x0)/nx1
dx2 = (x2 - x1)/nx2
h1 = nodes1.nodesPerSmoothingScale*dx1
h2 = nodes2.nodesPerSmoothingScale*dx2
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
SpheralGnuPlotUtilities.plotAnswer(answer, control.time(),
                                   rhoPlot, velPlot, epsPlot, PPlot, HPlot)
