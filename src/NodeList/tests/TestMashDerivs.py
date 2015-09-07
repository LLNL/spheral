from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *

################################################################################
# Generic problem parameters
n1 = 100
rho1 = 1.0
x0, x1 = -0.5, 0.5

m1 = (x1 - x0)*rho1/n1

eps0, eps1 = 0.0, 1.0
v0, v1 = 1.0, 0.0

gamma = 2.0
mu = 1.0
Cl = 0.75
Cq = 2.0
epsilon2 = 1e-4
HsmoothMin, HsmoothMax = 0.0001, 0.1

neighborSearchType = 3 # GatherScatter
#neighborSearchType = 2 # Scatter
#neighborSearchType = 1 # Gather
numGridLevels = 10
topGridCellSize = 0.25
origin = Vector1d(0.0)

sumForMassDensity = 2

################################################################################
title('1-D MASH derivatives test.')

eos = GammaLawGasMKS1d(gamma, mu)

# Create two empty NodeLists
nodes1 = MashNodeList1d(n1, eos)
output('nodes1.numNodes')

# Set node positions.
dx = (x1 - x0)/n1
for i in xrange(n1):
    nodes1.positions[i] = x0 + i*dx

# Set node masses
nodes1.mass[:] = [m1]*n1

# Set node specific thermal energies
deps = (eps1 - eps0)/n1
for i in xrange(n1):
    nodes1.specificThermalEnergy[i] = eps0 + i*deps

# Set node densities.
nodes1.massDensity[:] = [rho1]*n1

# Set node velocities
dv = (v1 - v0)/n1
for i in xrange(n1):
    nodes1.velocity[i].x = v0 + i*dv

# Set node weights.
nodes1.weight[:] = [1.0]*n1

# Set the smoothing scales.
dx1 = (x1 - x0)/n1
h1 = 1.0/(2.01*dx1)
for H in nodes1.Hfield:
    H.xx = h1

# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
W = BSplineKernel1d()
#W = W4SplineKernel1d()
#W = GaussianKernel1d()
#W = SuperGaussianKernel1d()
#WPi = PiGaussianKernel1d(1.0)
WPi = W
output('W')
output('WPi')
kernelExtent = W.kernelExtent

# Set the table kernel.
WT = TableKernel1d()
WT.setTableData(W, 100)
output('WT')
WTPi = TableKernel1d()
WTPi.setTableData(WPi, 100)
output('WTPi')

# Construct the neighbor object and associate it with the node list.
neighborTimer = SpheralTimer('Neighbor initialization.')
neighborTimer.start()
neighbor1 = NestedGridNeighbor1d(nodes1,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes1.neighbor = neighbor1
neighborTimer.stop()
neighborTimer.printStatus()

# Create boundary conditions.  We need at least this much to create the initial
# mass density field.
xPlane0 = Plane1d((x0), ( 1.0))
xPlane1 = Plane1d((x1), (-1.0))
xbc0 = ReflectingBoundary1d(xPlane0)
xbc1 = ReflectingBoundary1d(xPlane1)
bcList = [xbc0, xbc1]

# Construct a DataBase to hold our node list
db = DataBase1d()
output('db')
output('db.appendNodeList(nodes1)')
output('db.numNodeLists')
output('db.numFluidNodeLists')

# Construct a standard Monaghan-Gingold artificial viscosity.
q = MonaghanGingoldViscosity1d(Cl, Cq)
q.epsilon2 = epsilon2
# q = VonNeumanViscosity1d(Cl, Cq)
output('q')
output('q.Cl')
output('q.Cq')

# Construct the hydro physics object.
hydro = Hydro1d(q)
hydro.kernel = WT
hydro.PiKernel = WTPi
output('hydro')
output('hydro.kernel')
output('hydro.PiKernel')
output('hydro.cfl')
output('hydro.valid')

# Construct a synchronous RK2 integrator, and add the one physics package and
# boundary condtion.
integrator = CheapSynchronousRK2Integrator1d(db)
output('integrator')
integrator.appendPhysicsPackage(hydro)
output('integrator.havePhysicsPackage(hydro)')
output('integrator.valid')
integrator.HsmoothMin = HsmoothMin
integrator.HsmoothMax = HsmoothMax
output('integrator.HsmoothMin')
output('integrator.HsmoothMax')
integrator.kernel = WT
output('integrator.kernel')
integrator.sumForMassDensity = sumForMassDensity
output('integrator.sumForMassDensity')

control = SpheralController(integrator, WT,  boundaryConditions=[xbc0, xbc1])
output('control')

##################################################################################

# Enforce boundaries on the specific thermal energy and velocities.
depsdx = (eps1 - eps0)/(x1 - x0)
dvdx = (v1 - v0)/(x1 - x0)
for i in xrange(nodes1.firstGhostNode, nodes1.numNodes):
    x = nodes1.positions[i].x
    nodes1.specificThermalEnergy[i] = eps0 + (x - x0)*depsdx
    nodes1.velocity[i].x = v0 + (x - x0)*dvdx

# Determine the hydrodynamic derivatives.
rhoSum = ScalarFieldList1d()
maxQ = ScalarFieldList1d()
DrhoDt = ScalarFieldList1d()
DvDt = VectorFieldList1d()
DepsDt = ScalarFieldList1d()
DvDx = TensorFieldList1d()
DHDt = SymTensorFieldList1d()

rhoSum.appendField(nodes1.massDensitySum)
maxQ.appendField(nodes1.maxViscousPressure)
DrhoDt.appendField(nodes1.DmassDensityDt)
DvDt.appendField(nodes1.DvelocityDt)
DepsDt.appendField(nodes1.DspecificThermalEnergyDt)
DvDx.appendField(nodes1.DvelocityDx)
DHDt.appendField(nodes1.DHDt)

hydro.evaluateDerivatives(db, 1.0, 1.0, rhoSum, maxQ,
                          DrhoDt, DvDt, DepsDt, DvDx, DHDt)

# Plot the mass derivative.
DrhoPlot = plotFieldList(DrhoDt,
                         winTitle='Mass density derivative')
DvelPlot = plotFieldList(DvDt,
                         yFunction='%s.x',
                         winTitle='Velocity derivative')
DepsPlot = plotFieldList(DepsDt,
                         winTitle='Specific thermal energy derivative')
