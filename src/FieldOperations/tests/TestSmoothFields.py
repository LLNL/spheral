from Spheral import *
from SpheralTestUtilities import *


def plotFieldList(fieldList, plotGhost=0, color='black'):
    for field in fieldList:
        nodeList = field.nodeList
        if plotGhost:
            numPoints = nodeList.numNodes
        else:
            numPoints = nodeList.numInternalNodes
        x = []
        y = []
        print numPoints
        for nodeID in xrange(numPoints):
            x.append(nodeList.positions[nodeID].x)
            y.append(field[nodeID])
        xarray = array(x)
        yarray = array(y)
        print x
        print y
        plg(yarray, xarray, color=color)

################################################################################
n1, n2 = 100, 25
rho1, rho2 = 1.0, 0.25
m1, m2 = 0.5*rho1/n1, 0.5*rho2/n2
eps1, eps2 = 1.0, 0.1
v1, v2 = Vector1d(1.0), Vector1d(-1.0)
x0, x1 = 0.0, 1.0

gamma = 2.0
mu = 1.0

neighborSearchType = 3 # GatherScatter
numGridLevels = 10
topGridCellSize = 2.0
origin = Vector1d(0.0)

################################################################################
title('1-D FieldList Test')

eos = GammaLawGasMKS1d(gamma, mu)

nodes1 = SphNodeList1d(n1, eos)
nodes2 = SphNodeList1d(n2, eos)
output('nodes1.numNodes')
output('nodes2.numNodes')

W = BSplineKernel1d()
#W = W4SplineKernel1d()
#W = GaussianKernel1d()
#W = SuperGaussianKernel1d()
#W = PiGaussianKernel1d(1.0)
output('W')
kernelExtent = W.kernelExtent

# Set the table kernel for the FieldList divergence.
WT = TableKernel1d()
WT.setTableData(W, 100)
output('WT')

dx1 = 0.5*(x1 - x0)/n1
dx2 = 0.5*(x1 - x0)/n2
for i in xrange(n1):
    nodes1.positions[i] = (i + 0.5)*dx1
for i in xrange(n2):
    nodes2.positions[i] = 0.5 + (i + 0.5)*dx2

output('nodes1.positions[:]')
output('nodes2.positions[:]')

nodes1.mass[:] = [m1]*nodes1.numNodes
nodes2.mass[:] = [m2]*nodes2.numNodes

nodes1.specificThermalEnergy[:] = [eps1]*nodes1.numNodes
nodes2.specificThermalEnergy[:] = [eps2]*nodes2.numNodes

nodes1.velocity[:] = [v1]*nodes1.numNodes
nodes2.velocity[:] = [v2]*nodes2.numNodes

h1 = 1.0/(2.01*dx1)
h2 = 1.0/(2.01*dx2)
for H in nodes1.Hfield:
    H.xx = h1
for H in nodes2.Hfield:
    H.xx = h2
output('nodes1.Hfield[:]')
output('nodes2.Hfield[:]')

neighbor1 = NestedGridNeighbor1d(nodes1,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes1.neighbor = neighbor1

neighbor2 = NestedGridNeighbor1d(nodes2,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes2.neighbor = neighbor2

nodes1.massDensity[:] = [rho1]*nodes1.numNodes
nodes2.massDensity[:] = [rho2]*nodes2.numNodes
output('nodes1.massDensity[:]')
output('nodes2.massDensity[:]')

##output('nodes.updateMassDensity(W)')
output('nodes1.updateWeight()')
output('nodes2.updateWeight()')

db = DataBase1d()
output('db')
output('db.appendNodeList(nodes1)')
output('db.appendNodeList(nodes2)')
output('db.numNodeLists')
output('db.numFluidNodeLists')

output('db.globalMass[:]')
output('db.fluidMass[:]')
output('db.globalPosition[:]')
output('db.fluidMassDensity[:]')
output('db.fluidSpecificThermalEnergy[:]')
output('db.fluidVelocity[:]')
output('db.fluidWeight[:]')
output('db.fluidHfield[:]')

pressure = db.fluidPressure
output('pressure[:]')
output('pressure[0][:]')
output('pressure[1][:]')

################################################################################
# Test the smoothing function of the FieldList
title('1-D FieldList.smoothFields() test')

fluidPosition = db.fluidPosition
fluidWeight = db.fluidWeight
fluidHfield = db.fluidHfield
smoothPressure = pressure.smoothFields(fluidPosition, fluidWeight, fluidHfield, WT)
monotonicPressure = pressure.smoothFieldsMash(fluidPosition, fluidWeight, fluidHfield, WT)

# Plot the results of the smoothing operation.
from gist import *
window(0)
plotFieldList(pressure, color='black')
plotFieldList(smoothPressure, color='red')
plotFieldList(monotonicPressure, color='green')
pltitle('1-D FieldList smoothFields -- no Boundary conditions')

################################################################################
# Try applying boundary conditions to the pressure, and redo the smoothing.
title('1-D FieldList.smoothFields() test with reflecting boundaries.')

xPlane0, xPlane1 = Plane1d((0.0), (1.0)), Plane1d((1.0), (-1.0))
xbc0, xbc1 = ReflectingBoundary1d(xPlane0), ReflectingBoundary1d(xPlane1)
output('xbc0, xbc1')
output('xbc0.setGhostNodes(nodes1)')
output('xbc0.applyScalarBoundary(pressure[0])')
output('xbc0.applySymTensorBoundary(fluidHfield[0])')
output('xbc0.applyScalarBoundary(fluidWeight[0])')
output('xbc1.setGhostNodes(nodes1)')
output('xbc1.applyScalarBoundary(pressure[0])')
output('xbc1.applySymTensorBoundary(fluidHfield[0])')
output('xbc1.applyScalarBoundary(fluidWeight[0])')
output('xbc0.setGhostNodes(nodes2)')
output('xbc0.applyScalarBoundary(pressure[1])')
output('xbc0.applySymTensorBoundary(fluidHfield[1])')
output('xbc0.applyScalarBoundary(fluidWeight[1])')
output('xbc1.setGhostNodes(nodes2)')
output('xbc1.applyScalarBoundary(pressure[1])')
output('xbc1.applySymTensorBoundary(fluidHfield[1])')
output('xbc1.applyScalarBoundary(fluidWeight[1])')

output('nodes1.neighbor.updateNodes()')
output('nodes2.neighbor.updateNodes()')

smoothPressure = pressure.smoothFields(fluidPosition, fluidWeight, fluidHfield, WT)
monotonicPressure = pressure.smoothFieldsMash(fluidPosition, fluidWeight, fluidHfield, WT)

# Plot the results of the smoothing operation.
from gist import *
window(1)
plotFieldList(pressure, color='black')
plotFieldList(smoothPressure, color='red')
plotFieldList(monotonicPressure, color='green')
pltitle('Reflecting Boundary conditions')
