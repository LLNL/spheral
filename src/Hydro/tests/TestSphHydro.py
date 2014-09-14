from Spheral import *
from SpheralTestUtilities import *

def plotFieldList(fieldList, plotGhost=0):
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
            y.append(field[nodeID].x)
        xarray = array(x)
        yarray = array(y)
        print x
        print y
        plg(yarray, xarray)

################################################################################
nx, ny = 50, 50
n1 = nx*ny
rho1 = 1.0
m1 = rho1/n1
eps0 = 0.0
epsMultiplier = 1.0
vr = -1.0
x0, x1 = 0.0, 1.0
y0, y1 = 0.0, 1.0

gamma = 5.0/3.0
mu = 1.0
Cl = 1.0
Cq = 1.0

neighborSearchType = 3 # GatherScatter
numGridLevels = 10
topGridCellSize = 0.5
origin = Vector2d(0.0)

copyFields = 1

################################################################################
title('2-D Hydro derivative Test')

eos = GammaLawGasMKS2d(gamma, mu)

nodes1 = SphVarGradNodeList2d(n1, eos)
output('nodes1.numNodes')

W = BSplineKernel2d()
#W = W4SplineKernel2d()
#W = GaussianKernel2d()
#W = SuperGaussianKernel2d()
#W = PiGaussianKernel2d(1.0)
output('W')
kernelExtent = W.kernelExtent

# Set the table kernel.
WT = TableKernel2d()
WT.setTableData(W, 100)
output('WT')

# Set node positions
dx, dy = (x1 - x0)/nx, (y1 - y0)/ny
for iy in xrange(ny):
    for ix in xrange(nx):
        nodeID = ix + nx*iy
        nodes1.positions[nodeID] = ((ix + 0.5)*dx, (iy + 0.5)*dy)

# Set node masses
nodes1.mass[:] = [m1]*nodes1.numNodes

# Set node specific thermal energies
for nodeID in xrange(nodes1.numNodes):
    nodes1.specificThermalEnergy[nodeID] = eps0 + epsMultiplier*(nodes1.positions[nodeID].magnitude2())

# Set node velocities
for nodeID in xrange(nodes1.numNodes):
    unit = nodes1.positions[nodeID]
    mag = unit.magnitude() + 1.0e-10
    unit = Vector2d(unit.x/mag, unit.y/mag)
    nodes1.velocity[nodeID] = (vr*unit.x, vr*unit.y)

# Set the smoothing scales.
h = 1.0/(2.01*dx)
for H in nodes1.Hfield:
    H.xx = h
    H.yy = h

# Set the mass densities.
nodes1.massDensity[:] = [rho1]*nodes1.numNodes

# Update the nodal weights
updateTimer = SpheralTimer('Update node weights')
updateTimer.start()
nodes1.updateWeight()
updateTimer.stop()
updateTimer.printStatus()

# Construct neighbor objects and associate them with the node lists.
neighborTimer = SpheralTimer('Neighbor initialization.')
neighborTimer.start()
neighbor1 = NestedGridNeighbor2d(nodes1,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes1.neighbor = neighbor1
neighborTimer.stop()
neighborTimer.printStatus()

# Construct a DataBase to hold our two node lists
db = DataBase2d()
output('db')
output('db.appendNodeList(nodes1)')
output('db.numNodeLists')
output('db.numFluidNodeLists')

# Plot the velocity field.
from SpheralGistUtilities import *
from gist import *
xnodes1, ynodes1 = nodePositions2Mesh2d(nodes1, nx, ny)
vxnodes1, vynodes1 = vectorField2NumericArrays2d(nodes1.velocity, nx, ny)
window(0)
plv(vynodes1, vxnodes1, ynodes1, xnodes1, color='blue')
pltitle('Velocity Field.')

################################################################################
# Apply reflecting boundary conditions along the x=0 and y=0 lines.
title('Apply boundary conditions.')
xPlane0 = Plane2d((0.0, 0.0), (1.0, 0.0))
yPlane0 = Plane2d((0.0, 0.0), (0.0, 1.0))
xbc0 = ReflectingBoundary2d(xPlane0)
ybc0 = ReflectingBoundary2d(yPlane0)
output('xbc0, ybc0')

bcTimer = SpheralTimer('Apply boundary conditions')
bcTimer.start()
for bc in [xbc0, ybc0]:
    output('bc.setGhostNodes(nodes1)')
    output('bc.applyScalarBoundary(nodes1.mass)')
    output('bc.applyScalarBoundary(nodes1.massDensity)')
    output('bc.applyVectorBoundary(nodes1.velocity)')
    output('bc.applyScalarBoundary(nodes1.specificThermalEnergy)')
    output('bc.applyScalarBoundary(nodes1.weight)')
    output('bc.applySymTensorBoundary(nodes1.Hfield)')
bcTimer.stop()
bcTimer.printStatus()

neighborTimer.start()
output('nodes1.neighbor.updateNodes()')
neighborTimer.stop()
neighborTimer.printStatus()

################################################################################
# Test the hydro packages derivative
title('2-D Hydro test')

q = MonaghanGingoldViscosity2d(Cl, Cq)
output('q')

hydro = Hydro2d(WT, WT, q)
output('hydro')
output('hydro.kernel')
output('hydro.PiKernel')
output('hydro.valid')

constructFieldsTimer = SpheralTimer('Construct FieldLists')
constructFieldsTimer.start()

print 'Build FieldLists'
rhoSum = ScalarFieldList2d()
maxQPressure = ScalarFieldList2d()
DrhoDt = ScalarFieldList2d()
DvDt = VectorFieldList2d()
DepsDt = ScalarFieldList2d()
DvDx = TensorFieldList2d()
DHDt = SymTensorFieldList2d()
rhoSum.copyFields()
maxQPressure.copyFields()
DrhoDt.copyFields()
DvDt.copyFields()
DepsDt.copyFields()
DvDx.copyFields()
DHDt.copyFields()
##rhoSum = ScalarFieldList2d(copyFields)
##maxQPressure = ScalarFieldList2d(copyFields)
##DrhoDt = ScalarFieldList2d(copyFields)
##DvDt = VectorFieldList2d(copyFields)
##DepsDt = ScalarFieldList2d(copyFields)
##DvDx = TensorFieldList2d(copyFields)
##DHDt = SymTensorFieldList2d(copyFields)
print 'done.'

print 'Add Fields to FieldLists'
rhoSum.appendField(ScalarField2d(nodes1))
maxQPressure.appendField(ScalarField2d(nodes1))
DrhoDt.appendField(ScalarField2d(nodes1))
DvDt.appendField(VectorField2d(nodes1))
DepsDt.appendField(ScalarField2d(nodes1))
DvDx.appendField(TensorField2d(nodes1))
DHDt.appendField(SymTensorField2d(nodes1))
constructFieldsTimer.stop()
constructFieldsTimer.printStatus()
print 'done.'

derivTimer = SpheralTimer('Calculate derivatives')
derivTimer.start()
hydro.evaluateDerivatives(db, 0.0, 0.0, rhoSum, maxQPressure, DrhoDt, DvDt, DepsDt, DvDx, DHDt)
derivTimer.stop()
derivTimer.printStatus()

# Plot Drho/Dt as a function of radius
window(1)
r1, DrhoDtr = radialProfile(DrhoDt[0])
plg(DrhoDtr[:nodes1.numInternalNodes], r1[:nodes1.numInternalNodes],
    type=0, marker='\4', color='blue')
pltitle('Drho/Dt as a function of radius.')

# Plot the analytic solution for Drho/Dt, which in this case is rho0/r.
answer1 = rho1/r1
plg(answer1[:nodes1.numInternalNodes], r1[:nodes1.numInternalNodes], type=0, marker='\2', color='black')

# Plot DH/Dt as a function of radius
window(2)
DHDtr = array(map(lambda x: x.xx, DHDt[0]))
plg(DHDtr[:nodes1.numInternalNodes], r1[:nodes1.numInternalNodes],
    type=0, marker='\4', color='blue')
pltitle('DH/Dt as a function of radius.')

# Plot the analytic solution for DH/Dt, which in this case is H0/2r.
answer2 = 0.5*h/r1
plg(answer2[:nodes1.numInternalNodes], r1[:nodes1.numInternalNodes], type=0, marker='\2', color='black')

# Plot ||DvDx|| as a function of radius.
window(3)
DetDvDxr = array([0.0]*nodes1.numNodes)
for i in xrange(nodes1.numNodes):
    DetDvDxr[i] = DvDx[0][i].Determinant()
plg(DetDvDxr[:nodes1.numInternalNodes], r1[:nodes1.numInternalNodes], type=0, marker='\2', color='blue')
pltitle('||DvDx|| as a function of radius.')

# Plot |DvDt| as a function of radius.
window(4)
DvDtMag = array(map(lambda x: x.magnitude(), DvDt[0]))
for i in xrange(nodes1.numNodes):
    DvDtMag[i] = DvDt[0][i].magnitude()
plg(DvDtMag[:nodes1.numInternalNodes], r1[:nodes1.numInternalNodes], type=0, marker='\2', color='blue')
pltitle('|DvDt| as a function of radius.')
