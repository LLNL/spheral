from Spheral import *
from SpheralTestUtilities import *
from NeighborTestUtilities import *

################################################################################
title('3-D Sph node list')
nxNodes = 10
nyNodes = 10
nzNodes = 10
nNodes = nxNodes*nyNodes*nzNodes
h = 1.0/1.01

gamma = 5.0/3.0
mu = 1.0

neighborSearchType = 3 # GatherScatter
numGridLevels = 10
topGridCellSize = 2.5
origin = Vector3d(0.0, 0.0, 0.0)
kernelExtent = 2.0

eos = GammaLawGasCGS3d(gamma, mu)
nodes = SphNodeList3d(eos, nNodes)

dx = 1.0
dy = 1.0
dz = 1.0
nxyNodes = nxNodes*nyNodes
for i in xrange(nzNodes):
    for j in xrange(nyNodes):
        for k in xrange(nxNodes):
            iNode = i*nxyNodes + j*nxNodes + k
            nodes.positions[iNode] = Vector3d((k + 0.5)*dx,
                                              (j + 0.5)*dy,
                                              (i + 0.5)*dz)
            nodes.Hfield[iNode].xx = h
            nodes.Hfield[iNode].yy = h
            nodes.Hfield[iNode].zz = h

nestedNeighbor = NestedGridNeighbor3d(nodes,
                                      neighborSearchType,
                                      numGridLevels,
                                      topGridCellSize,
                                      origin,
                                      kernelExtent)

# Create reflecting boundary conditions.
xPlane0, xPlane1 = Plane3d((0,0,0), (1,0,0)), Plane3d((10,0,0), (-1,0,0))
yPlane0, yPlane1 = Plane3d((0,0,0), (0,1,0)), Plane3d((0,10,0), (0,-1,0))
zPlane0, zPlane1 = Plane3d((0,0,0), (0,0,1)), Plane3d((0,0,10), (0,0,-1))
xbc0, xbc1 = ReflectingBoundary3d(xPlane0), ReflectingBoundary3d(xPlane1)
ybc0, ybc1 = ReflectingBoundary3d(yPlane0), ReflectingBoundary3d(yPlane1)
zbc0, zbc1 = ReflectingBoundary3d(zPlane0), ReflectingBoundary3d(zPlane1)

output('xbc0, xbc1')
output('ybc0, ybc1')
output('zbc0, zbc1')

# Apply the boundary conditions to the node list.
print 'Applying boundary conditions to node list.'
output('nodes.numNodes, nodes.numInternalNodes, nodes.numGhostNodes')
output('xbc0.setNodeListGhostNodes(nodes)')
output('nodes.numNodes, nodes.numInternalNodes, nodes.numGhostNodes')
output('xbc1.setNodeListGhostNodes(nodes)')
output('nodes.numNodes, nodes.numInternalNodes, nodes.numGhostNodes')
output('ybc0.setNodeListGhostNodes(nodes)')
output('nodes.numNodes, nodes.numInternalNodes, nodes.numGhostNodes')
output('ybc1.setNodeListGhostNodes(nodes)')
output('nodes.numNodes, nodes.numInternalNodes, nodes.numGhostNodes')
output('zbc0.setNodeListGhostNodes(nodes)')
output('nodes.numNodes, nodes.numInternalNodes, nodes.numGhostNodes')
output('zbc1.setNodeListGhostNodes(nodes)')
output('nodes.numNodes, nodes.numInternalNodes, nodes.numGhostNodes')

neighborTimer = SpheralTimer("3-D nested neighbor timer")
neighborTimer.start()
output('nestedNeighbor.setMasterList(0)')
neighborTimer.stop()
output('nestedNeighbor.masterList[:]')
output('nestedNeighbor.coarseNeighborList[:]')
output('nestedNeighbor.setRefineNeighborList(0)')
output('nestedNeighbor.refineNeighborList[:]')

# Find the correct answer, and make sure all the actual neighbor nodes are
# included in the refined list.
checkTimer = SpheralTimer("3-D check neighbor timer")
checkTimer.start()
test = checkNeighbors(nestedNeighbor.refineNeighborList,
                      findNeighborNodes(nodes.positions[0],
                                        kernelExtent/nodes.Hfield[0].xx,
                                        nodes))
checkTimer.stop()

if test:
    print '3-D Sph Nested Neighbor test PASSED'
else:
    print '3-D Sph Nested Neighbor test FAILED'


neighborTimer.printStatus()
checkTimer.printStatus()

################################################################################
title('2-D Sph node list')
nxNodes = 10
nyNodes = 10
nNodes = nxNodes*nyNodes
h = 1.0/1.01

gamma = 5.0/3.0
mu = 1.0

neighborSearchType = 3 # GatherScatter
numGridLevels = 10
topGridCellSize = 2.5
origin = Vector2d(0.0, 0.0)
kernelExtent = 2.0

eos = GammaLawGasCGS2d(gamma, mu)
nodes = SphNodeList2d(eos, nNodes)

dx = 1.0
dy = 1.0
for i in xrange(nyNodes):
    for j in xrange(nxNodes):
        iNode = i*nxNodes + j
        nodes.positions[iNode] = Vector2d((j + 0.5)*dx,
                                          (i + 0.5)*dy)
        nodes.Hfield[iNode].xx = h
        nodes.Hfield[iNode].yy = h

nestedNeighbor = NestedGridNeighbor2d(nodes,
                                      neighborSearchType,
                                      numGridLevels,
                                      topGridCellSize,
                                      origin,
                                      kernelExtent)

# Create reflecting boundary conditions.
xPlane0, xPlane1 = Plane2d((0,0), (1,0)), Plane2d((10,0), (-1,0))
yPlane0, yPlane1 = Plane2d((0,0), (0,1)), Plane2d((0,10), (0,-1))
xbc0, xbc1 = ReflectingBoundary2d(xPlane0), ReflectingBoundary2d(xPlane1)
ybc0, ybc1 = ReflectingBoundary2d(yPlane0), ReflectingBoundary2d(yPlane1)

output('xbc0, xbc1')
output('ybc0, ybc1')

# Apply the boundary conditions to the node list.
print 'Applying boundary conditions to node list.'
output('nodes.numNodes, nodes.numInternalNodes, nodes.numGhostNodes')
output('xbc0.setNodeListGhostNodes(nodes)')
output('nodes.numNodes, nodes.numInternalNodes, nodes.numGhostNodes')
plotNodes(nodes, 0, 'After applying xbc0')
output('xbc1.setNodeListGhostNodes(nodes)')
output('nodes.numNodes, nodes.numInternalNodes, nodes.numGhostNodes')
plotNodes(nodes, 1, 'After applying xbc1')
output('ybc0.setNodeListGhostNodes(nodes)')
output('nodes.numNodes, nodes.numInternalNodes, nodes.numGhostNodes')
plotNodes(nodes, 2, 'After applying ybc0')
output('ybc1.setNodeListGhostNodes(nodes)')
output('nodes.numNodes, nodes.numInternalNodes, nodes.numGhostNodes')
plotNodes(nodes, 3, 'After applying ybc1')

neighborTimer = SpheralTimer("2-D nested neighbor timer")
neighborTimer.start()
output('nestedNeighbor.setMasterList(0)')
neighborTimer.stop()
output('nestedNeighbor.masterList[:]')
output('nestedNeighbor.coarseNeighborList[:]')
output('nestedNeighbor.setRefineNeighborList(0)')
output('nestedNeighbor.refineNeighborList[:]')

# Find the correct answer, and make sure all the actual neighbor nodes are
# included in the refined list.
checkTimer = SpheralTimer("2-D check neighbor timer")
checkTimer.start()
test = checkNeighbors(nestedNeighbor.refineNeighborList,
                      findNeighborNodes(nodes.positions[0],
                                        kernelExtent/nodes.Hfield[0].xx,
                                        nodes))
checkTimer.stop()

if test:
    print '2-D Sph Nested Neighbor test PASSED'
else:
    print '2-D Sph Nested Neighbor test FAILED'

neighborTimer.printStatus()
checkTimer.printStatus()

