from math import *
from Spheral import *

import random
g = random.Random()

#===============================================================================
# Calculate one over the smoothing scale for the given number of nodes and
# volume.
#===============================================================================
def determineH(nGlobal, xRange,
               nNodesPerh = 2.01):
    assert nGlobal > 0
    
    vol = (xRange[1][0] - xRange[0][0])*(xRange[1][1] - xRange[0][1])
    assert vol > 0.0
    dV = vol/nGlobal
    dx = sqrt(dV)
    hi = 1.0/(nNodesPerh*dx)
    Hi = SymTensor2d(hi, 0.0,
                     0.0, hi)
    return Hi

#===============================================================================
# Distribute nodes randomly in the given volume.
#===============================================================================
def randomDistribute(nNodesGlobal,   # global number of nodes in this nodelist
                     xyRangeTotal):  # total simulation volume

    nodePositions = []
    for globalNodeID in xrange(nNodesGlobal):
        nodePositions.append(Vector2d(g.uniform(xyRangeTotal[0][0],
                                                xyRangeTotal[1][0]),
                                      g.uniform(xyRangeTotal[0][1],
                                                xyRangeTotal[1][1])))

    return nodePositions

#===============================================================================
# Create a set of NodeLists and a DataBase for use in the tests.
#===============================================================================
# Generic parameters for 2-D tests.
n1 = 1000
n2 = 2500
n3 = 500

range1 = [(0.0, 0.0), (1.0, 1.0)]
range2 = [(1.0, 0.0), (1.5, 1.0)]
range3 = [(1.5, 0.0), (2.0, 1.0)]

neighborSearchType = Neighbor2d.NeighborSearchType.GatherScatter
numGridLevels = 10
topGridCellSize = 1.0
origin = Vector2d(0.0, 0.0)
kernelExtent = 2.0

# Construct the NodeLists to be distributed
nodes1 = SphNodeList2d()
nodes2 = SphNodeList2d()
nodes3 = SphNodeList2d()
for thpt in ((nodes1, n1, range1),
             (nodes2, n2, range2),
             (nodes3, n3, range3)):
    nodes = thpt[0]
    nGlobal = thpt[1]
    globalRange = thpt[2]
    xyNodes = randomDistribute(nGlobal,
                               globalRange)
    n = nGlobal
    nodes.numInternalNodes = n
    Hi = determineH(nGlobal, globalRange)
    for i in xrange(n):
        nodes.mass()[i] = 1.0
        nodes.positions()[i] = xyNodes[i]
        nodes.Hfield()[i] = Hi

neighbor1 = NestedGridNeighbor2d(nodes1,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
neighbor2 = NestedGridNeighbor2d(nodes2,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
neighbor3 = NestedGridNeighbor2d(nodes3,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes1.registerNeighbor(neighbor1)
nodes2.registerNeighbor(neighbor2)
nodes3.registerNeighbor(neighbor3)

# Put the NodeLists into a DataBase.
dataBase = DataBase2d()
dataBase.appendNodeList(nodes1)
dataBase.appendNodeList(nodes2)
dataBase.appendNodeList(nodes3)

