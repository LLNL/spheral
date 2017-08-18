from SpheralModules.Spheral import Vector1d, Vector2d, Vector3d
from SpheralModules.Spheral.NodeSpace import *
from SpheralModules.Spheral.NeighborSpace import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The generic VoidNodeList pattern.
#-------------------------------------------------------------------------------
VoidNodeListFactoryString = """
def makeVoidNodeList%(dim)s(name,
                            numInternal = 0,
                            numGhost = 0,
                            hmin = 1.0e-20,
                            hmax = 1.0e20,
                            hminratio = 0.1,
                            nPerh = 2.01,
                            maxNumNeighbors = 500,
                            searchType = GatherScatter,
                            numGridLevels = 31,
                            topGridCellSize = 100.0,
                            origin = Vector%(dim)s.zero,
                            kernelExtent = 2.0,
                            gridCellInfluenceRadius = 1):
    result = NodeList%(dim)s(name, numInternal, numGhost, 
                             hmin, hmax, hminratio, 
                             nPerh, maxNumNeighbors)
    result._neighbor = NestedGridNeighbor%(dim)s(result, searchType, 
                                                 numGridLevels, topGridCellSize, 
                                                 origin, kernelExtent, 
                                                 gridCellInfluenceRadius)
    result.registerNeighbor(result._neighbor)
    return result
"""

#-------------------------------------------------------------------------------
# Create the different dimension implementations.
#-------------------------------------------------------------------------------
for dim in dims:
    exec(VoidNodeListFactoryString % {"dim"    : "%id" % dim})
