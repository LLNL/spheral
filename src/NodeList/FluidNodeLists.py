from SpheralModules.Spheral import Vector1d, Vector2d, Vector3d
from SpheralModules.Spheral.NodeSpace import *
from SpheralModules.Spheral.NeighborSpace import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The generic FluidNodeList pattern.
#-------------------------------------------------------------------------------
FluidNodeListFactoryString = """
def makeFluidNodeList%(dim)s(name,
                             eos,
                             numInternal = 0,
                             numGhost = 0,
                             hmin = 1.0e-20,
                             hmax = 1.0e20,
                             hminratio = 0.1,
                             nPerh = 2.01,
                             maxNumNeighbors = 500,
                             rhoMin = 1.0e-10,
                             rhoMax = 1e10,
                             #NeighborType = TreeNeighbor%(dim)s,
                             searchType = GatherScatter,
                             #numGridLevels = 31,
                             #topGridCellSize = 100.0,
                             #origin = Vector%(dim)s.zero,
                             kernelExtent = 2.0,
                             #gridCellInfluenceRadius = 1,
                             xmin = Vector%(dim)s.one * -10.0,
                             xmax = Vector%(dim)s.one *  10.0):
    result = FluidNodeList%(dim)s(name, eos, numInternal, numGhost, 
                                  hmin, hmax, hminratio, 
                                  nPerh, maxNumNeighbors,
                                  rhoMin, rhoMax)
        # if NeighborType == NestedGridNeighbor%(dim)s:
        #     print "makeFluidNodeList Deprecation Warning: NestedGridNeighbor is deprecated and will no longer be used.  Using TreeNeighbor instead."
        # result._neighbor = NestedGridNeighbor%(dim)s(result, searchType, 
        #                                              numGridLevels, topGridCellSize, 
        #                                              origin, kernelExtent, 
        #                                              gridCellInfluenceRadius)
        # else:
    result._neighbor = TreeNeighbor%(dim)s(result, searchType, kernelExtent, xmin, xmax)
    result.registerNeighbor(result._neighbor)
    result.eos = eos
    return result
"""

#-------------------------------------------------------------------------------
# Create the different dimension implementations.
#-------------------------------------------------------------------------------
for dim in dims:
    exec(FluidNodeListFactoryString % {"dim"    : "%id" % dim})
