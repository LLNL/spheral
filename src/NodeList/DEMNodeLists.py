from SpheralCompiledPackages import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The generic SolidNodeList defintion.
#-------------------------------------------------------------------------------
DEMNodeListFactoryString = """

def makeDEMNodeList%(dim)s(name,
                           neighborSearchBuffer = 0.1,
                           numInternal = 0,
                           numGhost = 0,
                           hmin = 1.0e-20,
                           hmax = 1.0e20,
                           hminratio = 0.1,
                           nPerh = 1.01,
                           maxNumNeighbors = 500,

                           # Neighboring stuff
                           NeighborType = TreeNeighbor%(dim)s,
                           searchType = GatherScatter,
                           kernelExtent = 1.0,

                           # Parameters only for NestedGridNeighbor (deprecated)
                           # numGridLevels = 31,
                           # topGridCellSize = 100.0,
                           # origin = Vector%(dim)s.zero,
                           # gridCellInfluenceRadius = 1,

                           # Parameters for TreeNeighbor
                           xmin = Vector%(dim)s.one * -10.0,
                           xmax = Vector%(dim)s.one *  10.0):
    result = DEMNodeList%(dim)s(name, numInternal, numGhost, 
                                  hmin, hmax, hminratio, 
                                  nPerh, neighborSearchBuffer, maxNumNeighbors)

    if NeighborType == NestedGridNeighbor%(dim)s:
        print "makeDEMNodeList Deprecation Warning: NestedGridNeighbor is deprecated: suggest using TreeNeighbor."
        result._neighbor = NestedGridNeighbor%(dim)s(result, searchType, 
                                                     kernelExtent = kernelExtent)
                                                     #numGridLevels, topGridCellSize, 
                                                     #origin, kernelExtent, 
                                                     #gridCellInfluenceRadius)
    else:
        result._neighbor = TreeNeighbor%(dim)s(result, searchType, kernelExtent, xmin, xmax)
    result.registerNeighbor(result._neighbor)
    return result
"""

#-------------------------------------------------------------------------------
# Create the different SolidNodeLists.
#-------------------------------------------------------------------------------
for dim in dims:
    exec(DEMNodeListFactoryString % {"dim" : "%id" % dim})
