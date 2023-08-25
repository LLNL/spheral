from SpheralCompiledPackages import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The generic SolidNodeList defintion.
#-------------------------------------------------------------------------------
def DEMNodeListFactory(ndim):
    suffix = "{}d".format(ndim)
    DEMNodeList = eval("DEMNodeList" + suffix)
    TreeNeighbor = eval("TreeNeighbor" + suffix)
    NestedGridNeighbor = eval("NestedGridNeighbor" + suffix)
    Vector = eval("Vector" + suffix)
    def factory(name,
                numInternal = 0,
                numGhost = 0,
                hmin = 1.0e-20,
                hmax = 1.0e20,
                hminratio = 0.1,
                nPerh = 1.01,
                maxNumNeighbors = 500,
                neighborSearchBuffer=0.1,

                # Neighboring stuff
                NeighborType = TreeNeighbor,
                searchType = GatherScatter,
                kernelExtent = 1.0,

                # Parameters for TreeNeighbor
                xmin = Vector.one * -10.0,
                xmax = Vector.one *  10.0):
        result = DEMNodeList(name, numInternal, numGhost, 
                             hmin, hmax, hminratio, 
                             nPerh, neighborSearchBuffer, maxNumNeighbors)

        if NeighborType == TreeNeighbor:
            result._neighbor = TreeNeighbor(result, searchType, kernelExtent, xmin, xmax)
        elif NeighborType == NestedGridNeighbor:
            print("DEMNodeList Deprecation Warning: NestedGridNeighbor is deprecated: suggest using TreeNeighbor.")
            result._neighbor = NestedGridNeighbor(result, searchType, 
                                                  kernelExtent = kernelExtent)
        else:
            raise ValueError("Unknown NeighborType")
        result.registerNeighbor(result._neighbor)
        return result
    return factory

#-------------------------------------------------------------------------------
# Create the different SolidNodeLists.
#-------------------------------------------------------------------------------
for ndim in dims:
    exec("makeDEMNodeList{ndim}d = DEMNodeListFactory({ndim})".format(ndim=ndim))
