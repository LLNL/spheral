from SpheralCompiledPackages import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The generic VoidNodeList pattern.
#-------------------------------------------------------------------------------
def VoidNodeListFactory(ndim):
    suffix = "{}d".format(ndim)
    TreeNeighbor = eval("TreeNeighbor" + suffix)
    NestedGridNeighbor = eval("NestedGridNeighbor" + suffix)
    Vector = eval("Vector" + suffix)
    NodeList = eval("NodeList" + suffix)
    def factory(name,
                numInternal = 0,
                numGhost = 0,
                hmin = 1.0e-20,
                hmax = 1.0e20,
                hminratio = 0.1,
                nPerh = 2.01,
                maxNumNeighbors = 500,

                # Neighboring stuff
                NeighborType = TreeNeighbor,
                searchType = GatherScatter,
                kernelExtent = 2.0,

                # Parameters for TreeNeighbor
                xmin = Vector.one * -10.0,
                xmax = Vector.one *  10.0):
        result = NodeList(name, numInternal, numGhost, 
                          hmin, hmax, hminratio, 
                          nPerh, maxNumNeighbors)
        if NeighborType == TreeNeighbor:
            result._neighbor = TreeNeighbor(result, searchType, kernelExtent, xmin, xmax)
        elif NeighborType == NestedGridNeighbor:
            print("VoidNodeList Deprecation Warning: NestedGridNeighbor is deprecated: suggest using TreeNeighbor.")
            result._neighbor = NestedGridNeighbor(result, searchType, 
                                                  kernelExtent = kernelExtent)
        else:
            raise ValueError("Unknown NeighborType")
        result.registerNeighbor(result._neighbor)
        return result
    return factory

#-------------------------------------------------------------------------------
# Create the different dimension implementations.
#-------------------------------------------------------------------------------
for ndim in dims:
    exec("makeVoidNodeList{ndim}d = VoidNodeListFactory({ndim})".format(ndim=ndim))
