from SpheralCompiledPackages import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The generic SolidNodeList defintion.
#-------------------------------------------------------------------------------
def SolidNodeListFactory(ndim):
    suffix = "{}d".format(ndim)
    SolidNodeList = eval("SolidNodeList" + suffix)
    TreeNeighbor = eval("TreeNeighbor" + suffix)
    NestedGridNeighbor = eval("NestedGridNeighbor" + suffix)
    Vector = eval("Vector" + suffix)
    NullStrength = eval("NullStrength" + suffix)
    def factory(name,
                eos,
                strength = NullStrength(),
                numInternal = 0,
                numGhost = 0,
                hmin = 1.0e-20,
                hmax = 1.0e20,
                hminratio = 0.1,
                nPerh = 2.01,
                maxNumNeighbors = 500,
                rhoMin = 1.0e-10,
                rhoMax = 1e10,

                # Neighboring stuff
                NeighborType = TreeNeighbor,
                searchType = GatherScatter,
                kernelExtent = 2.0,

                # Parameters for TreeNeighbor
                xmin = Vector.one * -10.0,
                xmax = Vector.one *  10.0):
        result = SolidNodeList(name = name,
                               eos = eos,
                               strength = strength,
                               numInternal = numInternal,
                               numGhost = numGhost,
                               hmin = hmin,
                               hmax = hmax,
                               hminratio = hminratio, 
                               nPerh = nPerh,
                               maxNumNeighbors = maxNumNeighbors,
                               rhoMin = rhoMin,
                               rhoMax = rhoMax)
        result.eos = eos
        result.strength = strength
        if NeighborType == TreeNeighbor:
            result._neighbor = TreeNeighbor(result, searchType, kernelExtent, xmin, xmax)
        elif NeighborType == NestedGridNeighbor:
            print("SolidNodeList Deprecation Warning: NestedGridNeighbor is deprecated: suggest using TreeNeighbor.")
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
    exec("makeSolidNodeList{ndim}d = SolidNodeListFactory({ndim})".format(ndim=ndim))
