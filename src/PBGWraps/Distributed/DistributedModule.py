from pybindgen import *

from PBGutils import *
from BoundaryModule import generateBoundaryVirtualBindings
from CXXTypesModule import generateStdVectorBindings, generateStdPairBindings

#-------------------------------------------------------------------------------
# The class for wrapping this module.
#-------------------------------------------------------------------------------
class Distributed:

    #---------------------------------------------------------------------------
    # Add the types to the module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir, dims):

        self.dims = dims

        # Includes.
        mod.add_include('"%s/DistributedTypes.hh"' % srcdir)

        # Namespaces.
        self.space = mod.add_cpp_namespace("Spheral")

        # Expose types.
        for dim in self.dims:
            exec('''
Boundary%(dim)id = findObject(self.space, "Boundary%(dim)id")

self.DistributedBoundary%(dim)id = addObject(self.space, "DistributedBoundary%(dim)id", parent=Boundary%(dim)id, allow_subclassing=True)

self.NestedGridDistributedBoundary%(dim)id = addObject(self.space, "NestedGridDistributedBoundary%(dim)id", parent=self.DistributedBoundary%(dim)id, allow_subclassing=True, is_singleton=True)

self.BoundingVolumeDistributedBoundary%(dim)id = addObject(self.space, "BoundingVolumeDistributedBoundary%(dim)id", parent=self.DistributedBoundary%(dim)id, allow_subclassing=True, is_singleton=True)

self.TreeDistributedBoundary%(dim)id = addObject(self.space, "TreeDistributedBoundary%(dim)id", parent=self.DistributedBoundary%(dim)id, allow_subclassing=True, is_singleton=True)

self.DomainBoundaryNodes%(dim)id = addObject(self.space, "DomainBoundaryNodes", outer_class=self.DistributedBoundary%(dim)id, allow_subclassing=True)

self.DomainNode%(dim)id = addObject(self.space, "DomainNode%(dim)id", allow_subclassing=True)

self.RedistributeNodes%(dim)id = addObject(self.space, "RedistributeNodes%(dim)id", allow_subclassing=True)

#self.NestedGridRedistributeNodes%(dim)id = addObject(self.space, "NestedGridRedistributeNodes%(dim)id", parent=self.RedistributeNodes%(dim)id, allow_subclassing=True)

self.SpaceFillingCurveRedistributeNodes%(dim)id = addObject(self.space, "SpaceFillingCurveRedistributeNodes%(dim)id", parent=self.RedistributeNodes%(dim)id, allow_subclassing=True)

self.MortonOrderRedistributeNodes%(dim)id = addObject(self.space, "MortonOrderRedistributeNodes%(dim)id", parent=self.SpaceFillingCurveRedistributeNodes%(dim)id, allow_subclassing=True)

self.PeanoHilbertOrderRedistributeNodes%(dim)id = addObject(self.space, "PeanoHilbertOrderRedistributeNodes%(dim)id", parent=self.SpaceFillingCurveRedistributeNodes%(dim)id, allow_subclassing=True)

self.SortAndDivideRedistributeNodes%(dim)id = addObject(self.space, "SortAndDivideRedistributeNodes%(dim)id", parent=self.RedistributeNodes%(dim)id, allow_subclassing=True)

self.VoronoiRedistributeNodes%(dim)id = addObject(self.space, "VoronoiRedistributeNodes%(dim)id", parent=self.RedistributeNodes%(dim)id, allow_subclassing=True)

self.vector_of_DomainNode%(dim)id = addObject(mod, "vector_of_DomainNode%(dim)id", allow_subclassing=True)

self.pair_ULL_DomainNode%(dim)id = addObject(mod, "pair_ULL_DomainNode%(dim)id", allow_subclassing=True)

self.vector_of_pair_ULL_DomainNode%(dim)id = addObject(mod, "vector_of_pair_ULL_DomainNode%(dim)id", allow_subclassing=True)
''' % {"dim" : dim})

        if 1 in self.dims:
            self.DistributeByXPosition1d = addObject(self.space, "DistributeByXPosition1d", allow_subclassing=True)
        if 2 in self.dims:
            self.DistributeByXPosition2d = addObject(self.space, "DistributeByXPosition2d", allow_subclassing=True)

        return

    #---------------------------------------------------------------------------
    # Generate the bindings.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        for dim in self.dims:
            exec('''
self.distributedBoundaryBindings(self.DistributedBoundary%(dim)id, %(dim)i)
self.nestedGridDistributedBoundaryBindings(self.NestedGridDistributedBoundary%(dim)id, %(dim)i)
self.boundingVolumeDistributedBoundaryBindings(self.BoundingVolumeDistributedBoundary%(dim)id, %(dim)i)
self.treeDistributedBoundaryBindings(self.TreeDistributedBoundary%(dim)id, %(dim)i)
self.domainBoundaryNodesBindings(self.DomainBoundaryNodes%(dim)id, %(dim)i)
self.domainNodeBindings(self.DomainNode%(dim)id, %(dim)i)
self.redistributeNodesBindings(self.RedistributeNodes%(dim)id, %(dim)i)
#self.nestedGridRedistributeNodesBindings(self.NestedGridRedistributeNodes%(dim)id, %(dim)i)
self.spaceFillingCurveRedistributeNodesBindings(self.SpaceFillingCurveRedistributeNodes%(dim)id, %(dim)i)
self.genericOrderRedistributeNodesBindings(self.MortonOrderRedistributeNodes%(dim)id, %(dim)i)
self.genericOrderRedistributeNodesBindings(self.PeanoHilbertOrderRedistributeNodes%(dim)id, %(dim)i)
self.genericSortAndDivideRedistributeNodesBindings(self.SortAndDivideRedistributeNodes%(dim)id, %(dim)i)
self.voronoiRedistributeNodesBindings(self.VoronoiRedistributeNodes%(dim)id, %(dim)i)

generateStdPairBindings(self.pair_ULL_DomainNode%(dim)id, "uint64_t", "Spheral::DomainNode%(dim)id", "pair_ULL_DomainNode%(dim)id")

generateStdVectorBindings(self.vector_of_DomainNode%(dim)id, "Spheral::DomainNode%(dim)id", "vector_of_DomainNode%(dim)id", True)

generateStdVectorBindings(self.vector_of_pair_ULL_DomainNode%(dim)id, "pair_ULL_DomainNode%(dim)id", "vector_of_pair_ULL_DomainNode%(dim)id", True)
''' % {"dim" : dim})

        if 1 in self.dims:
            self.distributeByXPositionBindings(self.DistributeByXPosition1d, 1)
        if 2 in self.dims:
            self.distributeByXPositionBindings(self.DistributeByXPosition2d, 2)

        return

    #---------------------------------------------------------------------------
    # DistributedBoundary
    #---------------------------------------------------------------------------
    def distributedBoundaryBindings(self, x, ndim):

        me = "DistributedBoundary%id" % ndim
        dim = "Spheral::Dim<%i> " % ndim
        vector = "Vector%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        domainboundarynodes = "Spheral::DistributedBoundary%id::DomainBoundaryNodes" % ndim
        intfield = "Spheral::IntField%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        tensorfield = "Spheral::TensorField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        thirdranktensorfield = "Spheral::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::VectorDoubleField%id" % ndim
        vectorvectorfield = "Spheral::VectorVectorField%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        nestedgridneighbor = "Spheral::NestedGridNeighbor%id" % ndim

        # Constructors.
        x.add_constructor([])

        # Methods.
        x.add_method("domainBoundaryNodes", domainboundarynodes, [constrefparam(nodelist, "nodeList"), param("int", "neighborDomainID")], is_const=True)
        x.add_method("communicatedProcs", None, [refparam("vector_of_int", "sendProcs"), refparam("vector_of_int", "recvProcs")], is_const=True)

        x.add_method("beginExchangeField", None, [refparam(intfield, "field")], is_const=True)
        x.add_method("beginExchangeField", None, [refparam(scalarfield, "field")], is_const=True)
        x.add_method("beginExchangeField", None, [refparam(vectorfield, "field")], is_const=True)
        x.add_method("beginExchangeField", None, [refparam(tensorfield, "field")], is_const=True)
        x.add_method("beginExchangeField", None, [refparam(symtensorfield, "field")], is_const=True)
        x.add_method("beginExchangeField", None, [refparam(thirdranktensorfield, "field")], is_const=True)
        x.add_method("beginExchangeFieldVariableSize", None, [refparam(vectordoublefield, "field")], is_const=True)
        x.add_method("beginExchangeFieldVariableSize", None, [refparam(vectorvectorfield, "field")], is_const=True)

        x.add_method("setAllGhostNodes", None, [refparam(database, "dataBase")], is_virtual=True, is_pure_virtual=True)
        x.add_method("finalizeGhostBoundary", None, [], is_const=True, is_virtual=True)
        x.add_method("finalizeExchanges", None, [])
        x.add_method("cullGhostNodes", None, [constrefparam(intfieldlist, "flagSet"),
                                              refparam(intfieldlist, "old2newIndexMap"),
                                              refparam("vector_of_int", "numNodesRemoved")], is_virtual=True)

        x.add_method("reset", None, [constrefparam(database, "dataBase")], visibility="protected", is_virtual=True)

        # x.add_function_as_method("getNestedGridNeighbor",
        #                          retval(ptr(nestedgridneighbor), reference_existing_object=True),
        #                          [param(me, "self"), constrefparam(nodelist, "nodeList")],
        #                          template_parameters = [dim],
        #                          custom_name = "getNestedGridNeighbor")

        # Override the Boundary abstract methods.
        generateBoundaryVirtualBindings(x, ndim, False)

        # Attributes.
        x.add_instance_attribute("domainID", "int", getter="domainID", is_const=True)
        x.add_instance_attribute("numDomains", "int", getter="numDomains", is_const=True)

        return

    #---------------------------------------------------------------------------
    # NestedGridDistributedBoundary
    #---------------------------------------------------------------------------
    def nestedGridDistributedBoundaryBindings(self, x, ndim):

        me = "Spheral::NestedGridDistributedBoundary%id" % ndim
        dim = "Spheral::Dim<%i> " % ndim
        vector = "Vector%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        domainboundarynodes = "DomainBoundaryNodes%id" % ndim
        intfield = "Spheral::IntField%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        tensorfield = "Spheral::TensorField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        thirdranktensorfield = "Spheral::ThirdRankTensorField%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        gridcellindex = "Spheral::GridCellIndex%id" % ndim
        vector_of_gridcellindex = "vector_of_GridCellIndex%id" % ndim
        vector_of_vector_of_gridcellindex = "vector_of_vector_of_GridCellIndex%id" % ndim

        # No constructors -- just the instance method since this is a singleton.
        x.add_method("instancePtr", retval(ptr(me), caller_owns_return=True), [], is_static=True, custom_name="instance")

        # Methods.
        x.add_method("setAllGhostNodes", None, [refparam(database, "dataBase")], is_virtual=True)
        x.add_method("reset", None, [refparam(database, "dataBase")], is_virtual=True)
        x.add_method("maxNumGridLevels", "int", [constrefparam(database, "dataBase")], is_const=True)
        x.add_method("setGridCellInfluenceRadius", "int", [constrefparam(database, "dataBase"),
                                                           param("int", "newGridCellInfluenceRadius")], is_const=True)
        x.add_method("flattenOccupiedGridCells", None, [constrefparam(database, "dataBase"),
                                                        refparam(vector_of_vector_of_gridcellindex, "gridCells")], is_const=True)
        x.add_method("packGridCellIndices", None, [constrefparam(vector_of_vector_of_gridcellindex, "gridCells"),
                                                    refparam("vector_of_int", "packedGridCells")], is_const=True)
        x.add_method("unpackGridCellIndices", None, [constrefparam("vector_of_int", "packedGridCells"),
                                                      constrefparam("vector_of_int", "gridCellDimensions"),
                                                      refparam(vector_of_vector_of_gridcellindex, "gridCells")], is_const=True)
        
        # Attributes.
        x.add_instance_attribute("boxCulling", "bool", getter="boxCulling", setter="boxCulling")
        x.add_instance_attribute("gridCellInfluenceRadius", "int", getter="gridCellInfluenceRadius", setter="gridCellInfluenceRadius")

        return

    #---------------------------------------------------------------------------
    # BoundingVolumeDistributedBoundary
    #---------------------------------------------------------------------------
    def boundingVolumeDistributedBoundaryBindings(self, x, ndim):

        me = "Spheral::BoundingVolumeDistributedBoundary%id" % ndim
        dim = "Spheral::Dim<%i> " % ndim
        database = "Spheral::DataBase%id" % ndim

        # No constructors -- just the instance method since this is a singleton.
        x.add_method("instancePtr", retval(ptr(me), caller_owns_return=True), [], is_static=True, custom_name="instance")

        # Methods.
        x.add_method("setAllGhostNodes", None, [refparam(database, "dataBase")], is_virtual=True)

        return

    #---------------------------------------------------------------------------
    # TreeDistributedBoundary
    #---------------------------------------------------------------------------
    def treeDistributedBoundaryBindings(self, x, ndim):

        me = "Spheral::TreeDistributedBoundary%id" % ndim
        dim = "Spheral::Dim<%i> " % ndim
        database = "Spheral::DataBase%id" % ndim

        # No constructors -- just the instance method since this is a singleton.
        x.add_method("instancePtr", retval(ptr(me), caller_owns_return=True), [], is_static=True, custom_name="instance")

        # Methods.
        x.add_method("setAllGhostNodes", None, [refparam(database, "dataBase")], is_virtual=True)

        return

    #---------------------------------------------------------------------------
    # DistributedBoundary::DomainBoundaryNodes
    #---------------------------------------------------------------------------
    def domainBoundaryNodesBindings(self, x, ndim):

        me = "DistributedBoundary%id::DomainBoundaryNodes%id" % (ndim, ndim)
        vector = "Vector%id" % ndim

        # Attributes.
        x.add_instance_attribute("sendNodes", "vector_of_int")
        x.add_instance_attribute("receiveNodes", "vector_of_int")

        return

    #---------------------------------------------------------------------------
    # DomainNode
    #---------------------------------------------------------------------------
    def domainNodeBindings(self, x, ndim):

        me = "DomainNode%id" % ndim
        vector = "Vector%id" % ndim

        # Constructors.
        x.add_constructor([])
        x.add_constructor([param(me, "rhs")])

        # Attributes.
        x.add_instance_attribute("localNodeID", "int")
        x.add_instance_attribute("uniqueLocalNodeID", "int")
        x.add_instance_attribute("globalNodeID", "int")
        x.add_instance_attribute("nodeListID", "int")
        x.add_instance_attribute("domainID", "int")
        x.add_instance_attribute("work", "double")
        x.add_instance_attribute("position", vector)

        # Methods
        x.add_binary_comparison_operator("==")
        x.add_binary_comparison_operator("!=")

        x.add_method("packSize", "int", [], is_static=True)
        x.add_method("pack", "vector_of_double", [], is_const=True)

        return

    #---------------------------------------------------------------------------
    # RedistributeNodes
    #---------------------------------------------------------------------------
    def redistributeNodesBindings(self, x, ndim):

        me = "RedistributeNodes%id" % ndim
        vector = "Vector%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        domainboundarynodes = "DomainBoundaryNodes%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        vector_of_domainnode = "vector_of_DomainNode%id" % ndim
        vector_of_boundary = "vector_of_Boundary%id" % ndim

        # Constructors.
        x.add_constructor([])

        # Methods.
        x.add_method("redistributeNodes", None, [refparam(database, "dataBase"),
                                                 param(vector_of_boundary, "boundaries", default_value="%s()" % vector_of_boundary)],
                     is_virtual = True,
                     is_pure_virtual = True)

        x.add_method("numGlobalNodes", "int", [constrefparam(database, "dataBase")], is_const=True)
        x.add_method("currentDomainDecomposition", vector_of_domainnode, [constrefparam(database, "dataBase"),
                                                                          constrefparam(intfieldlist, "globalIDs")], is_const=True)
        x.add_method("currentDomainDecomposition", vector_of_domainnode, [constrefparam(database, "dataBase"),
                                                                          constrefparam(intfieldlist, "globalIDs"),
                                                                          constrefparam(scalarfieldlist, "workPerNode")], is_const=True)
        x.add_method("enforceDomainDecomposition", None, [constrefparam(vector_of_domainnode, "nodeDistribution"),
                                                          refparam(database, "dataBase")], is_const=True)
        x.add_method("validDomainDecomposition", "bool", [constrefparam(vector_of_domainnode, "nodeDistribution"),
                                                          constrefparam(database, "dataBase")], is_const=True)
        x.add_method("workPerNode", scalarfieldlist, [constrefparam(database, "dataBase"), param("double", "Hextent")], is_const=True)
        x.add_method("gatherDomainDistributionStatistics", "std::string", [constrefparam(scalarfieldlist, "work")], is_const=True)

        # Attributes.
        x.add_instance_attribute("domainID", "int", getter="domainID", is_const=True)
        x.add_instance_attribute("numDomains", "int", getter="numDomains", is_const=True)
        x.add_instance_attribute("computeWork", "bool", getter="computeWork", setter="computeWork")

        return

    #---------------------------------------------------------------------------
    # DistributeByXPosition
    #---------------------------------------------------------------------------
    def distributeByXPositionBindings(self, x, ndim):

        database = "Spheral::DataBase%id" % ndim
        vector_of_boundary = "vector_of_Boundary%id" % ndim
        
        # Constructors.
        x.add_constructor([])

        # Methods.
        x.add_method("redistributeNodes", None, [refparam(database, "dataBase"),
                                                 param(vector_of_boundary, "boundaries", default_value="%s()" % vector_of_boundary)],
                     is_virtual = True)

        return


    #---------------------------------------------------------------------------
    # NestedGridRedistributeNodes
    #---------------------------------------------------------------------------
    def nestedGridRedistributeNodesBindings(self, x, ndim):

        me = "NestedGridRedistributeNodes%id" % ndim
        vector = "Vector%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        domainboundarynodes = "DomainBoundaryNodes%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        vector_of_domainnode = "vector_of_DomainNode%id" % ndim
        vector_of_boundary = "vector_of_Boundary%id" % ndim
        gridcell = "Spheral::GridCellIndex%id" % ndim

        # Constructors.
        x.add_constructor([param("double", "Hextent")])

        # Methods.
        x.add_method("redistributeNodes", None, [refparam(database, "dataBase"),
                                                 param(vector_of_boundary, "boundaries", default_value="%s()" % vector_of_boundary)],
                     is_virtual = True)

        x.add_method("setMasterNodeLists", None, [refparam(database, "dataBase"),
                                                  constrefparam(gridcell, "gridCell"),
                                                  param("int", "gridLevel"),
                                                  refparam("vector_of_vector_of_int", "masterLists"),
                                                  refparam("vector_of_vector_of_int", "coarseNeighbors")], is_const=True)
        x.add_method("gatherAvailableCoarseNodes", None, [constrefparam(database, "dataBase"),
                                                          constrefparam(vector_of_domainnode, "nodeDistribution"),
                                                          constrefparam(scalarfieldlist, "work"),
                                                          constrefparam("vector_of_vector_of_int", "localCoarseNeighbors"),
                                                          refparam("vector_of_int", "globalNodeIndices"),
                                                          refparam("vector_of_double", "globalNodeWork")], is_const=True)

        # Attributes.
        x.add_instance_attribute("Hextent", "double", getter="Hextent", setter="Hextent")

        return

    #---------------------------------------------------------------------------
    # SpaceFillingCurveRedistributeNodes
    #---------------------------------------------------------------------------
    def spaceFillingCurveRedistributeNodesBindings(self, x, ndim):

        me = "SpaceFillingCurveRedistributeNodes%id" % ndim
        vector = "Vector%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        domainboundarynodes = "DomainBoundaryNodes%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        vector_of_domainnode = "vector_of_DomainNode%id" % ndim
        vector_of_boundary = "vector_of_Boundary%id" % ndim
        gridcell = "Spheral::GridCellIndex%id" % ndim
        ullfieldlist = "Spheral::ULLFieldList%id" % ndim
        pair_vector_vector = "pair_Vector%id_Vector%id" % (ndim, ndim)
        vector_of_pair_ull_domainnode = "vector_of_pair_ULL_DomainNode%id" % ndim
        key = "uint64_t"
        
        # Constructors.
        x.add_constructor([param("double", "dummy"),
                           param("double", "minNodesPerDomainFraction"),
                           param("double", "maxNodesPerDomainFraction"),
                           param("bool", "workBalance", default_value="false"),
                           param("bool", "localReorderOnly", default_value="false")])

        # Methods.
        x.add_method("computeHashedIndices", ullfieldlist, [constrefparam(database, "dataBase")],
                     is_const = True,
                     is_virtual = True,
                     is_pure_virtual = True)

        x.add_method("redistributeNodes", None, [refparam(database, "dataBase"),
                                                 param(vector_of_boundary, "boundaries", default_value="%s()" % vector_of_boundary)],
                     is_virtual = True)
        x.add_method("computeStepSize", vector, [constrefparam(pair_vector_vector, "box")], is_const=True)
        x.add_method("buildIndex2IDPairs", vector_of_pair_ull_domainnode, [constrefparam(ullfieldlist, "indices"),
                                                                           constrefparam(vector_of_domainnode, "domainNodes")], is_const=True)
##         x.add_method("findUpperKey", key, [constrefparam("vector_of_ULL", "indices"),
##                                            constrefparam("vector_of_int", "count"),
##                                            constrefparam("vector_of_double", "work"),
##                                            param(key, "lowerBound"),
##                                            param(key, "maxUpperBound"),
##                                            param("double", "workTarget"),
##                                            param("int", "minNodes"),
##                                            param("int", "maxNodes"),
##                                            refparam(key, "upperKey"),
##                                            refparam("int", "numNodes")], is_const=True)
        x.add_method("numIndicesInRange", "int", [constrefparam("vector_of_ULL", "indices"),
                                                   constrefparam("vector_of_int", "count"),
                                                   param(key, "lowerBound"),
                                                   param(key, "upperBound")], is_const=True)
        x.add_method("workInRange", "double", [constrefparam("vector_of_ULL", "indices"),
                                               constrefparam("vector_of_double", "work"),
                                               param(key, "lowerBound"),
                                               param(key, "upperBound")], is_const=True)
        x.add_method("workAndNodesInRange", None, [constrefparam("vector_of_ULL", "indices"),
                                                   constrefparam("vector_of_int", "count"),
                                                   constrefparam("vector_of_double", "work"),
                                                   param(key, "lowerBound"),
                                                   param(key, "upperBound"),
                                                   refparam("int", "countInRange"),
                                                   refparam("double", "workInRange")], is_const=True)
        x.add_method("targetNumNodes", "int", [param("int", "numGlobal"),
                                               param("int", "numProcs"),
                                               param("int", "targetProc")], is_const=True)
        x.add_method("findNextIndex", key, [constrefparam("vector_of_ULL", "indices"),
                                            param(key, "index"),
                                            param(key, "maxIndex")], is_const=True)
        x.add_method("domainForIndex", "int", [param(key, "index"),
                                               constrefparam("vector_of_pair_ULL_ULL", "indexRanges")], is_const=True)

        # Attributes.
        x.add_instance_attribute("workBalance", "bool", getter="workBalance", setter="workBalance")
        x.add_instance_attribute("localReorderOnly", "bool", getter="localReorderOnly", setter="localReorderOnly")
        x.add_instance_attribute("minNodesPerDomainFraction", "double", getter="minNodesPerDomainFraction", setter="minNodesPerDomainFraction")
        x.add_instance_attribute("maxNodesPerDomainFraction", "double", getter="maxNodesPerDomainFraction", setter="maxNodesPerDomainFraction")

        return

    #---------------------------------------------------------------------------
    # MortonOrderRedistributeNodes, PeanoHilbertOrderRedistributeNodes
    #---------------------------------------------------------------------------
    def genericOrderRedistributeNodesBindings(self, x, ndim):

        vector = "Vector%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        domainboundarynodes = "DomainBoundaryNodes%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        vector_of_domainnode = "vector_of_DomainNode%id" % ndim
        vector_of_boundary = "vector_of_Boundary%id" % ndim
        gridcell = "Spheral::GridCellIndex%id" % ndim
        ullfieldlist = "Spheral::ULLFieldList%id" % ndim
        pair_vector_vector = "pair_Vector%id_Vector%id" % (ndim, ndim)
        vector_of_pair_ull_domainnode = "vector_of_pair_ULL_DomainNode%id" % ndim
        key = "uint64_t"
        
        # Constructors.
        x.add_constructor([param("double", "dummy"),
                           param("double", "minNodesPerDomainFraction", default_value="0.5"),
                           param("double", "maxNodesPerDomainFraction", default_value="1.5"),
                           param("bool", "workBalance", default_value="false"),
                           param("bool", "localReorderOnly", default_value="false")])

        # Methods.
        x.add_method("computeHashedIndices", ullfieldlist, [constrefparam(database, "dataBase")],
                     is_const = True,
                     is_virtual = True)

        return
    
    #---------------------------------------------------------------------------
    # SortAndDivideRedistributeNodes
    #---------------------------------------------------------------------------
    def genericSortAndDivideRedistributeNodesBindings(self, x, ndim):

        vector = "Vector%id" % ndim
        eigenstruct = "EigenStruct%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        domainboundarynodes = "DomainBoundaryNodes%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        vector_of_domainnode = "vector_of_DomainNode%id" % ndim
        vector_of_boundary = "vector_of_Boundary%id" % ndim

        # Constructors.
        x.add_constructor([param("double", "Hextent")])

        # Methods.
        x.add_method("redistributeNodes", None, [refparam(database, "dataBase"),
                                                 param(vector_of_boundary, "boundaries", default_value="%s()" % vector_of_boundary)],
                     is_virtual = True)
        x.add_method("shapeTensor", eigenstruct, [constrefparam(vector_of_domainnode, "domainNodes")], is_const=True)
        x.add_method("rotateIntoShapeTensorFrame", None, [constrefparam(eigenstruct, "shapeTensor"),
                                                          refparam(vector_of_domainnode, "domainNodes")], is_const=True)
        x.add_method("reduceDomainNodes", vector_of_domainnode, [constrefparam(vector_of_domainnode, "domainNodes"),
                                                                 param("int", "targetProc")], is_const=True)
        x.add_method("broadcastDomainNodes", vector_of_domainnode, [constrefparam(vector_of_domainnode, "domainNodes"),
                                                                    param("int", "sendProc")], is_const=True)

        # Attributes.
        x.add_instance_attribute("Hextent", "double", getter="Hextent", setter="Hextent")

        return

    def sortAndDivideRedistributeNodesBindings2d(self, x, ndim):

        eigenstruct = "EigenStruct%id" % ndim

        # Generic methods.
        self.genericSortAndDivideRedistributeNodesBindings(x, ndim)

        # Methods.
        x.add_method("domainsPerChunk", "vector_of_int", [constrefparam(eigenstruct, "shapeTensor")], is_const=True)

        return

    def sortAndDivideRedistributeNodesBindings3d(self, x, ndim):

        eigenstruct = "EigenStruct%id" % ndim

        # Generic methods.
        self.genericSortAndDivideRedistributeNodesBindings(x, ndim)

        # Methods.
        x.add_method("domainsPerChunk", "vector_of_vector_of_int", [constrefparam(eigenstruct, "shapeTensor")], is_const=True)

        return

    #---------------------------------------------------------------------------
    # VoronoiRedistributeNodes
    #---------------------------------------------------------------------------
    def voronoiRedistributeNodesBindings(self, x, ndim):

        vector = "Vector%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        domainboundarynodes = "DomainBoundaryNodes%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        vector_of_vector = "vector_of_%s" % vector
        vector_of_domainnode = "vector_of_DomainNode%id" % ndim
        vector_of_boundary = "vector_of_Boundary%id" % ndim
        gridcell = "Spheral::GridCellIndex%id" % ndim
        ullfieldlist = "Spheral::ULLFieldList%id" % ndim
        pair_vector_vector = "pair_Vector%id_Vector%id" % (ndim, ndim)
        vector_of_pair_ull_domainnode = "vector_of_pair_ULL_DomainNode%id" % ndim
        key = "uint64_t"
        
        # Constructors.
        x.add_constructor([param("double", "dummy"),
                           param("bool", "workBalance", default_value="false"),
                           param("bool", "balanceGenerators", default_value="true"),
                           param("double", "tolerance", default_value="1.0e-2"),
                           param("unsigned int", "maxIterations", default_value="200")])

        # Methods.
        x.add_method("redistributeNodes", None, [refparam(database, "dataBase"),
                                                 param(vector_of_boundary, "boundaries", default_value="%s()" % vector_of_boundary)],
                     is_virtual = True)
##         x.add_method("assignNodesToGenerators", None, [refparam(vector_of_domainnode, "nodes"),
##                                                        constrefparam(vector_of_vector, "generators"),
##                                                        refparam("vector_of_double", "generatorWork"),
##                                                        refparam(vector_of_vector, "newGenerators"),
##                                                        refparam("double", "maxDeltaGenerator"),
##                                                        refparam("double", "minWork"),
##                                                        refparam("double", "maxWork"),
##                                                        refparam("unsigned int", "minNodes"),
##                                                        refparam("unsigned int", "maxNodes")], is_const=True)

        # Attributes.
        x.add_instance_attribute("workBalance", "bool", getter="workBalance", setter="workBalance")
        x.add_instance_attribute("balanceGenerators", "bool", getter="balanceGenerators", setter="balanceGenerators")
        x.add_instance_attribute("tolerance", "double", getter="tolerance", setter="tolerance")
        x.add_instance_attribute("maxIterations", "int", getter="maxIterations", setter="maxIterations")

        return
    
