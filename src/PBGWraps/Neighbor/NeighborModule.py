from pybindgen import *

from PBGutils import *
from ref_return_value import *

from CXXTypesModule import generateStdVectorBindings

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class Neighbor:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir, dims):

        self.dims = dims

        # Includes.
        mod.add_include('"%s/NeighborTypes.hh"' % srcdir)
    
        # Namespace.
        space = mod.add_cpp_namespace("Spheral")

        # Expose types.
        self.NeighborSearchType = space.add_enum("NeighborSearchType", [("None", "Spheral::NeighborSearchType::None"),
                                                                        ("Gather", "Spheral::NeighborSearchType::Gather"),
                                                                        ("Scatter", "Spheral::NeighborSearchType::Scatter"),
                                                                        ("GatherScatter", "Spheral::NeighborSearchType::GatherScatter")])

        for ndim in self.dims:
            exec("""
self.GridCellIndex%(ndim)id = addObject(space, "GridCellIndex%(ndim)id")
self.GridCellPlane%(ndim)id = addObject(space, "GridCellPlane%(ndim)id")
self.Neighbor%(ndim)id = addObject(space, "Neighbor%(ndim)id", allow_subclassing=True)
self.NestedGridNeighbor%(ndim)id = addObject(space, "NestedGridNeighbor%(ndim)id", parent=self.Neighbor%(ndim)id)
self.TreeNeighbor%(ndim)id = addObject(space, "TreeNeighbor%(ndim)id", parent=self.Neighbor%(ndim)id)
self.ConnectivityMap%(ndim)id = addObject(space, "ConnectivityMap%(ndim)id", allow_subclassing=True)
self.vector_of_GridCellIndex%(ndim)id = addObject(mod, "vector_of_GridCellIndex%(ndim)id", allow_subclassing=True)
self.vector_of_vector_of_GridCellIndex%(ndim)id = addObject(mod, "vector_of_vector_of_GridCellIndex%(ndim)id", allow_subclassing=True)
""" % {"ndim" : ndim})

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        for ndim in self.dims:
            exec("""
self.generateGridCellIndexBindings(self.GridCellIndex%(ndim)id, %(ndim)i)
self.generateGridCellPlaneBindings(self.GridCellPlane%(ndim)id, %(ndim)i)
self.generateNeighborBindings(self.Neighbor%(ndim)id, %(ndim)i)
self.generateNestedGridNeighborBindings(self.NestedGridNeighbor%(ndim)id, %(ndim)i)
self.generateTreeNeighborBindings(self.TreeNeighbor%(ndim)id, %(ndim)i)
self.generateConnectivityMapBindings(self.ConnectivityMap%(ndim)id, %(ndim)i)
generateStdVectorBindings(self.vector_of_GridCellIndex%(ndim)id, "Spheral::GridCellIndex%(ndim)id", "vector_of_GridCellIndex%(ndim)id", indexAsPointer=True)
generateStdVectorBindings(self.vector_of_vector_of_GridCellIndex%(ndim)id, "vector_of_GridCellIndex%(ndim)id", "vector_of_vector_of_GridCellIndex%(ndim)id", indexAsPointer=True)
""" % {"ndim" : ndim})

        return

    #---------------------------------------------------------------------------
    # GridCellIndex
    #---------------------------------------------------------------------------
    def generateGridCellIndexBindings(self, x, ndim):

        # Objects.
        me = "Spheral::GridCellIndex%id" % ndim
        gridcellplane = "Spheral::GridCellPlane%id" % ndim

        # Constructors.
        x.add_constructor([])
        x.add_constructor([param("int", "xIndex")])
        x.add_constructor([param("int", "xIndex"),
                           param("int", "yIndex")])
        x.add_constructor([param("int", "xIndex"),
                           param("int", "yIndex"),
                           param("int", "zIndex")])
        x.add_constructor([constrefparam(me, "rhs")])

        # Methods.
        x.add_method("setIndices", None, [param("int", "xIndex")])
        x.add_method("setIndices", None, [param("int", "xIndex"),
                                           param("int", "yIndex")])
        x.add_method("setIndices", None, [param("int", "xIndex"),
                                           param("int", "yIndex"),
                                           param("int", "zIndex")])
        x.add_method("operator()", "int", [param("int", "index")], is_const=True, custom_name="__call__")
        x.add_method("operator-", me, [], is_const=True, custom_name="__neg__")
        x.add_method("dot", "int", [constrefparam(me, "rhs")], is_const=True)
        x.add_method("compare", "int", [constrefparam(me, "rhs")], is_const=True)
        x.add_method("inRange", "bool", [constrefparam(me, "minGridCell"), constrefparam(me, "maxGridCell")], is_const=True)
        x.add_method("magnitude", "double", [], is_const=True)
        x.add_method("minElement", "int", [], is_const=True)
        x.add_method("maxElement", "int", [], is_const=True)
        x.add_method("sumElements", "int", [], is_const=True)
        x.add_method("productElements", "int", [], is_const=True)
        x.add_method("indexMin", "int", [], is_const=True)
        x.add_method("indexMax", "int", [], is_const=True)

        x.add_binary_numeric_operator("+")
        x.add_binary_numeric_operator("-")

        x.add_binary_numeric_operator("+", right="int")
        x.add_binary_numeric_operator("-", right="int")
        x.add_binary_numeric_operator("*", right="int")
        x.add_binary_numeric_operator("/", right="int")

        x.add_binary_comparison_operator("==")
        x.add_binary_comparison_operator("!=")
        x.add_binary_comparison_operator("<")
        x.add_binary_comparison_operator(">")
        x.add_binary_comparison_operator("<=")
        x.add_binary_comparison_operator(">=")

##         x.add_binary_comparison_operator("<", right_cppclass=gridcellplane)
##         x.add_binary_comparison_operator(">", right_cppclass=gridcellplane)
##         x.add_binary_comparison_operator("<=", right_cppclass=gridcellplane)
##         x.add_binary_comparison_operator(">=", right_cppclass=gridcellplane)

        # Attributes.
        x.add_instance_attribute("xIndex", "int", getter="xIndex", setter="xIndex")
        x.add_instance_attribute("yIndex", "int", getter="yIndex", setter="yIndex")
        x.add_instance_attribute("zIndex", "int", getter="zIndex", setter="zIndex")

        return

    #---------------------------------------------------------------------------
    # GridCellPlane
    #---------------------------------------------------------------------------
    def generateGridCellPlaneBindings(self, x, ndim):

        # Objects.
        me = "Spheral::GridCellPlane%id" % ndim
        gridcellindex = "Spheral::GridCellIndex%id" % ndim

        # Constructors.
        x.add_constructor([])
        x.add_constructor([constrefparam(me, "rhs")])
        x.add_constructor([constrefparam(gridcellindex, "point"),
                           constrefparam(gridcellindex, "normal")])

        # Methods.
        x.add_method("minimumDistance", "double", [constrefparam(gridcellindex, "point")], is_const=True)
        x.add_method("coplanar", "bool", [constrefparam(gridcellindex, "point")], is_const=True)
        x.add_method("parallel", "bool", [constrefparam(me, "point")], is_const=True)
        x.add_method("valid", "bool", [], is_const=True)

        x.add_binary_comparison_operator("==")
        x.add_binary_comparison_operator("!=")

##         x.add_binary_comparison_operator(">", right_cppclass=gridcellindex)
##         x.add_binary_comparison_operator("<", right_cppclass=gridcellindex)
##         x.add_binary_comparison_operator(">=", right_cppclass=gridcellindex)
##         x.add_binary_comparison_operator("<=", right_cppclass=gridcellindex)

        # Attributes.
        x.add_instance_attribute("point", gridcellindex, getter="point", setter="setPoint")
        x.add_instance_attribute("normal", gridcellindex, getter="normal", setter="setNormal")

        return

    #---------------------------------------------------------------------------
    # Neighbor
    #---------------------------------------------------------------------------
    def generateNeighborBindings(self, x, ndim):

        # Objects.
        me = "Spheral::Neighbor%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        gridcellplane = "Spheral::GridCellPlane%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        plane = "Plane%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim

        # Constructors.
        x.add_constructor([refparam(nodelist, "nodeList"), 
                           param("NeighborSearchType", "searchType"),
                           param("double", "kernelExtent")])

        # Methods.
        const_ref_return_value(x, me, "%s::nodeExtentField" % me, vectorfield, [], "nodeExtentField")
        const_ref_return_value(x, me, "%s::nodeList" % me, nodelist, [], "nodeList")
        x.add_method("nodeList", None, [refparam(nodelist, "nodeList")], custom_name="nodeList")
        x.add_method("unregisterNodeList", None, [])
        x.add_method("nodeExtent", vector, [param("int", "nodeID")], is_const=True)
        x.add_method("setNodeExtents", None, [])
        x.add_method("setNodeExtents", None, [constrefparam("vector_of_int", "nodeIDs")])
        x.add_method("setInternalNodeExtents", None, [])
        x.add_method("setGhostNodeExtents", None, [])
        x.add_method("setMasterList", None, [param("int", "nodeID"),
                                             refparam("vector_of_int", "masterList"),
                                             refparam("vector_of_int", "coarseNeighbors")], is_const=True, is_virtual=True)
        x.add_method("setRefineNeighborList", None, [param("int", "nodeID"),
                                                     constrefparam("vector_of_int", "coarseNeighbors"),
                                                     refparam("vector_of_int", "refineNeighbors")], is_const=True, is_virtual=True)
        x.add_method("precullList", "vector_of_int", [constrefparam(vector, "minMasterPosition"),
                                                      constrefparam(vector, "maxMasterPosition"),
                                                      constrefparam(vector, "minMasterExtent"),
                                                      constrefparam(vector, "maxMasterExtent"),
                                                      constrefparam("vector_of_int", "coarseList")], is_const=True)
        x.add_method("reinitialize", None,
                     [constrefparam(vector, "xmin"), constrefparam(vector, "xmax"), param("double", "htarget")],
                     is_virtual = True)
        x.add_method("valid", "bool", [], is_const=True, is_virtual=True)
        x.add_method("HExtent", vector, [param("double", "H"), param("double", "kernelExtent")], is_static=True)
        x.add_method("HExtent", vector, [constrefparam(symtensor, "H"), param("double", "kernelExtent")], is_static=True)

        # Attributes.
        x.add_instance_attribute("neighborSearchType", "NeighborSearchType", getter="neighborSearchType", setter="neighborSearchType")
        x.add_instance_attribute("kernelExtent", "double", getter="kernelExtent", setter="kernelExtent")

        # Add the pure virtual methods.
        self.generateNeighborVirtualBindings(x, ndim, True)

        return

    #---------------------------------------------------------------------------
    # Neighbor virtual methods
    #---------------------------------------------------------------------------
    def generateNeighborVirtualBindings(self, x, ndim, pureVirtual):

        # Objects.
        dim = "Spheral::Dim<%i>" % ndim
        gridcellplane = "Spheral::GridCellPlane%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        plane = "Plane%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim

        # Methods.
        x.add_method("setMasterList", None, [constrefparam(vector, "position"), constrefparam("double", "H"),
                                             refparam("vector_of_int", "masterList"),
                                             refparam("vector_of_int", "coarseNeighbors")], is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
        x.add_method("setMasterList", None, [constrefparam(vector, "position"), constrefparam(symtensor, "H"),
                                             refparam("vector_of_int", "masterList"),
                                             refparam("vector_of_int", "coarseNeighbors")], is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
        x.add_method("setRefineNeighborList", None, [constrefparam(vector, "position"), constrefparam("double", "H"),
                                                     constrefparam("vector_of_int", "coarseNeighbors"),
                                                     refparam("vector_of_int", "refineNeighbors")], is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
        x.add_method("setRefineNeighborList", None, [constrefparam(vector, "position"), constrefparam(symtensor, "H"),
                                                     constrefparam("vector_of_int", "coarseNeighbors"),
                                                     refparam("vector_of_int", "refineNeighbors")], is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
        x.add_method("setMasterList", None, [constrefparam(vector, "position"),
                                             refparam("vector_of_int", "masterList"),
                                             refparam("vector_of_int", "coarseNeighbors")], is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
        x.add_method("setRefineNeighborList", None, [constrefparam(vector, "position"),
                                                     constrefparam("vector_of_int", "coarseNeighbors"),
                                                     refparam("vector_of_int", "refineNeighbors")], is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
        x.add_method("setMasterList", None, [constrefparam(plane, "enterPlane"), constrefparam(plane, "exitPlane"), 
                                             refparam("vector_of_int", "masterList"),
                                             refparam("vector_of_int", "coarseNeighbors")], is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
        x.add_method("updateNodes", None, [], is_virtual=True, is_pure_virtual=pureVirtual)
        x.add_method("updateNodes", None, [constrefparam("vector_of_int", "nodeIDs")], is_virtual=True, is_pure_virtual=pureVirtual)

        return

    #---------------------------------------------------------------------------
    # NestedGridNeighbor
    #---------------------------------------------------------------------------
    def generateNestedGridNeighborBindings(self, x, ndim):

        # Objects.
        me = "Spheral::NestedGridNeighbor%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        gridcellplane = "Spheral::GridCellPlane%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        plane = "Plane%id" % ndim
        gridcellindex = "Spheral::GridCellIndex%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        vector_of_gridcell = "vector_of_GridCellIndex%id" % ndim
        vector_of_vector_of_gridcell = "vector_of_vector_of_GridCellIndex%id" % ndim

        # Constructors.
        x.add_constructor([refparam(nodelist, "nodeList"),
                           param("NeighborSearchType", "searchType", default_value="Spheral::NeighborSearchType::GatherScatter"),
                           param("int", "numGridLevels", default_value="31"),
                           param("double", "topGridCellSize", default_value="100.0"),
                           param(vector, "origin", default_value="%s()" % vector),
                           param("double", "kernelExtent", default_value="2.0"),
                           param("int", "gridCellInfluenceRadius", default_value="1")])

        # Methods.
        x.add_method("gridLevel", "int", [param("int", "nodeID")], is_const=True)
        x.add_method("gridCellIndex", gridcellindex, [param("int", "nodeID"), param("int", "gridLevel")], is_const=True)
        x.add_method("gridCellIndex", gridcellindex, [constrefparam(vector, "position"), param("int", "gridLevel")], is_const=True)
        x.add_method("translateGridCellRange", None, [constrefparam(gridcellindex, "gridCellMin"),
                                                      constrefparam(gridcellindex, "gridCellMax"),
                                                      param("int", "gridLevel"),
                                                      param("int", "targetGridLevel"),
                                                      refparam(gridcellindex, "targetMin"),
                                                      refparam(gridcellindex, "targetMax")], is_const=True)
##         x.add_method("cellInTree", "bool", [constrefparam(gridcellindex, "gridCell"), param("int", "gridLevel")], is_const=True)
        x.add_method("cellOccupied", "bool", [constrefparam(gridcellindex, "gridCell"), param("int", "gridLevel")], is_const=True)
##         x.add_method("daughterCells", vector_of_gridcell, [constrefparam(gridcellindex, "gridCell"), param("int", "gridLevel")], is_const=True)
        x.add_method("occupiedGridCells", vector_of_vector_of_gridcell, [], is_const=True)
        x.add_method("occupiedGridCells", vector_of_gridcell, [param("int", "gridLevel")], is_const=True)
        x.add_method("headOfGridCell", "int", [constrefparam(gridcellindex, "gridCell"), param("int", "gridLevel")], is_const=True)
        x.add_method("nextNodeInCell", "int", [param("int", "nodeID")], is_const=True)
        x.add_method("internalNodesInCell", "vector_of_int", [constrefparam(gridcellindex, "gridCell"),
                                                              param("int", "gridLevel")], is_const=True)
        x.add_method("nodesInCell", "vector_of_int", [constrefparam(gridcellindex, "gridCell"),
                                                      param("int", "gridLevel")], is_const=True)
        x.add_method("appendNodesInCell", None, [constrefparam(gridcellindex, "gridCell"),
                                                 param("int", "gridLevel"),
                                                 refparam("vector_of_int", "nodes")], is_const=True)
        x.add_method("occupiedGridCellsInRange", None, [refparam(vector_of_gridcell, "gridCells"),
                                                        constrefparam(gridcellindex, "minGridCell"),
                                                        constrefparam(gridcellindex, "maxGridCell"),
                                                        param("int", "gridLevel")], is_const=True)
        x.add_method("gridNormal", gridcellindex, [constrefparam(vector, "normal")], is_const=True)
        x.add_method("mapGridCell", gridcellindex, [constrefparam(gridcellindex, "gridCell"),
                                                    param("int", "gridLevel"),
                                                    constrefparam(plane, "enterPlane"),
                                                    constrefparam(plane, "exitPlane")], is_const=True)
        x.add_method("valid", "bool", [], is_const=True, is_virtual=True)
        x.add_method("setNestedMasterList", None, [constrefparam(gridcellindex, "gridCell"), param("int", "gridLevel"),
                                                   refparam("vector_of_int", "masterList"),
                                                   refparam("vector_of_int", "coarseNeighbors")])
        x.add_method("findNestedNeighbors", "vector_of_int", [constrefparam(gridcellindex, "gridCell"), param("int", "gridLevel")], is_const=True)

        # Attributes.
        x.add_instance_attribute("numGridLevels", "int", getter="numGridLevels", setter="numGridLevels")
##         x.add_instance_attribute("firstParentGridLevel", "int", getter="firstParentGridLevel", is_const=True)
        x.add_instance_attribute("numOccupiedGridLevels", "int", getter="numOccupiedGridLevels", is_const=True)
        x.add_instance_attribute("occupiedGridLevels", "vector_of_int", getter="occupiedGridLevels", is_const=True)
        x.add_instance_attribute("origin", vector, getter="origin", setter="origin")
        x.add_instance_attribute("topGridSize", "double", getter="topGridSize", setter="topGridSize")
        x.add_instance_attribute("gridCellSizeInv", "vector_of_double", getter="gridCellSizeInv", is_const=True)
        x.add_instance_attribute("nodeInCell", vector_of_vector_of_gridcell, getter="nodeInCell", is_const=True)
        x.add_instance_attribute("endOfLinkList", "int", getter="endOfLinkList", is_const=True)

        # Add the base virtual methods.
        self.generateNeighborVirtualBindings(x, ndim, False)

        return
    
    #---------------------------------------------------------------------------
    # TreeNeighbor
    #---------------------------------------------------------------------------
    def generateTreeNeighborBindings(self, x, ndim):

        # Objects.
        me = "Spheral::TreeNeighbor%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        plane = "Plane%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim

        # Constructors.
        x.add_constructor([refparam(nodelist, "nodeList"),
                           param("NeighborSearchType", "searchType"),
                           param("double", "kernelExtent"),
                           constrefparam(vector, "xmin"),
                           constrefparam(vector, "xmax")])

        # Methods.
        x.add_method("gridLevel", "unsigned int", [param("double", "h")], is_const=True)
        x.add_method("gridLevel", "unsigned int", [param(symtensor, "H")], is_const=True)
        x.add_method("dumpTree", "std::string", [param("bool", "globalTree")], is_const=True)
        x.add_method("dumpTreeStatistics", "std::string", [param("bool", "globalTree")], is_const=True)
        x.add_method("valid", "bool", [], is_const=True)

        # Attributes.
        x.add_instance_attribute("xmin", vector, getter="xmin", is_const=True)
        x.add_instance_attribute("xmax", vector, getter="xmax", is_const=True)
        x.add_instance_attribute("boxLength", "double", getter="boxLength", is_const=True)

        # Add the base virtual methods.
        self.generateNeighborVirtualBindings(x, ndim, False)

        return
    
    #---------------------------------------------------------------------------
    # ConnectivityMap
    #---------------------------------------------------------------------------
    def generateConnectivityMapBindings(self, x, ndim):

        # Objects.
        me = "Spheral::ConnectivityMap%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        gridcellindex = "Spheral::GridCellIndex%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        fluidnodelist = "Spheral::FluidNodeList%id" % ndim
        vector_of_const_nodelist = "vector_of_const_nodelist%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim

        # Constructor.
        x.add_constructor([])

        # Methods.
        x.add_method("patchConnectivity", None, [constrefparam(intfieldlist, "flags"),
                                                 constrefparam(intfieldlist, "old2new")])
        x.add_method("connectivityForNode", "vector_of_vector_of_int", [param(ptr("const " + nodelist), "nodeList"), param("int", "nodeID")], is_const=True)
        x.add_method("connectivityForNode", "vector_of_vector_of_int", [param("int", "nodeListID"), param("int", "nodeID")], is_const=True)
        x.add_method("connectivityIntersectionForNodes", "vector_of_vector_of_int", 
                     [param("int", "nodeListi"), param("int", "i"),
                      param("int", "nodeListj"), param("int", "j")], is_const=True)
        x.add_method("connectivityUnionForNodes", "vector_of_vector_of_int", 
                     [param("int", "nodeListi"), param("int", "i"),
                      param("int", "nodeListj"), param("int", "j")], is_const=True)
        x.add_method("numNeighborsForNode", "int", [param(ptr("const " + nodelist), "nodeList"), param("int", "nodeID")], is_const=True)
        x.add_method("calculatePairInteraction", "bool", [param("int", "nodeListi"),
                                                          param("int", "i"),
                                                          param("int", "nodeListj"),
                                                          param("int", "j"),
                                                          param("int", "firstGhostNodej")], is_const=True)
        x.add_function_as_method("nodeListFromConnectivityMap",
                                 retval(ptr(nodelist), reference_existing_object=True),
                                 [param(me, "self"), param("int", "index")],
                                 template_parameters = [dim],
                                 custom_name = "nodeList")
        x.add_function_as_method("numNodeListsInConnectivityMap",
                                 "int",
                                 [param(me, "self")],
                                 template_parameters = [dim],
                                 custom_name = "numNodeLists")
        x.add_method("numNodes", "int", [param("int", "nodeList")], is_const=True)
        x.add_method("ithNode", "int", [param("int", "nodeList"), param("int", "index")], is_const=True)
        x.add_method("valid", "bool", [], is_const=True)

        # Attributes.
        x.add_instance_attribute("buildGhostConnectivity", "bool", getter="buildGhostConnectivity", is_const=True)

        return
