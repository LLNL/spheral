#-------------------------------------------------------------------------------
# Neighbor
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11template("Dimension")
class Neighbor:

    typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    typedef NodeList<%(Dimension)s> NodeListType;
    typedef GeomPlane<%(Dimension)s> Plane;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               nodeList = "NodeListType&",
               searchType = "const NeighborSearchType",
               kernelExtent = "const double"):
        "Construct a Neighbor with minimal info"

    #...........................................................................
    # Methods
    @PYB11returnpolicy("reference_internal")
    @PYB11const
    def nodeExtentField(self):
        "The internal field of node extents"
        return "const Field<%(Dimension)s, Vector>&"

    @PYB11const
    def nodeList(self):
        "The NodeList associated with this Neighbor"
        return "const NodeListType&"

    @PYB11pycppname("nodeList")
    def setnodeList(self, val="NodeListType&"):
        "Set the NodeList for this Neighbor"
        return "void"

    def unregisterNodeList(self):
        "Disassociate with the currently registered NodeList."
        return "void"

    @PYB11const
    def nodeExtent(self, nodeID="int"):
        "Return the node extent for the requested node ID in the associated NodeList"
        return "Vector"
                           
    def setNodeExtents(self):
        "Set the node extent field based on the current NodeList"
        return "void"

    @PYB11pycppname("setNodeExtents")
    def setNodeExtents1(self, nodeIDs="const std::vector<int>&"):
        "Set the node extent field for the given nodes in the current NodeList"
        return "void"

    def setInternalNodeExtents(self):
        "Set the node extent field for the internal nodes in the current NodeList"
        return "void"

    def setGhostNodeExtents(self):
        "Set the node extent field for the ghost nodes in the current NodeList"
        return "void"

    @PYB11pycppname("setMasterList")
    @PYB11const
    def setMasterList0(self,
                       nodeID = "int",
                       masterList = "std::vector<int>&",
                       coarseNeighbors = "std::vector<int>&"):
        "Fill the given arrays with (master, coarse) neighbor info for the given node"
        return "void"

    @PYB11pycppname("setRefineNeighborList")
    @PYB11const
    def setRefineNeighborList0(self,
                               nodeID = "int",
                               coarseNeighbors = "const std::vector<int>&",
                               refineNeighbors = "std::vector<int>&"):
        "Fill the final array with refine neighbor info for the given node based on the coarse set"
        return "void"

    @PYB11const
    def precullList(self,
                    minMasterPosition = "const Vector&",
                    maxMasterPosition = "const Vector&",
                    minMasterExtent = "const Vector&",
                    maxMasterExtent = "const Vector&",
                    coarseNeighbors = "const std::vector<int>&"):
        "Return a culled list of potential neighbors based on the (min,max) info given"
        return "std::vector<int>"

    #...........................................................................
    # Abstract interface
    @PYB11pure_virtual
    @PYB11const
    @PYB11pycppname("setMasterList")
    def setMasterList1(self,
                       position = "const Vector&",
                       H = "const Scalar&",
                       masterList = "std::vector<int>&",
                       coarseNeighbors = "std::vector<int>&"):
        "Fill the given arrays with (master, coarse) neighbor info for the given (position, H)"
        return "void"

    @PYB11pure_virtual
    @PYB11const
    @PYB11pycppname("setMasterList")
    def setMasterList2(self,
                       position = "const Vector&",
                       H = "const SymTensor&",
                       masterList = "std::vector<int>&",
                       coarseNeighbors = "std::vector<int>&"):
        "Fill the given arrays with (master, coarse) neighbor info for the given (position, H)"
        return "void"

    @PYB11pure_virtual
    @PYB11const
    @PYB11pycppname("setRefineNeighborList")
    def setRefineNeighborList1(self,
                               position = "const Vector&",
                               H = "const Scalar&",
                               coarseNeighbors = "const std::vector<int>&",
                               refineList = "std::vector<int>&"):
        "Fill the given arrays with (coarse, refine) neighbor info for the given (position, H)"
        return "void"

    @PYB11pure_virtual
    @PYB11const
    @PYB11pycppname("setRefineNeighborList")
    def setRefineNeighborList2(self,
                               position = "const Vector&",
                               H = "const SymTensor&",
                               coarseNeighbors = "const std::vector<int>&",
                               refineList = "std::vector<int>&"):
        "Fill the given arrays with (coarse, refine) neighbor info for the given (position, H)"
        return "void"

    @PYB11pure_virtual
    @PYB11const
    @PYB11pycppname("setMasterList")
    def setMasterList3(self,
                       position = "const Vector&",
                       masterList = "std::vector<int>&",
                       coarseNeighbors = "std::vector<int>&"):
        "Fill the given arrays with (master, coarse) neighbor info for the given position"
        return "void"

    @PYB11pure_virtual
    @PYB11const
    @PYB11pycppname("setRefineNeighborList")
    def setRefineNeighobrList3(self,
                               position = "const Vector&",
                               coarseNeighbors = "const std::vector<int>&",
                               refineList = "std::vector<int>&"):
        "Fill the given arrays with (coarse, refine) neighbor info for the given position"
        return "void"

    @PYB11pure_virtual
    @PYB11const
    @PYB11pycppname("setMasterList")
    def setMasterList4(self,
                       enterPlane = "const Plane&",
                       exitPlane = "const Plane&",
                       masterList = "std::vector<int>&",
                       coarseNeighbors = "std::vector<int>&"):
        "Fill the given arrays with (master, coarse) neighbor info for the given (enter, exit) plane proximity"
        return "void"

    @PYB11pure_virtual
    def updateNodes(self):
        "Update the internal connectivity information based on the state of associated NodeList"
        return "void"

    @PYB11pure_virtual
    @PYB11pycppname("updateNodes")
    def updateNodes1(self,
                     nodeIDs = "const std::vector<int>&"):
        "Update the internal connectivity information for the given nodes based on the state of associated NodeList"
        return "void"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def reinitialize(self,
                     xmin = "const Vector&",
                     xmax = "const Vector&",
                     htarget = "const Scalar"):
        "Reinitialize to possibly more efficient based on the specified box (xmin,xmax) and htarget size"
        return "void"

    @PYB11virtual
    @PYB11const
    def valid(self):
        "Test if the Neighbor is valid, i.e., ready to be queried for connectivity information."
        return "bool"

    #...........................................................................
    # Static methods
    @PYB11static
    @PYB11pycppname("HExtent")
    def HExtent0(self,
                 H = "const Scalar&",
                 kernelExtent = "const double"):
        "Return the maximum extent in each Cartesian direction for the given H"
        return "Vector"

    @PYB11static
    @PYB11pycppname("HExtent")
    def HExtent1(self,
                 H = "const SymTensor&",
                 kernelExtent = "const double"):
        "Return the maximum extent in each Cartesian direction for the given H"
        return "Vector"

    #...........................................................................
    # Protected methods
    @PYB11implementation("&NeighborPublicist<%(Dimension)s>::accessNodeExtentField")
    def accessNodeExtentField(self):
        "Read/write access the protected node extent field member"
        return "Field<%(Dimension)s, Vector>&"

    #...........................................................................
    # Properties
    neighborSearchType = PYB11property("NeighborSearchType", "neighborSearchType", "neighborSearchType", doc="The search algorithm for nodes interacting")
    kernelExtent = PYB11property("double", "kernelExtent", "kernelExtent", doc="The kernel extent for nodes interacting in eta space")
