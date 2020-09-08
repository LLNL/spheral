#-------------------------------------------------------------------------------
# Neighbor
#-------------------------------------------------------------------------------
from PYB11Generator import *
from NeighborAbstractMethods import *

@PYB11template("Dimension")
class Neighbor:

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
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
    @PYB11virtual
    @PYB11const
    def setMasterList0(self,
                       nodeID = "int",
                       masterList = "std::vector<int>&",
                       coarseNeighbors = "std::vector<int>&",
                       ghostConnectivity = ("const bool", "false")):
        "Fill the given arrays with (master, coarse) neighbor info for the given node"
        return "void"

    @PYB11pycppname("setRefineNeighborList")
    @PYB11virtual
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
    # Virtual methods
    @PYB11virtual
    def reinitialize(self):
        "Reinitialize in case something changed."
        return "void"

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
    @PYB11protected
    def accessNodeExtentField(self):
        "Read/write access the protected node extent field member"
        return "Field<%(Dimension)s, Vector>&"

    #...........................................................................
    # Properties
    neighborSearchType = PYB11property("NeighborSearchType", "neighborSearchType", "neighborSearchType", doc="The search algorithm for nodes interacting")
    kernelExtent = PYB11property("double", "kernelExtent", "kernelExtent", doc="The kernel extent for nodes interacting in eta space")

#-------------------------------------------------------------------------------
# Add the abstract interface
#-------------------------------------------------------------------------------
PYB11inject(NeighborAbstractMethods, Neighbor, pure_virtual=True)
