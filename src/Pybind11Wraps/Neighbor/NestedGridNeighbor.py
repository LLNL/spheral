#-------------------------------------------------------------------------------
# NestedGridNeighbor
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Neighbor import *

@PYB11template("Dimension")
class NestedGridNeighbor(Neighbor):

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
               searchType = ("const NestedGridNeighborSearchType", "NeighborSearchType::GatherScatter"),
               numGridLevels = ("int", "31"),
               topGridCellSize = ("double", "100.0"),
               origin = ("Vector", "Vector::zero"),
               kernelExtent = ("const double", "2.0"),
               gridCellInfluenceRadius = ("int", "1")):
        "Construct a NestedGridNeighbor"

    #...........................................................................
    # Virtual interface
    @PYB11virtual
    @PYB11const
    @PYB11pycppname("setMasterList")
    def setMasterList1(self,
                       position = "const Vector&",
                       H = "const Scalar&",
                       masterList = "std::vector<int>&",
                       coarseNestedGridNeighbors = "std::vector<int>&"):
        "Fill the given arrays with (master, coarse) neighbor info for the given (position, H)"
        return "void"

    @PYB11virtual
    @PYB11const
    @PYB11pycppname("setMasterList")
    def setMasterList2(self,
                       position = "const Vector&",
                       H = "const SymTensor&",
                       masterList = "std::vector<int>&",
                       coarseNestedGridNeighbors = "std::vector<int>&"):
        "Fill the given arrays with (master, coarse) neighbor info for the given (position, H)"
        return "void"

    @PYB11virtual
    @PYB11const
    @PYB11pycppname("setRefineNestedGridNeighborList")
    def setRefineNestedGridNeighborList1(self,
                               position = "const Vector&",
                               H = "const Scalar&",
                               coarseNestedGridNeighbors = "const std::vector<int>&",
                               refineList = "std::vector<int>&"):
        "Fill the given arrays with (coarse, refine) neighbor info for the given (position, H)"
        return "void"

    @PYB11virtual
    @PYB11const
    @PYB11pycppname("setRefineNestedGridNeighborList")
    def setRefineNestedGridNeighborList2(self,
                               position = "const Vector&",
                               H = "const SymTensor&",
                               coarseNestedGridNeighbors = "const std::vector<int>&",
                               refineList = "std::vector<int>&"):
        "Fill the given arrays with (coarse, refine) neighbor info for the given (position, H)"
        return "void"

    @PYB11virtual
    @PYB11const
    @PYB11pycppname("setMasterList")
    def setMasterList3(self,
                       position = "const Vector&",
                       masterList = "std::vector<int>&",
                       coarseNestedGridNeighbors = "std::vector<int>&"):
        "Fill the given arrays with (master, coarse) neighbor info for the given position"
        return "void"

    @PYB11virtual
    @PYB11const
    @PYB11pycppname("setRefineNestedGridNeighborList")
    def setRefineNeighobrList3(self,
                               position = "const Vector&",
                               coarseNestedGridNeighbors = "const std::vector<int>&",
                               refineList = "std::vector<int>&"):
        "Fill the given arrays with (coarse, refine) neighbor info for the given position"
        return "void"

    @PYB11virtual
    @PYB11const
    @PYB11pycppname("setMasterList")
    def setMasterList4(self,
                       enterPlane = "const Plane&",
                       exitPlane = "const Plane&",
                       masterList = "std::vector<int>&",
                       coarseNestedGridNeighbors = "std::vector<int>&"):
        "Fill the given arrays with (master, coarse) neighbor info for the given (enter, exit) plane proximity"
        return "void"

    @PYB11virtual
    def updateNodes(self):
        "Update the internal connectivity information based on the state of associated NodeList"
        return "void"

    @PYB11virtual
    @PYB11pycppname("updateNodes")
    def updateNodes1(self,
                     nodeIDs = "const std::vector<int>&"):
        "Update the internal connectivity information for the given nodes based on the state of associated NodeList"
        return "void"

    @PYB11virtual
    @PYB11const
    def valid(self):
        "Test if the NestedGridNeighbor is valid, i.e., ready to be queried for connectivity information."
        return "bool"

    #...........................................................................
    # Methods
    @PYB11const
    def gridLevel(self, nodeID="const int", gridLevel="const int"):
        "Find the gridlevel for the given nodeID"
        return "int"

    @PYB11pycppname("gridLevel")
    @PYB11const
    def gridLevel1(self, H="const SymTensor&"):
        "Find the gridlevel for the given H"
        return "int"

    @PYB11pycppname("gridLevel")
    @PYB11const
    def gridLevel2(self, H="const Scalar&"):
        "Find the gridlevel for the given H"
        return "int"

    @PYB11const
    def gridCellIndex(self, nodeID="const int", gridLevel="const int"):
        "Find the GridCellIndex for the given node on the given level"
        return "GridCellIndexType"

    @PYB11pycppname("gridCellIndex")
    @PYB11const
    def gridCellIndex1(self, position="const Vector&", gridLevel="const int"):
        "Find the GridCellIndex for the given position on the given level"
        return "GridCellIndexType"

    def translateGridCellRange(self):
        return
                           
    def cellOccupied(self):
        "Test if the given (grid cell, grid level) is occupied"
        return

    @PYB11returnpolicy("reference_internal")
    @PYB11const
    def occupiedGridCells(self):
        "The full set of occupied gridcells on all gridlevels"
        return "std::vector<std::vector<int>>&"

    @PYB11returnpolicy("reference_internal")
    @PYB11pycppname("occupiedGridCells")
    @PYB11const
    def occupiedGridCells1(self, gridLevel="const int"):
        "The set of occupied gridcells on the given gridlevel"
        return "std::vector<int>&"

    #...........................................................................
    # Properties
    @PYB11ignore
    @PYB11pycppname("numGridLevels")
    @PYB11const
    def getnumGridLevels(self):
        return "int"

    @PYB11ignore
    @PYB11pycppname("numGridLevels")
    def setnumGridLevels(self, val="const int"):
        return "void"

    @PYB11ignore
    @PYB11pycppname("numOccupiedGridLevels")
    @PYB11const
    def getnumOccupiedGridLevels(self):
        return "int"

    @PYB11ignore
    @PYB11pycppname("occupiedGridLevels")
    @PYB11const
    def getoccupiedGridLevels(self):
        return "std::vector<int>"

    numGridLevels = property(getnumGridLevels, setnumGridLevels, doc="The number of gridlevels")
    numOccupiedGridLevels = property(getnumOccupiedGridLevels, doc="The number of occupied gridlevels")
    occupiedGridLevels = property(getoccupiedGridLevels, doc="The occupied grid levels")
