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
    typedef GridCellIndex<%(Dimension)s> GridCellIndexType;
    typedef GeomPlane<%(Dimension)s> Plane;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               nodeList = "NodeListType&",
               searchType = ("const NeighborSearchType", "NeighborSearchType::GatherScatter"),
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
                       coarseNeighbors = "std::vector<int>&"):
        "Fill the given arrays with (master, coarse) neighbor info for the given (position, H)"
        return "void"

    @PYB11virtual
    @PYB11const
    @PYB11pycppname("setMasterList")
    def setMasterList2(self,
                       position = "const Vector&",
                       H = "const SymTensor&",
                       masterList = "std::vector<int>&",
                       coarseNeighbors = "std::vector<int>&"):
        "Fill the given arrays with (master, coarse) neighbor info for the given (position, H)"
        return "void"

    @PYB11virtual
    @PYB11const
    @PYB11pycppname("setRefineNeighborList")
    def setRefineNeighborList1(self,
                               position = "const Vector&",
                               H = "const Scalar&",
                               coarseNeighbors = "const std::vector<int>&",
                               refineNeighbors = "std::vector<int>&"):
        "Fill the given arrays with (coarse, refine) neighbor info for the given (position, H)"
        return "void"

    @PYB11virtual
    @PYB11const
    @PYB11pycppname("setRefineNeighborList")
    def setRefineNeighborList2(self,
                               position = "const Vector&",
                               H = "const SymTensor&",
                               coarseNeighbors = "const std::vector<int>&",
                               refineNeighbors = "std::vector<int>&"):
        "Fill the given arrays with (coarse, refine) neighbor info for the given (position, H)"
        return "void"

    @PYB11virtual
    @PYB11const
    @PYB11pycppname("setMasterList")
    def setMasterList3(self,
                       position = "const Vector&",
                       masterList = "std::vector<int>&",
                       coarseNeighbors = "std::vector<int>&"):
        "Fill the given arrays with (master, coarse) neighbor info for the given position"
        return "void"

    @PYB11virtual
    @PYB11const
    @PYB11pycppname("setRefineNeighborList")
    def setRefineNeighobrList3(self,
                               position = "const Vector&",
                               coarseNeighbors = "const std::vector<int>&",
                               refineNeighbors = "std::vector<int>&"):
        "Fill the given arrays with (coarse, refine) neighbor info for the given position"
        return "void"

    @PYB11virtual
    @PYB11const
    @PYB11pycppname("setMasterList")
    def setMasterList4(self,
                       enterPlane = "const Plane&",
                       exitPlane = "const Plane&",
                       masterList = "std::vector<int>&",
                       coarseNeighbors = "std::vector<int>&"):
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
        "Test if the Neighbor is valid, i.e., ready to be queried for connectivity information."
        return "bool"

    #...........................................................................
    # Methods
    @PYB11const
    def gridLevel(self, nodeID="const int"):
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
        return "const std::vector<std::vector<GridCellIndexType>>&"

    @PYB11returnpolicy("reference_internal")
    @PYB11pycppname("occupiedGridCells")
    @PYB11const
    def occupiedGridCells1(self, gridLevel="const int"):
        "The set of occupied gridcells on the given gridlevel"
        return "const std::vector<GridCellIndexType>&"

    def headOfGridCell(self):
        "Return the head of the chain for (grid cell, grid level)"

    def nextNodeInCell(self):
        "Find the next node in the chain from a given node"

    def internalNodesInCell(self):
        "Return a list of the internal nodes in the given (grid cell, grid level)"

    def nodesInCell(self):
        "Return a list of the nodes in the given (grid cell, grid level)"

    def appendNodesInCell(self):
        "Add to the chain of nodes for a given (grid cell, grid level)"

    def occupiedGridCellsInRange(self):
        "Find the occupied grid cells given (min, max) cells and grid level"

    def gridNormal(self):
        "Convert a coordinate vector to an integer normal"

    def mapGridCell(self):
        "Map a (grid cell, grid level) through a pair of planes"

    @PYB11const
    def setNestedMasterList(self,
                            gridCell = "const GridCellIndexType&",
                            gridLevel = "const int",
                            masterList = "std::vector<int>&",
                            coarseNeighbors = "std::vector<int>&"):
        "Worker method used to set master/coarse information"
        return "void"

    def findNestedNeighbors(self):
        "Return the neighbors for the given (grid cell, grid level)"
    
    #...........................................................................
    # Properties
    numGridLevels = PYB11property("int", "numGridLevels", "numGridLevels", doc="The maximum number of grid levels allowed")
    numOccupiedGridLevels = PYB11property("int", "numOccupiedGridLevels", doc="The number of grid levels populated by nodes")
    occupiedGridLevels = PYB11property("std::vector<int>", "occupiedGridLevels", doc="Array of the occupied grid levels")
    origin = PYB11property("const Vector&", "origin", "origin", doc="The origin for computing the GridCellIndex of a coordinate Vector")
    topGridSize = PYB11property("const double", "topGridSize", "topGridSize", doc="The cell size on the coarsest (top) grid level")
    gridCellInfluenceRadius = PYB11property("const int", "gridCellInfluenceRadius", "gridCellInfluenceRadius", doc="The radius in grid cells on a level a cell can interact with")
    gridCellSizeInv = PYB11property("const std::vector<double>&", "gridCellSizeInv", doc="The array of 1/grid cell size for each level")
    nodeInCell = PYB11property("const std::vector<std::vector<GridCellIndexType>>&", "nodeInCell", doc="The cell each node is in")
    masterGridLevel = PYB11property("int", "masterGridLevel", doc="The current master grid level")
    masterGridCellIndex = PYB11property("GridCellIndexType", "masterGridCellIndex", doc="The current master grid cell index")
    endOfLinkList = PYB11property("int", "endOfLinkList", doc="Value used to terminate a link list chain")
