#-------------------------------------------------------------------------------
# TreeNeighbor
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Neighbor import *
from NeighborAbstractMethods import *

@PYB11template("Dimension")
class TreeNeighbor(Neighbor):

    PYB11typedefs = """
    typedef typename TreeNeighbor<%(Dimension)s>::LevelKey LevelKey;
    typedef typename TreeNeighbor<%(Dimension)s>::CellKey CellKey;
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
               kernelExtent = ("const double", "2.0"),
               xmin = "const Vector&",
               xmax = "const Vector&"):
        "Construct a TreeNeighbor"

    #...........................................................................
    # Methods
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

    @PYB11const
    def gridLevel(self, H="const SymTensor&"):
        "Find the tree level appropriate for H (units of 1/length)"
        return "unsigned"
                  
    @PYB11const
    def gridLevel(self, h="const double&"):
        "Find the tree level appropriate for h (units of length)"
        return "unsigned"
                  
    @PYB11const
    def dumpTree(self, globalTree="const bool"):
        "Return a dump of the tree structure as a string"
        return "std::string"

    @PYB11const
    def dumpTreeStatistics(self, globalTree="const bool"):
        "Return a string describing the overall statistics of the tree"
        return "std::string"
    
    @PYB11const
    def cellSize(self,
                 levelID = "const LevelKey"):
        "Cell size on the given level"
        return "double"

    @PYB11const
    def serialize(self, buffer="std::vector<char>&"):
        "Serialize the current tree state to a buffer"
        return "void"

    @PYB11const
    def nearestCellCenter(self,
                          xi = "const Vector&",
                          Hi = "const SymTensor&"):
        "Cell size on the given level"
        return "Vector"
                          
    @PYB11const
    def nearestCellCenter(self,
                          xi = "const Vector&",
                          hi = "const double"):
        "Cell size on the given level"
        return "Vector"
                          
    @PYB11const
    def occupied(self,
                 xi = "const Vector&",
                 Hi = "const SymTensor&"):
        """Test if the given position and extent is occupied, in the sense that
any tree cell is occupied at any coarser level down to it's native
level."""
        return "bool"

    @PYB11const
    def occupied(self,
                 xi = "const Vector&",
                 hi = "const double"):
        """Test if the given position and extent is occupied, in the sense that
any tree cell is occupied at any coarser level down to it's native
level."""
        return "bool"

    @PYB11const
    def setTreeMasterList(self,
                          levelID = "const LevelKey",
                          cellID = "const CellKey",
                          masterList = "std::vector<int>&",
                          coarseNeighbors = "std::vector<int>&",
                          ghostConnectivity = "const bool"):
        "For our parallel algorithm it is useful to be able to set the master/coarse information based on the given (level, cell)."
        return "void"

    @PYB11virtual
    @PYB11const
    def valid(self):
        "Test if the Neighbor is valid, i.e., ready to be queried for connectivity information."
        return "bool"

    #...........................................................................
    # Properties
    xmin = PYB11property("const Vector&", "xmin", doc="The minimum coordinate for the simulation bounding box")
    xmax = PYB11property("const Vector&", "xmax", doc="The maximum coordinate for the simulation bounding box")
    boxLength = PYB11property("double", "boxLength", doc="The maximum current cardinal coordinate distance across the bounding box")
    occupiedCells = PYB11property("std::vector<std::vector<CellKey>>", "occupiedCells", doc="The encoded cell key hashes for the cells currently occupied by nodes")

#-------------------------------------------------------------------------------
# Add the virtual interface
#-------------------------------------------------------------------------------
PYB11inject(NeighborAbstractMethods, TreeNeighbor, virtual=True)
