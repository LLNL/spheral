#-------------------------------------------------------------------------------
# NestedGridDistributedBoundary
#-------------------------------------------------------------------------------
from PYB11Generator import *

from DistributedBoundary import *

@PYB11template("Dimension")
#@PYB11singleton
class NestedGridDistributedBoundary(DistributedBoundary):
    """NestedGridDistributedBoundary -- Implementation of the Distributed Boundary
condition for use with NestedGridNeighbor based NodeLists."""

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    typedef typename %(Dimension)s::ThirdRankTensor ThirdRankTensor;
    typedef typename %(Dimension)s::FourthRankTensor FourthRankTensor;
    typedef typename %(Dimension)s::FifthRankTensor FifthRankTensor;
    typedef typename %(Dimension)s::FacetedVolume FacetedVolume;

    typedef typename DistributedBoundary<%(Dimension)s>::DomainBoundaryNodes DomainBoundaryNodes;
    typedef std::map<int, DomainBoundaryNodes> DomainBoundaryNodeMap;
    typedef std::map<NodeList<%(Dimension)s>*, DomainBoundaryNodeMap> NodeListDomainBoundaryNodeMap;
"""

    # The instance attribute.  We expose this as a property of the class.
    @PYB11static
    @PYB11cppname("instancePtr")
    @PYB11returnpolicy("take_ownership")
    def instance(self):
        return "NestedGridDistributedBoundary<%(Dimension)s>*"

    #...........................................................................
    # Methods
    @PYB11const
    def maxNumGridLevels(self,
                         dataBase = "const DataBase<%(Dimension)s>&"):
        "Determine the max number of occupied grid levels for all NodeLists in a DataBase."
        return "int"

    @PYB11const
    def setGridCellInfluenceRadius(self,
                                   dataBase = "DataBase<%(Dimension)s>&",
                                   newGridCellInfluenceRadius = "const int"):
        "Determine the radius of influence in gridcells being used."
        return "int"

    @PYB11const
    def flattenOccupiedGridCells(self,
                                 dataBase = "const DataBase<%(Dimension)s>&",
                                 gridCells = "std::vector< std::vector< GridCellIndex<%(Dimension)s> > >&"):
        "Build the set of occupied grid cells for all NodeLists on this process."
        return "void"

    @PYB11const
    def packGridCellIndices(self,
                            gridCellSet = "const std::vector< std::vector< GridCellIndex<%(Dimension)s> > >&",
                            packedGridCellIndices = "std::vector<int>&"):
        "Pack the occupied grid cell set into a C style array syntax for messaging with MPI."
        return "void"

    @PYB11const
    def unpackGridCellIndices(self,
                              packedGridCellIndices = "const std::vector<int>&",
                              gridCellDimension = "const std::vector<int>&",
                              gridCellSet = "std::vector< std::vector< GridCellIndex<%(Dimension)s> > >&"):
        "Unpack the occupied grid cell set from C style array syntax to the more sensible set of occupied grid cells."
        return "void"

    #...........................................................................
    # Virtual methods
    @PYB11pure_virtual
    def setAllGhostNodes(self,
                         dataBase = "DataBase<%(Dimension)s>&"):
        "Descendent Distributed Neighbors are required to provide the setGhostNodes method for DataBases."
        return "void"

    @PYB11virtual
    def reset(self,
              dataBase = "const DataBase<%(Dimension)s>&"):
        "Override the Boundary method for clearing the maps."
        return "void"

    #...........................................................................
    # Properties
    occupiedGridCells = PYB11property("const std::vector<std::vector<std::vector<GridCellIndex<%(Dimension)s>>>>&",
                                      returnpolicy="reference_internal",
                                      doc="The last set of occupied grid cells (only the master process knows the full set")
    boxCulling = PYB11property("bool", "boxCulling", "boxCulling")
    gridCellInfluenceRadius = PYB11property("int", "gridCellInfluenceRadius", "gridCellInfluenceRadius",
                                            doc="The grid cell influence radius to use for communicating.")
