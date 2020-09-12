#-------------------------------------------------------------------------------
# TreeDistributedBoundary
#-------------------------------------------------------------------------------
from PYB11Generator import *

from DistributedBoundary import *

@PYB11template("Dimension")
#@PYB11singleton
class TreeDistributedBoundary(DistributedBoundary):
    """TreeDistributedBoundary -- Implementation of the Distributed Boundary
condition for use with TreeNeighbor based NodeLists."""

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
        return "TreeDistributedBoundary<%(Dimension)s>*"

    #...........................................................................
    # Virtual methods
    @PYB11pure_virtual
    def setAllGhostNodes(self,
                         dataBase = "DataBase<%(Dimension)s>&"):
        "Descendent Distributed Neighbors are required to provide the setGhostNodes method for DataBases."
        return "void"
