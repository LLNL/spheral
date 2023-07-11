#-------------------------------------------------------------------------------
# SortAndDivideRedistributeNodes2d
#-------------------------------------------------------------------------------
from PYB11Generator import *
from SortAndDivideRedistributeNodes import *

@PYB11template()
@PYB11template_dict({"Dimension": "Dim<2>"})
class SortAndDivideRedistributeNodes2d(SortAndDivideRedistributeNodes):
    """SortAndDivideRedistributeNodes2d -- 2-D implementation of the sort and 
divide algorithm for domain decomposition."""

    PYB11typedefs = """
    typedef typename KeyTraits::Key Key;
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               Hextent = "const double"):
        "Constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def redistributeNodes(self,
                          dataBase = "DataBase<%(Dimension)s>&",
                          boundaries = ("std::vector<Boundary<%(Dimension)s>*>", "std::vector<Boundary<%(Dimension)s>*>()")):
        """Given a Spheral++ data base of NodeLists, repartition it among the processors.
This is the method required of all descendent classes."""
        return "void"
