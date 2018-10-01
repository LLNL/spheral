#-------------------------------------------------------------------------------
# ConstantRVelocityBoundary
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Boundary import *
from ConstantVelocityBoundary import *
from RestartMethods import *

@PYB11template("Dimension")
class ConstantRVelocityBoundary(ConstantVelocityBoundary):
    """ConstantRVelocityBoundary -- A boundary condition to enforce a constant 
radial velocity component on a given set of nodes.

This boundary is very specialized -- it explicitly works on only one 
NodeList.
"""

    typedefs = """
    typedef %(Dimension)s DIM;
    typedef typename DIM::Scalar Scalar;
    typedef typename DIM::Vector Vector;
    typedef typename DIM::Tensor Tensor;
    typedef typename DIM::SymTensor SymTensor;
    typedef typename DIM::ThirdRankTensor ThirdRankTensor;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               nodeList = "const NodeList<DIM>&",
               nodeIndices = "const std::vector<int>&"):
        "Construct a constant radial velocity for the given nodes"

    #...........................................................................
    # Methods
    @PYB11virtual
    @PYB11const
    def enforceBoundary(self, field="Field<DIM, Vector>&"):
        "Apply the boundary condition to the violation node values in the given Field."
        return "void"

#-------------------------------------------------------------------------------
# Inject restart methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, ConstantRVelocityBoundary)
