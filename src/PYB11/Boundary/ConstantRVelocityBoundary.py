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

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    typedef typename %(Dimension)s::FacetedVolume FacetedVolume;
    typedef typename %(Dimension)s::ThirdRankTensor ThirdRankTensor;
    typedef typename %(Dimension)s::FourthRankTensor FourthRankTensor;
    typedef typename %(Dimension)s::FifthRankTensor FifthRankTensor;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               nodeList = "const NodeList<%(Dimension)s>&",
               nodeIndices = "const std::vector<int>&"):
        "Construct a constant radial velocity for the given nodes"

    #...........................................................................
    # Methods
    @PYB11virtual
    @PYB11const
    def enforceBoundary(self, field="Field<%(Dimension)s, Vector>&"):
        "Apply the boundary condition to the violation node values in the given Field."
        return "void"

#-------------------------------------------------------------------------------
# Inject restart methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, ConstantRVelocityBoundary)
