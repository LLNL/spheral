#-------------------------------------------------------------------------------
# ConstantXVelocityBoundary
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Boundary import *
from ConstantVelocityBoundary import *

@PYB11template("Dimension")
class ConstantXVelocityBoundary(ConstantVelocityBoundary):
    """ConstantXVelocityBoundary -- A boundary condition to enforce a constant 
x-component velocity on a given set of nodes.

This boundary is very specialized -- it explicitly works on only one 
NodeList.
"""

    typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    typedef typename %(Dimension)s::ThirdRankTensor ThirdRankTensor;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               nodeList = "const NodeList<%(Dimension)s>&",
               nodeIndices = "const std::vector<int>&"):
        "Construct a constant X velocity for the given nodes"

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(BoundaryAbstractMethods, ConstantVelocityBoundary, virtual=True, pure_virtual=False)
PYB11inject(RestartMethods, ConstantVelocityBoundary)
