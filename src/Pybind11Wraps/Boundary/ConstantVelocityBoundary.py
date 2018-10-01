#-------------------------------------------------------------------------------
# ConstantVelocityBoundary
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Boundary import *
from BoundaryAbstractMethods import *
from RestartMethods import *

@PYB11template("Dimension")
class ConstantVelocityBoundary(Boundary):
    """ConstantVelocityBoundary -- A boundary condition to enforce a constant 
velocity on a given set of nodes.

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
    typedef GeomPlane<DIM> Plane;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               nodeList = "const NodeList<DIM>&",
               nodeIndices = "const std::vector<int>&"):
        "Construct a constant velocity boundary for the specified nodes"

    #...........................................................................
    # Methods
    @PYB11virtual
    @PYB11const
    def valid(self):
        return "bool"

    #...........................................................................
    # Properties
    nodeList = PYB11property("const NodeList<DIM>&", "nodeList", doc="The NodeList this boundary applies to")
    nodeIndices = PYB11property("std::vector<int>", "nodeIndices", doc="The nodes this boundary is in control of")
    velocityCondition = PYB11property("std::vector<Vector>", "velocityCondition", doc="The velocities for the nodes we control")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(BoundaryAbstractMethods, ConstantVelocityBoundary, virtual=True, pure_virtual=False)
PYB11inject(RestartMethods, ConstantVelocityBoundary)
