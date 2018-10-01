#-------------------------------------------------------------------------------
# ConstantBoundary
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Boundary import *
from BoundaryAbstractMethods import *
from RestartMethods import *

@PYB11template("Dimension")
class ConstantBoundary(Boundary):
    """ConstantBoundary -- Take a snapshot of the state of a set of nodes on 
NodeList, and create ghost nodes with that state from that time on.

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
               nodeList = "NodeList<DIM>&",
               nodeIndices = "const std::vector<int>&",
               denialPlane = "const Plane&"):
        "Construct a constant boundary for the specified nodes, including plane nodes are not allowed through"

    #...........................................................................
    # Methods
    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary20(self,
                             field = "Field<DIM, std::vector<Scalar>>&"):
        return "void"

    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary21(self,
                             field = "Field<DIM, std::vector<Vector>>&"):
        return "void"

    @PYB11virtual
    def initializeProblemStartup(self):
        return "void"

    @PYB11virtual
    @PYB11const
    def valid(self):
        return "bool"

    #...........................................................................
    # Properties
    nodeList = PYB11property("const NodeList<DIM>&", "nodeList", doc="The NodeList this boundary applies to")
    nodeIndices = PYB11property("std::vector<int>", "nodeIndices", doc="The nodes this boundary is in control of")
    numConstantNodes = PYB11property("int", "numConstantNodes", doc="The number of nodes we are controlling")
    reflectOperator = PYB11property("const Tensor&", "reflectOperator", doc="The tensor reflection transformation")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(BoundaryAbstractMethods, ConstantBoundary, virtual=True, pure_virtual=False)
PYB11inject(RestartMethods, ConstantBoundary)
