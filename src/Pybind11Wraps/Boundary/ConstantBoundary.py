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

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    typedef typename %(Dimension)s::ThirdRankTensor ThirdRankTensor;
    typedef typename %(Dimension)s::FourthRankTensor FourthRankTensor;
    typedef typename %(Dimension)s::FifthRankTensor FifthRankTensor;
    typedef typename %(Dimension)s::FacetedVolume FacetedVolume;
    typedef GeomPlane<%(Dimension)s> Plane;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               nodeList = "NodeList<%(Dimension)s>&",
               nodeIndices = "const std::vector<int>&",
               denialPlane = "const Plane&"):
        "Construct a constant boundary for the specified nodes, including plane nodes are not allowed through"

    #...........................................................................
    # Methods
    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary20(self,
                             field = "Field<%(Dimension)s, std::vector<Scalar>>&"):
        return "void"

    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary21(self,
                             field = "Field<%(Dimension)s, std::vector<Vector>>&"):
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
    nodeList = PYB11property("const NodeList<%(Dimension)s>&", "nodeList", doc="The NodeList this boundary applies to")
    nodeIndices = PYB11property("std::vector<int>", "nodeIndices", doc="The nodes this boundary is in control of")
    numConstantNodes = PYB11property("int", "numConstantNodes", doc="The number of nodes we are controlling")
    reflectOperator = PYB11property("const Tensor&", "reflectOperator", doc="The tensor reflection transformation")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(BoundaryAbstractMethods, ConstantBoundary, virtual=True, pure_virtual=False)
PYB11inject(RestartMethods, ConstantBoundary)
