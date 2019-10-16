#-------------------------------------------------------------------------------
# InflowBoundary
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Boundary import *
from Physics import *
from BoundaryAbstractMethods import *
from PhysicsAbstractMethods import *
from RestartMethods import *

@PYB11template("Dimension")
class InflowBoundary(Boundary, Physics):
    """InflowBoundary -- creates inflow ghost images, which become internal nodes
as they cross the specified boundary plane."""

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    typedef typename %(Dimension)s::ThirdRankTensor ThirdRankTensor;
    typedef typename %(Dimension)s::FacetedVolume FacetedVolume;
    typedef typename Physics<%(Dimension)s>::TimeStepType TimeStepType;
    typedef GeomPlane<%(Dimension)s> Plane;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               nodeList = "NodeList<%(Dimension)s>&",
               plane = "const GeomPlane<%(Dimension)s>&"):
        "Constructor"

    #...........................................................................
    # Methods
    @PYB11virtual
    def initializeProblemStartup(self):
        "After physics have been initialized we take a snapshot of the node state."
        return "void"

    @PYB11virtual
    def finalize(self,
                 time = "const Scalar",
                 dt = "const Scalar",
                 dataBase = "DataBase<%(Dimension)s>&",
                 state = "State<%(Dimension)s>&",
                 derivs = "StateDerivatives<%(Dimension)s>&"):
        """Packages might want a hook to do some post-step finalizations.
Really we should rename this post-step finalize."""
        return "void"

    #...........................................................................
    # Properties
    numInflowNodes = PYB11property(doc="Number of nodes in inflow stencil")
    nodeList = PYB11property(doc="The NodeList we're allowing to flow into the problem")
    plane = PYB11property(doc="The inflow plane")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, InflowBoundary)
PYB11inject(BoundaryAbstractMethods, InflowBoundary, virtual=True, pure_virtual=False)
PYB11inject(PhysicsAbstractMethods, InflowBoundary, virtual=True, pure_virtual=False)
