#-------------------------------------------------------------------------------
# PeriodicBoundary
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Boundary import *
from PlanarBoundary import *
from BoundaryAbstractMethods import *
from RestartMethods import *

@PYB11template("Dimension")
class PeriodicBoundary(PlanarBoundary):

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
    def pyinit(self):
        "Default constructor"

    def pyinit1(self,
                plane1 = "const Plane&",
                plane2 = "const Plane&"):
        "Construct a periodic boundary mapping between the two (enter/exit) planes"

    #...........................................................................
    # Methods
    @PYB11virtual
    def cullGhostNodes(self,
                       flagSet = "const FieldList<DIM, int>&",
                       old2newIndexMap = "FieldList<DIM, int>&",
                       numNodesRemoved = "std::vector<int>&"):
        "Use a set of flags to cull out inactive ghost nodes."
        return "void"
    
    @PYB11virtual
    def reset(self,
              dataBase = "const DataBase<DIM>&"):
        "Overridable hook for clearing out the boundary condition."
        return "void"

    @PYB11virtual
    @PYB11const
    def label(self):
        "Label for restart files"
        return "std::string"

    #...........................................................................
    # Properties
    enterPlane = PYB11property("const Plane&", "enterPlane", "setEnterPlane", doc="The first plane for periodic wrapping")
    exitPlane = PYB11property("const Plane&", "exitPlane", "setExitPlane", doc="The second plane for periodic wrapping")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(BoundaryAbstractMethods, PeriodicBoundary, virtual=True, pure_virtual=False)
