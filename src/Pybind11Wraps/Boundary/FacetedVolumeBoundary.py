#-------------------------------------------------------------------------------
# ReflectingBoundary
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Boundary import *
from PlanarBoundary import *
from BoundaryAbstractMethods import *
from RestartMethods import *

@PYB11template("Dimension")
class FacetedVolumeBoundary(Boundary):

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
               poly = "const FacetedVolume&",
               interiorBoundary = "const bool",
               useGhosts = ("const bool", "false")):
        "Constructor"

    #...........................................................................
    # Methods
    @PYB11virtual
    def reset(self,
              dataBase = "const DataBase<%(Dimension)s>&"):
        "Overridable hook for clearing out the boundary condition."
        return "void"

    @PYB11virtual
    def cullGhostNodes(self,
                       flagSet = "const FieldList<%(Dimension)s, int>&",
                       old2newIndexMap = "FieldList<%(Dimension)s, int>&",
                       numNodesRemoved = "std::vector<int>&"):
        "Use a set of flags to cull out inactive ghost nodes."
        return "void"

    @PYB11const
    def reflectOperator(self, 
                        facetID = "unsigned"):
        "Return the effective reflection operator for the given facet"
        return "const Tensor&"

    #...........................................................................
    # Properties
    polyVolume = PYB11property(doc="The faceted volume defining this boundary")
    interiorBoundary = PYB11property(doc="Flag if this Boundary is interior or exterior to the problem")
    useGhosts = PYB11property(doc="Flag whether we should build ghost points or not")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(BoundaryAbstractMethods, FacetedVolumeBoundary, virtual=True, pure_virtual=False)
