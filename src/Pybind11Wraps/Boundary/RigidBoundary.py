#-------------------------------------------------------------------------------
# RigidBoundary
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Boundary import *
from PlanarBoundary import *
from BoundaryAbstractMethods import *
from RestartMethods import *

@PYB11template("Dimension")
class RigidBoundary(PlanarBoundary):

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
    def pyinit(self):
        "Default constructor"

    def pyinit1(self,
                plane = "const Plane&"):
        "Construct with the rigid plane"

    #...........................................................................
    # Properties
    reflectOperator = PYB11property("const Tensor&", "reflectOperator", doc="The tensor reflection transformation")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, RigidBoundary)
PYB11inject(BoundaryAbstractMethods, RigidBoundary, virtual=True, pure_virtual=False)
