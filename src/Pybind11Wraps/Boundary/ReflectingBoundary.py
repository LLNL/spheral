#-------------------------------------------------------------------------------
# ReflectingBoundary
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Boundary import *
from PlanarBoundary import *
from BoundaryAbstractMethods import *
from RestartMethods import *

@PYB11template("Dimension")
class ReflectingBoundary(PlanarBoundary):

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
                plane = "const Plane&"):
        "Construct with the reflecting plane"

    #...........................................................................
    # Properties
    reflectOperator = PYB11property("const Tensor&", "reflectOperator", doc="The tensor reflection transformation")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, ReflectingBoundary)
PYB11inject(BoundaryAbstractMethods, ReflectingBoundary, virtual=True, pure_virtual=False)
