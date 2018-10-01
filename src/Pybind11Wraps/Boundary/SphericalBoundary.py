#-------------------------------------------------------------------------------
# SphericalBoundary
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Boundary import *
from BoundaryAbstractMethods import *
from RestartMethods import *

@PYB11template()
class SphericalBoundary(Boundary):

    typedefs = """
    typedef Dim<3> DIM;
    typedef Dim<3>::Scalar Scalar;
    typedef Dim<3>::Vector Vector;
    typedef Dim<3>::Tensor Tensor;
    typedef Dim<3>::SymTensor SymTensor;
    typedef Dim<3>::ThirdRankTensor ThirdRankTensor;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               dataBase = "const DataBase<Dim<3>>&"):
        "Construct with the DataBase"

    #...........................................................................
    # Methods
    @PYB11static
    def reflectOperator(self,
                        r0 = "const Vector&",
                        r1 = "const Vector&"):
        return "Tensor"

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, SphericalBoundary)
PYB11inject(BoundaryAbstractMethods, SphericalBoundary, virtual=True, pure_virtual=False)
