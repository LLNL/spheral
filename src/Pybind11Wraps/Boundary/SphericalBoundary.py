#-------------------------------------------------------------------------------
# SphericalBoundary
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Boundary import *
from BoundaryAbstractMethods import *
from RestartMethods import *

@PYB11template()
@PYB11template_dict({"Dimension" : "Dim<3>"})
class SphericalBoundary(Boundary):

    PYB11typedefs = """
    typedef Dim<3>::Scalar Scalar;
    typedef Dim<3>::Vector Vector;
    typedef Dim<3>::Tensor Tensor;
    typedef Dim<3>::SymTensor SymTensor;
    typedef Dim<3>::ThirdRankTensor ThirdRankTensor;
    typedef typename %(Dimension)s::FourthRankTensor FourthRankTensor;
    typedef typename %(Dimension)s::FifthRankTensor FifthRankTensor;
    typedef Dim<3>::FacetedVolume FacetedVolume;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               dataBase = "const DataBase<Dim<3>>&"):
        "Construct with the DataBase"

    #...........................................................................
    # Methods
    @PYB11const
    def reflectOperator(self,
                        r0 = "const Vector&",
                        r1 = "const Vector&"):
        return "Tensor"

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, SphericalBoundary)
PYB11inject(BoundaryAbstractMethods, SphericalBoundary, virtual=True, pure_virtual=False)
