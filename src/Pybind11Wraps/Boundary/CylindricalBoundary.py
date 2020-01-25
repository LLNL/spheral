#-------------------------------------------------------------------------------
# CylindricalBoundary
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Boundary import *
from BoundaryAbstractMethods import *
from RestartMethods import *

@PYB11template()
@PYB11template_dict({"Dimension" : "Dim<3>"})
class CylindricalBoundary(Boundary):

    PYB11typedefs = """
    typedef Dim<3>::Scalar Scalar;
    typedef Dim<3>::Vector Vector;
    typedef Dim<3>::Tensor Tensor;
    typedef Dim<3>::SymTensor SymTensor;
    typedef Dim<3>::ThirdRankTensor ThirdRankTensor;
    typedef Dim<3>::FourthRankTensor FourthRankTensor;
    typedef Dim<3>::FifthRankTensor FifthRankTensor;
    typedef Dim<3>::FacetedVolume FacetedVolume;
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

    @PYB11static
    def angularSpacing(self,
                       ri = "const double",
                       hzi = "const double",
                       nodePerh = "const double",
                       kernelExtent = "const double"):
        return "double"

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, CylindricalBoundary)
PYB11inject(BoundaryAbstractMethods, CylindricalBoundary, virtual=True, pure_virtual=False)
