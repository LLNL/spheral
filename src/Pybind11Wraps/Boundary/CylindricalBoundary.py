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

    #............................................................................
    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary3(self,
                            field = "Field<%(Dimension)s, Vector>&"):
        "Apply the boundary condition to the ghost node values in the given Field."
        return "void"

    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary4(self,
                            field = "Field<%(Dimension)s, Tensor>&"):
        "Apply the boundary condition to the ghost node values in the given Field."
        return "void"

    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary5(self,
                            field = "Field<%(Dimension)s, SymTensor>&"):
        "Apply the boundary condition to the ghost node values in the given Field."
        return "void"

    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary6(self,
                            field = "Field<%(Dimension)s, ThirdRankTensor>&"):
        "Apply the boundary condition to the ghost node values in the given Field."
        return "void"

    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary7(self,
                            field = "Field<%(Dimension)s, FourthRankTensor>&"):
        "Apply the boundary condition to the ghost node values in the given Field."
        return "void"

    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary8(self,
                            field = "Field<%(Dimension)s, FifthRankTensor>&"):
        "Apply the boundary condition to the ghost node values in the given Field."
        return "void"
    
    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary9(self,
                            field = "Field<%(Dimension)s, FacetedVolume>&"):
        "Apply the boundary condition to the ghost node values in the given Field."
        return "void"

    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary10(self,
                             field = "Field<%(Dimension)s, RKCoefficients<%(Dimension)s>>&"):
        "Apply the boundary condition to the ghost node values in the given Field."
        return "void"

    #............................................................................
    @PYB11pycppname("enforceBoundary")
    @PYB11virtual
    @PYB11const
    def enforceBoundary3(self,
                         field = "Field<%(Dimension)s, Vector>&"):
        "Apply the boundary condition to the violation node values in the given Field."
        return "void"

    @PYB11pycppname("enforceBoundary")
    @PYB11virtual
    @PYB11const
    def enforceBoundary4(self,
                         field = "Field<%(Dimension)s, Tensor>&"):
        "Apply the boundary condition to the violation node values in the given Field."
        return "void"

    @PYB11pycppname("enforceBoundary")
    @PYB11virtual
    @PYB11const
    def enforceBoundary5(self,
                         field = "Field<%(Dimension)s, SymTensor>&"):
        "Apply the boundary condition to the violation node values in the given Field."
        return "void"

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
