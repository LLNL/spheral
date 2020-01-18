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

    #............................................................................
    # Constructors
    def pyinit(self,
               dataBase = "const DataBase<Dim<3>>&"):
        "Construct with the DataBase"

    #............................................................................
    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary1(self,
                            field = "Field<%(Dimension)s, int>&"):
        "Apply the boundary condition to the ghost node values in the given Field."
        return "void"

    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary2(self,
                            field = "Field<%(Dimension)s, Scalar>&"):
        "Apply the boundary condition to the ghost node values in the given Field."
        return "void"

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

    #............................................................................
    @PYB11pycppname("enforceBoundary")
    @PYB11virtual
    @PYB11const
    def enforceBoundary1(self,
                         field = "Field<%(Dimension)s, int>&"):
        "Apply the boundary condition to the violation node values in the given Field."
        return "void"

    @PYB11pycppname("enforceBoundary")
    @PYB11virtual
    @PYB11const
    def enforceBoundary2(self,
                         field = "Field<%(Dimension)s, Scalar>&"):
        "Apply the boundary condition to the violation node values in the given Field."
        return "void"

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

    @PYB11pycppname("enforceBoundary")
    @PYB11virtual
    @PYB11const
    def enforceBoundary6(self,
                         field = "Field<%(Dimension)s, ThirdRankTensor>&"):
        "Apply the boundary condition to the violation node values in the given Field."
        return "void"

    @PYB11pycppname("enforceBoundary")
    @PYB11virtual
    @PYB11const
    def enforceBoundary7(self,
                         field = "Field<%(Dimension)s, FourthRankTensor>&"):
        "Apply the boundary condition to the violation node values in the given Field."
        return "void"

    @PYB11pycppname("enforceBoundary")
    @PYB11virtual
    @PYB11const
    def enforceBoundary8(self,
                         field = "Field<%(Dimension)s, FifthRankTensor>&"):
        "Apply the boundary condition to the violation node values in the given Field."
        return "void"
    
    @PYB11pycppname("enforceBoundary")
    @PYB11virtual
    @PYB11const
    def enforceBoundary9(self,
                         field = "Field<%(Dimension)s, FacetedVolume>&"):
        "Apply the boundary condition to the violation node values in the given Field."
        return "void"

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
