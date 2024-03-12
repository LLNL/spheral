from PYB11Generator import *
from ASPHSmoothingScale import *
from SmoothingScaleAbstractMethods import *

#-------------------------------------------------------------------------------
# ASPHSmoothingScalev2
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
class ASPHSmoothingScalev2(ASPHSmoothingScale):

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using Tensor = typename %(Dimension)s::Tensor;
    using SymTensor = typename %(Dimension)s::SymTensor;
    using ScalarField = Field<%(Dimension)s, Scalar>;
    using VectorField = Field<%(Dimension)s, Vector>;
    using TensorField = Field<%(Dimension)s, Tensor>;
    using SymTensorField = Field<%(Dimension)s, SymTensor>;
"""

    def pyinit(self):
        "Default constructor"

    @PYB11const
    @PYB11virtual
    @PYB11pycppname("idealSmoothingScale")
    def idealSmoothingScale_points(self,
                                   H = "const SymTensor&", 
                                   pos = "const Vector&", 
                                   zerothMoment = "const Scalar", 
                                   secondMoment = "const SymTensor&", 
                                   W = "const TableKernel<%(Dimension)s>&", 
                                   hmin = "const typename %(Dimension)s::Scalar", 
                                   hmax = "const typename %(Dimension)s::Scalar", 
                                   hminratio = "const typename %(Dimension)s::Scalar", 
                                   nPerh = "const Scalar", 
                                   connectivityMap = "const ConnectivityMap<%(Dimension)s>&", 
                                   nodeListi = "const unsigned", 
                                   i = "const unsigned"):
        "Determine an 'ideal' H for the given moments."
        return "typename %(Dimension)s::SymTensor"
