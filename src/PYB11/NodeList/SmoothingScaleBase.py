from PYB11Generator import *
from SmoothingScaleAbstractMethods import *

#-------------------------------------------------------------------------------
# SmoothingScaleBase
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralNodeList")
class SmoothingScaleBase:

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
    def newSmoothingScaleAndDerivative(self,
                                       H = "const SymTensorField&",
                                       position = "const VectorField&",
                                       DvDx = "const TensorField&",
                                       zerothMoment = "const ScalarField&", 
                                       firstMoment = "const VectorField&", 
                                       secondMomentEta = "const SymTensorField&", 
                                       secondMomentLab = "const SymTensorField&", 
                                       connectivityMap = "const ConnectivityMap<%(Dimension)s>&", 
                                       W = "const TableKernel<%(Dimension)s>&", 
                                       hmin = "const Scalar", 
                                       hmax = "const Scalar", 
                                       hminratio = "const Scalar", 
                                       nPerh = "const Scalar", 
                                       DHDt = "SymTensorField&", 
                                       Hideal = "SymTensorField&"):
        "Compute the time derivative and ideal H simultaneously for a Field of H's."
        return "void"
    
    @PYB11const
    def hmax(self,
             Vi = "const Scalar",
             nPerh = "const Scalar"):
        "Compute an effective hmax given the volume and target nperh"
        return "Scalar"

#-------------------------------------------------------------------------------
# Add the abstract interface
#-------------------------------------------------------------------------------
PYB11inject(SmoothingScaleAbstractMethods, SmoothingScaleBase, pure_virtual=True)
