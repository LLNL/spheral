from PYB11Generator import *
from SmoothingScaleAbstractMethods import *

#-------------------------------------------------------------------------------
# SmoothingScaleBase
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
class SmoothingScaleBase:

    typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    typedef Field<%(Dimension)s, Scalar> ScalarField;
    typedef Field<%(Dimension)s, Vector> VectorField;
    typedef Field<%(Dimension)s, Tensor> TensorField;
    typedef Field<%(Dimension)s, SymTensor> SymTensorField;
"""

    def pyinit(self):
        "Default constructor"

    @PYB11const
    def newSmoothingScaleAndDerivative(self,
                                       H = "const SymTensorField&",
                                       position = "const VectorField&",
                                       DvDx = "const TensorField&",
                                       zerothMoment = "const ScalarField&", 
                                       secondMoment = "const SymTensorField&", 
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
    

#-------------------------------------------------------------------------------
# Add the abstract interface
#-------------------------------------------------------------------------------
PYB11inject(SmoothingScaleAbstractMethods, SmoothingScaleBase, pure_virtual=True)
