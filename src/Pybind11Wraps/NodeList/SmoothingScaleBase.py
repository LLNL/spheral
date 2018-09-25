from PYB11Generator import *

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
    
    @PYB11pure_virtual
    @PYB11const
    def smoothingScaleDerivative(self,
                                 H = "const SymTensor&",
                                 pos = "const Vector&", 
                                 DvDx = "const Tensor&", 
                                 hmin = "const Scalar", 
                                 hmax = "const Scalar", 
                                 hminratio = "const Scalar", 
                                 nPerh = "const Scalar"):
        "Time derivative of the smoothing scale."
        return "SymTensor"
  
    @PYB11pure_virtual
    @PYB11const
    def newSmoothingScale(self,
                          H = "const SymTensor&", 
                          pos = "const Vector&", 
                          zerothMoment = "const Scalar", 
                          secondMoment = "const SymTensor&", 
                          W = "const TableKernel<%(Dimension)s>&", 
                          hmin = "const Scalar", 
                          hmax = "const Scalar", 
                          hminratio = "const Scalar", 
                          nPerh = "const Scalar", 
                          connectivityMap = "const ConnectivityMap<%(Dimension)s>&", 
                          nodeListi = "const unsigned", 
                          i = "const unsigned"):
        "Return a new H, with limiting based on the old value."
        return "SymTensor"

    @PYB11pure_virtual
    @PYB11const
    def idealSmoothingScale(self,
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

    @PYB11pure_virtual
    @PYB11const
    @PYB11pycppname("idealSmoothingScale")
    def idealSmoothingScale1(self,
                             H = "const SymTensor&", 
                             mesh = "const Mesh<%(Dimension)s>&", 
                             zone = "const typename Mesh<%(Dimension)s>::Zone&", 
                             hmin = "const Scalar", 
                             hmax = "const Scalar", 
                             hminratio = "const Scalar", 
                             nPerh = "const Scalar"):
        "Compute the new H tensors for a tessellation."
        return "SymTensor"
