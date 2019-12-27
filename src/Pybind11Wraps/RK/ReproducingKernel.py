#-------------------------------------------------------------------------------
# ReproducingKernel
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11template("Dimension")
class ReproducingKernel:
    """Provides the reproducing kernel methods, analogous to the Kernel class for SPH

This is really just a convenient front-end for the methods in RKUtilities"""
    
    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::SymTensor SymTensor;
"""

    def pyinit(self,
               W = "const TableKernel<%(Dimension)s>&",
               order = "const RKOrder"):
        "Constructor"
        
    @PYB11const
    def evaluateBaseKernel(self,
                           x = "const Vector&",
                           H = "const SymTensor&"):
        "Evaluate base kernel"
        return "Scalar"

    @PYB11const
    def evaluateBaseGradient(self,
                             x = "const Vector&",
                             H = "const SymTensor&"):
        "Evaluate base gradient"
        return "Vector"

    @PYB11const
    def evaluateBaseHessian(self,
                            x = "const Vector&",
                            H = "const SymTensor&"):
        "Evaluate base hessian"
        return "SymTensor"


    @PYB11const
    def evaluateBaseKernelAndGradient(self,
                                      x = "const Vector&",
                                      H = "const SymTensor&"):
        "Return (W(x,h), gradW(x,h)) for the base (uncorrected) kernel"
        return "std::pair<%(Dimension)s::Scalar, %(Dimension)s::Vector>"

    @PYB11const
    def evaluateKernel(self,
                       x = "const Vector&",
                       H = "const SymTensor&",
                       corrections = "const std::vector<double>&"):
        "Evaluate base kernel"
        return "Scalar"

    @PYB11const
    def evaluateGradient(self,
                         x = "const Vector&",
                         H = "const SymTensor&",
                         corrections = "const std::vector<double>&"):
        "Evaluate base gradient"
        return "Vector"

    @PYB11const
    def evaluateHessian(self,
                        x = "const Vector&",
                        H = "const SymTensor&",
                        corrections = "const std::vector<double>&"):
        "Evaluate base hessian"
        return "SymTensor"
    
    @PYB11const
    def evaluateKernelAndGradient(self,
                                  x = "const Vector&",
                                  H = "const SymTensor&"):
        "Return (WR(x,h), gradWR(x,h)) for the reproducing kernel"
        return "std::pair<%(Dimension)s::Scalar, %(Dimension)s::Vector>"

    @PYB11const
    def evaluateKernelAndGradients(self,
                                   x = "const Vector&",
                                   H = "const SymTensor&"):
        "Return (WR(x,h), gradWR(x,h), ||gradW(x,h)||) for the reproducing kernel and base kernel"
        return "std::tuple<%(Dimension)s::Scalar, %(Dimension)s::Vector, %(Dimension)s::Scalar>"

    @PYB11const
    def computeCorrections(self,
                           connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
                           volume = "const FieldList<%(Dimension)s, Scalar>&",
                           position = "const FieldList<%(Dimension)s, Vector>&",
                           H = "const FieldList<%(Dimension)s, SymTensor>&",
                           needHessian = "const bool",
                           zerothCorrections = "FieldList<%(Dimension)s, std::vector<double>>&",
                           corrections = "FieldList<%(Dimension)s, std::vector<double>>&"):
        "Compute RK corrections"
        return "void"

    @PYB11const
    def computeNormal(self,
                      connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
                      volume = "const FieldList<%(Dimension)s, Scalar>&",
                      position = "const FieldList<%(Dimension)s, Vector>&",
                      H = "const FieldList<%(Dimension)s, SymTensor>&",
                      corrections = "const FieldList<%(Dimension)s, std::vector<double>>&",
                      surfaceArea = "FieldList<%(Dimension)s, Scalar>&",
                      normal = "FieldList<%(Dimension)s, Vector>&"):
        "Compute RK corrections"
        return "void"
    
    #..........................................................................
    # Attributes
    order = PYB11property(doc="order to which we are enforcing reproducibility")
    kernel = PYB11property(doc="The base (uncorrected) interpolation kernel")
