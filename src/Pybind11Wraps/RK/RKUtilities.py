#-------------------------------------------------------------------------------
# RKCorrections
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11template("Dimension", "RKOrder correctionOrder")
class RKUtilities:
    "Computes RK correction terms"
    
    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::SymTensor SymTensor;
"""

    def pyinit(self):
        "Constructor"
        
    @PYB11static
    def getPolynomials(self,
                       x = "const Vector&"):
        "Get polynomial vector"
        return "std::vector<double>"
    @PYB11static
    def getGradPolynomials(self,
                           x = "const Vector&"):
        "Get gradient polynomial vector"
        return "std::vector<double>"
    @PYB11static
    def getHessPolynomials(self, x = "const Vector&"):
        "Get hessian polynomial vector"
        return "std::vector<double>"
    
    @PYB11static
    def evaluateBaseKernel(self,
                           kernel = "const TableKernel<%(Dimension)s>&",
                           x = "const Vector&",
                           H = "const SymTensor&"):
        "Evaluate base kernel"
        return "Scalar"
    @PYB11static
    def evaluateBaseGradient(self,
                             kernel = "const TableKernel<%(Dimension)s>&",
                             x = "const Vector&",
                             H = "const SymTensor&"):
        "Evaluate base gradient"
        return "Vector"
    @PYB11static
    def evaluateBaseHessian(self,
                            kernel = "const TableKernel<%(Dimension)s>&",
                            x = "const Vector&",
                            H = "const SymTensor&"):
        "Evaluate base hessian"
        return "SymTensor"

    @PYB11static
    def evaluateKernel(self,
                       kernel = "const TableKernel<%(Dimension)s>&",
                       x = "const Vector&",
                       H = "const SymTensor&",
                       corrections = "const std::vector<double>&"):
        "Evaluate base kernel"
        return "Scalar"
    @PYB11static
    def evaluateGradient(self,
                         kernel = "const TableKernel<%(Dimension)s>&",
                         x = "const Vector&",
                         H = "const SymTensor&",
                         corrections = "const std::vector<double>&"):
        "Evaluate base gradient"
        return "Vector"
    @PYB11static
    def evaluateHessian(self,
                        kernel = "const TableKernel<%(Dimension)s>&",
                        x = "const Vector&",
                        H = "const SymTensor&",
                        corrections = "const std::vector<double>&"):
        "Evaluate base hessian"
        return "SymTensor"
    
    @PYB11static
    def computeCorrections(self,
                           connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
                           kernel = "const TableKernel<%(Dimension)s>&",
                           volume = "const FieldList<%(Dimension)s, Scalar>&",
                           position = "const FieldList<%(Dimension)s, Vector>&",
                           H = "const FieldList<%(Dimension)s, SymTensor>&",
                           needHessian = "const bool",
                           zerothCorrections = "FieldList<%(Dimension)s, std::vector<double>>&",
                           corrections = "FieldList<%(Dimension)s, std::vector<double>>&"):
        "Compute RK corrections"
        return "void"

    @PYB11static
    def computeNormal(self,
                      connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
                      kernel = "const TableKernel<%(Dimension)s>&",
                      volume = "const FieldList<%(Dimension)s, Scalar>&",
                      position = "const FieldList<%(Dimension)s, Vector>&",
                      H = "const FieldList<%(Dimension)s, SymTensor>&",
                      corrections = "const FieldList<%(Dimension)s, std::vector<double>>&",
                      surfaceArea = "FieldList<%(Dimension)s, Scalar>&",
                      normal = "FieldList<%(Dimension)s, Vector>&"):
        "Compute RK corrections"
        return "void"
    
    
