from PYB11Generator import *

#-------------------------------------------------------------------------------
# Helper for (re)generating abstract SmoothingScaleBase interface
#-------------------------------------------------------------------------------
@PYB11ignore
class SmoothingScaleAbstractMethods:

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
  
    @PYB11const
    def newSmoothingScale(self,
                          H = "const SymTensor&", 
                          pos = "const FieldList<%(Dimension)s, Vector>&", 
                          zerothMoment = "const Scalar", 
                          firstMoment = "const Vector&", 
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

    @PYB11const
    def idealSmoothingScale(self,
                            H = "const SymTensor&", 
                            pos = "const FieldList<%(Dimension)s, Vector>&", 
                            zerothMoment = "const Scalar", 
                            firstMoment = "const Vector&", 
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
