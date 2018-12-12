#-------------------------------------------------------------------------------
# ArtificialViscosity pure virtual interface
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11ignore
class ArtificialViscosityAbstractMethods:

    @PYB11const
    def Piij(self,
             nodeListi = "const unsigned",
             i = "const unsigned",
             nodeListj = "const unsigned",
             j = "const unsigned",
             xi = "const Vector&",
             etai = "const Vector&",
             vi = "const Vector&",
             rhoi = "const Scalar",
             csi = "const Scalar",
             Hi = "const SymTensor&",
             xj = "const Vector&",
             etaj = "const Vector&",
             vj = "const Vector&",
             rhoj = "const Scalar",
             csj = "const Scalar",
             Hj = "const SymTensor&"):
        "Require all descendents to return the artificial viscous Pi = P/rho^2 as a tensor. Scalar viscosities should just return a diagonal tensor with their value along the diagonal."
        return "std::pair<Tensor, Tensor>"
