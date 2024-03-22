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
"""

    def pyinit(self):
        "Default constructor"

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
