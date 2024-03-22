from PYB11Generator import *
from SmoothingScaleBase import *
from SmoothingScaleAbstractMethods import *

#-------------------------------------------------------------------------------
# ASPHSmoothingScale
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
class ASPHSmoothingScale(SmoothingScaleBase):

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using Tensor = typename %(Dimension)s::Tensor;
    using SymTensor = typename %(Dimension)s::SymTensor;
"""

    def pyinit(self):
        "Constructor: setting numPoints == 0 implies create lookup tables with same number of points as TableKernel W"

#-------------------------------------------------------------------------------
# Add the abstract interface
#-------------------------------------------------------------------------------
PYB11inject(SmoothingScaleAbstractMethods, ASPHSmoothingScale, virtual=True)
