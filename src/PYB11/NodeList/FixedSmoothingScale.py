from PYB11Generator import *
from SmoothingScaleBase import *
from SmoothingScaleAbstractMethods import *

#-------------------------------------------------------------------------------
# FixedSmoothingScale
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
class FixedSmoothingScale(SmoothingScaleBase):

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using Tensor = typename %(Dimension)s::Tensor;
    using SymTensor = typename %(Dimension)s::SymTensor;
"""

    def pyinit(self):
        "Default constructor"

#-------------------------------------------------------------------------------
# Add the abstract interface
#-------------------------------------------------------------------------------
PYB11inject(SmoothingScaleAbstractMethods, FixedSmoothingScale, virtual=True)
