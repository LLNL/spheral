from PYB11Generator import *
from SmoothingScaleBase import *
from SmoothingScaleAbstractMethods import *

#-------------------------------------------------------------------------------
# SPHSmoothingScale
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
class SPHSmoothingScale(SmoothingScaleBase):

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
PYB11inject(SmoothingScaleAbstractMethods, SPHSmoothingScale, virtual=True)
