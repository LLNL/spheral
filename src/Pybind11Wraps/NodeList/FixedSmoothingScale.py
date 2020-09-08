from PYB11Generator import *
from SmoothingScaleBase import *
from SmoothingScaleAbstractMethods import *

#-------------------------------------------------------------------------------
# FixedSmoothingScale
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
class FixedSmoothingScale(SmoothingScaleBase):

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
"""

    def pyinit(self):
        "Default constructor"

#-------------------------------------------------------------------------------
# Add the abstract interface
#-------------------------------------------------------------------------------
PYB11inject(SmoothingScaleAbstractMethods, FixedSmoothingScale, virtual=True)
