#-------------------------------------------------------------------------------
# WeightingFunctor
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11template("Dimension")
class WeightingFunctor:

    PYB11typedefs = """
typedef typename %(Dimension)s::Vector Vector;
typedef typename %(Dimension)s::FacetedVolume FacetedVolume;
"""

    def pyinit(self):
        "Default constructor"

    @PYB11virtual
    @PYB11const
    def __call__(self,
                 pos = "const Vector&",
                 boundary = "const FacetedVolume&"):
        "Override functor __call__ to provide weight"
        return "double"

                 
