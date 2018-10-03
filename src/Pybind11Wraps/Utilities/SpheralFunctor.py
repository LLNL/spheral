#-------------------------------------------------------------------------------
# The functors for integration and such
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11namespace("Spheral::PythonBoundFunctors")
@PYB11template("argT", "retT")
class SpheralFunctor:
    def pyinit(self):
        return

    @PYB11pure_virtual
    @PYB11const
    def __call__(self, x="%(argT)s"):
        "Required operator() to map %(argT)s --> %(retT)s"
        return "%(retT)s"
