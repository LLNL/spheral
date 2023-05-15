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

@PYB11namespace("Spheral::PythonBoundFunctors")
@PYB11template("argT1", "argT2", "retT")
class Spheral2ArgFunctor:
    def pyinit(self):
        return

    @PYB11pure_virtual
    @PYB11const
    def __call__(self,
                 x = "%(argT1)s",
                 y = "%(argT2)s"):
        "Required operator() to map %(argT1)s %(argT2)s --> %(retT)s"
        return "%(retT)s"
