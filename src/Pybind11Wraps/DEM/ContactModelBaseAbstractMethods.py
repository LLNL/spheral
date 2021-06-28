#-------------------------------------------------------------------------------
# ContactModel pure virtual interface
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11ignore
class ContactModelBaseAbstractMethods:

    @PYB11const
    def evaluateDerivatives(self,
                            time     = "const Scalar",
                            dt       = "const Scalar",
                            dataBase = "const DataBase<%(Dimension)s>&",
                            state    = "const State<%(Dimension)s>&",
                            derivs   = "StateDerivatives<%(Dimension)s>&"):
        "Increment the derivatives."
        return "void"

    @PYB11const
    def timeStep(self,
                 dataBase = "const DataBase<%(Dimension)s>&",
                 state    = "const State<%(Dimension)s>&",
                 derivs   = "const StateDerivatives<%(Dimension)s>&",
                 time     = "Scalar"):
        "All constact models must specify how time-step is calculated."
        return "Scalar"

