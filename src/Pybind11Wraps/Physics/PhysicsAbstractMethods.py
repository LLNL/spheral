#-------------------------------------------------------------------------------
# Physics pure virtual interface
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11ignore
class PhysicsAbstractMethods:

    @PYB11const
    def evaluateDerivatives(self,
                            time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<DIM>&",
                            state = "const State<DIM>&",
                            derivs = "StateDerivatives<DIM>&"):
        "Increment the derivatives."
        return "void"

    @PYB11const
    def dt(dataBase = "const DataBase<DIM>&", 
           state = "const State<DIM>&",
           derivs = "const StateDerivatives<DIM>&",
           currentTime = "const Scalar"):
        "Vote on a time step."
        return "TimeStepType"

    def registerState(self,
                      dataBase = "DataBase<DIM>&",
                      state = "State<DIM>&"):
        "Register the state you want carried around (and potentially evolved), as well as the policies for such evolution."
        return "void"

    def registerDerivatives(self,
                            dataBase = "DataBase<DIM>&",
                            derivs = "StateDerivatives<DIM>&"):
        "Register the derivatives/change fields for updating state."
        return "void"

    @PYB11const
    def label(self):
        "It's useful to have labels for Physics packages.  We'll require this to have the same signature as the restart label."
        return "std::string"
