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
                            dataBase = "const DataBase<%(Dimension)s>&",
                            state = "const State<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Increment the derivatives."
        return "void"

    # @PYB11implementation("""[](const Physics<%(Dimension)s>& self,
    #                            const DataBase<%(Dimension)s>& dataBase,
    #                            const State<%(Dimension)s>& state,
    #                            const StateDerivatives<%(Dimension)s>& derivs,
    #                            const Scalar currentTime) { auto result = self.dt(dataBase, state, derivs, currentTime);
    #                                                        return py::make_tuple(result.first, result.second);
    #                            }""")
    @PYB11const
    def dt(dataBase = "const DataBase<%(Dimension)s>&", 
           state = "const State<%(Dimension)s>&",
           derivs = "const StateDerivatives<%(Dimension)s>&",
           currentTime = "const Scalar"):
        "Vote on a time step."
        return "TimeStepType"

    def registerState(self,
                      dataBase = "DataBase<%(Dimension)s>&",
                      state = "State<%(Dimension)s>&"):
        "Register the state you want carried around (and potentially evolved), as well as the policies for such evolution."
        return "void"

    def registerDerivatives(self,
                            dataBase = "DataBase<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Register the derivatives/change fields for updating state."
        return "void"

    @PYB11const
    def label(self):
        "It's useful to have labels for Physics packages.  We'll require this to have the same signature as the restart label."
        return "std::string"
