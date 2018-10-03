#-------------------------------------------------------------------------------
# Integrator pure virtual interface
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11ignore
class IntegratorAbstractMethods:

    def step(self,
             maxTime = "Scalar",
             state = "State<DIM>&",
             derivs = "StateDerivatives<DIM>&"):
        "Master method to take a step, i.e., advance in time"
        return "void"
