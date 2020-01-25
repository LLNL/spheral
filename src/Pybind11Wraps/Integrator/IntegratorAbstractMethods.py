#-------------------------------------------------------------------------------
# Integrator pure virtual interface
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11ignore
class IntegratorAbstractMethods:

    def step(self,
             maxTime = "Scalar",
             state = "State<%(Dimension)s>&",
             derivs = "StateDerivatives<%(Dimension)s>&"):
        "Master method to take a step, i.e., advance in time"
        return "bool"
