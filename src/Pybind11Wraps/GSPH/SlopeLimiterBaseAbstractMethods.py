#-------------------------------------------------------------------------------
# wave Speed pure virtual interface
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11ignore
class SlopeLimiterBaseAbstractMethods:

    @PYB11const
    def fluxLimiter(self,
                    x = "const Scalar"):
        "flux limiter."
        return "Scalar"
