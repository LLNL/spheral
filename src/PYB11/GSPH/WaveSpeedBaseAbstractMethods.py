#-------------------------------------------------------------------------------
# wave Speed pure virtual interface
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11ignore
class WaveSpeedBaseAbstractMethods:

    @PYB11const
    def waveSpeed(self,
                  rhoi = "const Scalar",
                  rhoj = "const Scalar",
                  ci   = "const Scalar",
                  cj   = "const Scalar",
                  ui   = "const Scalar",
                  uj   = "const Scalar",
                  Si   = "Scalar&",
                  Sj   = "Scalar&"):
        "calculate wave speed."
        return "void"
