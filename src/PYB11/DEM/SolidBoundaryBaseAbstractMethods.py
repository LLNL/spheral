#-------------------------------------------------------------------------------
# pure virtual for solid boundary class lineage
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11ignore
class SolidBoundaryBaseAbstractMethods:

    @PYB11const
    def velocity(self,
                 position = "const Vector&"):
        "velocity of bc."
        return "Vector"

    @PYB11const
    def distance(self,
                 position = "const Vector&"):
        "distance vector to bc."
        return "Vector"

    def update(self,
               multiplier = "const double",
               t = "const double",
               dt = "const double",):
        "distance vector to bc."
        return "void"
