from PYB11Generator import *

class RKFieldNames:

    @PYB11static
    def rkCorrections(order = "const RKOrder"):
        return "std::string"

    @PYB11static
    def reproducingKernel(order = "const RKOrder"):
        return "std::string"
