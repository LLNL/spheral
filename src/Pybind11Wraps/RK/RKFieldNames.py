from PYB11Generator import *

class RKFieldNames:

    rkOrders = PYB11readonly(static=True)

    @PYB11static
    def rkCorrections(order = "const RKOrder"):
        return "std::string"

    @PYB11static
    def reproducingKernel(order = "const RKOrder"):
        return "std::string"
