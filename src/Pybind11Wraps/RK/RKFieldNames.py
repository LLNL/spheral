from PYB11Generator import *

class RKFieldNames:

    rkOrders = PYB11readonly(static=True)
    rkCorrectionsBase = PYB11readonly(static=True)
    reproducingKernelBase = PYB11readonly(static=True)

    @PYB11static
    def rkCorrections(order = "const RKOrder"):
        return "std::string"

    @PYB11static
    def reproducingKernel(order = "const RKOrder"):
        return "std::string"

    @PYB11static
    def correctionOrder(x = "const std::string&"):
        return "RKOrder"
