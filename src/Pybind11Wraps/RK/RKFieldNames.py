from PYB11Generator import *

class RKFieldNames:

    rkOrders = PYB11readonly(static=True, returnpolicy="copy")
    rkCorrectionsBase = PYB11readonly(static=True, returnpolicy="copy")
    reproducingKernelBase = PYB11readonly(static=True, returnpolicy="copy")

    @PYB11static
    def rkCorrections(order = "const RKOrder"):
        return "std::string"

    @PYB11static
    def reproducingKernel(order = "const RKOrder"):
        return "std::string"

    @PYB11static
    def correctionOrder(x = "const std::string&"):
        return "RKOrder"
