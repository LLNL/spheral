#-------------------------------------------------------------------------------
# Generic Kernel bindings.
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11template("Dimension", "Descendent")
class Kernel:

    typedefs = """
    typedef %(Dimension)s::Vector Vector;
    typedef %(Dimension)s::SymTensor SymTensor;
"""

    def pyinit(self):
        "Default constructor"

    @PYB11const
    def __call__(self,
                 etaMagnitude = "double",
                 H = "const SymTensor&"):
        "Return the kernel value at the given normalized radius"
        return "double"

    @PYB11const
    def kernelValue(self,
                    etaMagnitude = "double",
                    H = "const SymTensor&"):
        "Return the kernel value at the given normalized radius"
        return "double"

    @PYB11const
    def kernelValue(self,
                    eta = "const Vector&",
                    H = "const SymTensor&"):
        "Return the kernel value at the given normalized position"
        return "double"

    @PYB11const
    def volumeNormalization(self):
        "The volume normalization constant."
        return "double"
