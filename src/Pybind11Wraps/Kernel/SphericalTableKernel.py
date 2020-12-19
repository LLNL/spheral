#-------------------------------------------------------------------------------
# Generic Kernel bindings.
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Kernel import Kernel

@PYB11template()
@PYB11template_dict({"Dimension" : "Dim<1>",
                     "Descendant" : "SphericalTableKernel"})
class SphericalTableKernel(Kernel):

    PYB11typedefs = """
    typedef typename Dim<1>::Vector Vector;
"""

    def pyinit(self,
               kernel = "const TableKernel<Dim<3>>&"):
        "Construct a tabulated Spherical kernel based on an orindary 3D TableKernel."

    def pyinit_copy(self,
                    rhs = "const SphericalTableKernel&"):
        "Copy constructor"

    @PYB11const
    def __call__(self,
                 etaj = "const Vector&",
                 etai = "const Vector&"):
        "Return the kernel value at the given (rj/h, ri/h) == (etaj, etai) pair"
        return "double"

    #---------------------------------------------------------------------------
    # Attributes
    kernel = PYB11property(doc="The base 3D kernel")
    retamax = PYB11property(doc="The maximum interpolation r/h value")
