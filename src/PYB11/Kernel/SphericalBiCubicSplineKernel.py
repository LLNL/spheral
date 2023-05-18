#-------------------------------------------------------------------------------
# Generic Kernel bindings.
#-------------------------------------------------------------------------------
from PYB11Generator import *

class SphericalBiCubicSplineKernel:

    PYB11typedefs = """
    using Scalar = Dim<1>::Scalar;
    using Vector = Dim<1>::Vector;
    using SymTensor = Dim<1>::SymTensor;
"""

    def pyinit00(self,
                 numKernel = ("const unsigned", "200u")):
        "Construct a SphericalBiCubicSpline kernel."

    def pyinit_copy(self,
                    rhs = "const SphericalBiCubicSplineKernel&"):
        "Copy constructor"

    @PYB11const
    def __call__(self,
                 etaj = "const Vector&",
                 etai = "const Vector&",
                 Hdeti = "const Scalar"):
        "Return the kernel value at the given (rj/h, ri/h) == (etaj, etai) pair"
        return "double"

    @PYB11const
    def grad(self,
             etaj = "const Vector&",
             etai = "const Vector&",
             H = "const SymTensor&"):
        "Return the kernel gradient value at the given (rj/h, ri/h) == (etaj, etai) pair"
        return "Vector"

    @PYB11const
    @PYB11implementation("""[](const SphericalBiCubicSplineKernel& self, const Vector& etaj, const Vector& etai, const SymTensor& H) -> py::tuple {
        double W, deltaWsum;
        Vector gradW;
        self.kernelAndGrad(etaj, etai, H, W, gradW, deltaWsum);
        return py::make_tuple(W, gradW, deltaWsum);
      }""")
    def kernelAndGrad(self,
                      etaj = "const Vector&",
                      etai = "const Vector&",
                      H = "const SymTensor&"):
        "Simultaneously compute the W(etaj, etai), gradW(etaj, etai), deltaWsum(etaj, etai) -- returns a tuple of those results."
        return "py::tuple"

    #---------------------------------------------------------------------------
    # Attributes
    baseKernel3d = PYB11property(doc="The base 3D kernel")
    baseKernel1d = PYB11property(doc="The base 1D kernel (mostly for IdealH usage)")
    etamax = PYB11property(doc="The maximum kernel extent of the base 3D kernel")
