#-------------------------------------------------------------------------------
# Generic Kernel bindings.
#-------------------------------------------------------------------------------
from PYB11Generator import *

class SphericalKernelOslo:

    PYB11typedefs = """
    using Scalar = Dim<1>::Scalar;
    using Vector = Dim<1>::Vector;
    using SymTensor = Dim<1>::SymTensor;
"""

    def pyinit00(self,
                 kernel = "const TableKernel<Dim<3>>&",
                 numIntegral = ("const unsigned", "5000u"),
                 numKernel = ("const unsigned", "200u"),
                 useInterpolation = ("const bool", "true")):
        "Construct a tabulated Spherical kernel."

    def pyinit01(self,
                 kernel = "const BSplineKernel<Dim<3>>&",
                 numIntegral = ("const unsigned", "5000u"),
                 numKernel = ("const unsigned", "200u"),
                 useInterpolation = ("const bool", "true")):
        "Construct a tabulated Spherical kernel."

    def pyinit02(self,
                 kernel = "const NBSplineKernel<Dim<3>>&",
                 numIntegral = ("const unsigned", "5000u"),
                 numKernel = ("const unsigned", "200u"),
                 useInterpolation = ("const bool", "true")):
        "Construct a tabulated Spherical kernel."

    def pyinit03(self,
                 kernel = "const W4SplineKernel<Dim<3>>&",
                 numIntegral = ("const unsigned", "5000u"),
                 numKernel = ("const unsigned", "200u"),
                 useInterpolation = ("const bool", "true")):
        "Construct a tabulated Spherical kernel."

    def pyinit04(self,
                 kernel = "const GaussianKernel<Dim<3>>&",
                 numIntegral = ("const unsigned", "5000u"),
                 numKernel = ("const unsigned", "200u"),
                 useInterpolation = ("const bool", "true")):
        "Construct a tabulated Spherical kernel."

    def pyinit05(self,
                 kernel = "const SuperGaussianKernel<Dim<3>>&",
                 numIntegral = ("const unsigned", "5000u"),
                 numKernel = ("const unsigned", "200u"),
                 useInterpolation = ("const bool", "true")):
        "Construct a tabulated Spherical kernel."

    def pyinit06(self,
                 kernel = "const HatKernel<Dim<3>>&",
                 numIntegral = ("const unsigned", "5000u"),
                 numKernel = ("const unsigned", "200u"),
                 useInterpolation = ("const bool", "true")):
        "Construct a tabulated Spherical kernel."

    def pyinit07(self,
                 kernel = "const SincKernel<Dim<3>>&",
                 numIntegral = ("const unsigned", "5000u"),
                 numKernel = ("const unsigned", "200u"),
                 useInterpolation = ("const bool", "true")):
        "Construct a tabulated Spherical kernel."

    def pyinit08(self,
                 kernel = "const NSincPolynomialKernel<Dim<3>>&",
                 numIntegral = ("const unsigned", "5000u"),
                 numKernel = ("const unsigned", "200u"),
                 useInterpolation = ("const bool", "true")):
        "Construct a tabulated Spherical kernel."

    def pyinit09(self,
                 kernel = "const QuarticSplineKernel<Dim<3>>&",
                 numIntegral = ("const unsigned", "5000u"),
                 numKernel = ("const unsigned", "200u"),
                 useInterpolation = ("const bool", "true")):
        "Construct a tabulated Spherical kernel."

    def pyinit10(self,
                 kernel = "const QuinticSplineKernel<Dim<3>>&",
                 numIntegral = ("const unsigned", "5000u"),
                 numKernel = ("const unsigned", "200u"),
                 useInterpolation = ("const bool", "true")):
        "Construct a tabulated Spherical kernel."

    def pyinit11(self,
                 kernel = "const WendlandC2Kernel<Dim<3>>&",
                 numIntegral = ("const unsigned", "5000u"),
                 numKernel = ("const unsigned", "200u"),
                 useInterpolation = ("const bool", "true")):
        "Construct a tabulated Spherical kernel."

    def pyinit12(self,
                 kernel = "const WendlandC4Kernel<Dim<3>>&",
                 numIntegral = ("const unsigned", "5000u"),
                 numKernel = ("const unsigned", "200u"),
                 useInterpolation = ("const bool", "true")):
        "Construct a tabulated Spherical kernel."

    def pyinit13(self,
                 kernel = "const WendlandC6Kernel<Dim<3>>&",
                 numIntegral = ("const unsigned", "5000u"),
                 numKernel = ("const unsigned", "200u"),
                 useInterpolation = ("const bool", "true")):
        "Construct a tabulated Spherical kernel."

    def pyinit_copy(self,
                    rhs = "const SphericalKernelOslo&"):
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
    @PYB11implementation("""[](const SphericalKernelOslo& self, const Vector& etaj, const Vector& etai, const SymTensor& H) -> py::tuple {
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
    Winterpolator = PYB11property(doc="Interpolator for the kernel value")
    baseKernel3d = PYB11property(doc="The base 3D kernel")
    baseKernel1d = PYB11property(doc="The base 1D kernel (mostly for IdealH usage)")
    etamax = PYB11property(doc="The maximum kernel extent of the base 3D kernel")
    useInterpolation = PYB11property("bool", getter="useInterpolation", setter="useInterpolation",
                                     doc="Use interpolation to fit integral or numerically evaluate every time")
