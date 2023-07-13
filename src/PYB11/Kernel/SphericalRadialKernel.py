#-------------------------------------------------------------------------------
# SphericalRadialKernel -- for use in spherical coordinate specializations
#-------------------------------------------------------------------------------
from PYB11Generator import *

class SphericalRadialKernel:

    PYB11typedefs = """
    using Scalar = Dim<1>::Scalar;
    using Vector = Dim<1>::Vector;
    using SymTensor = Dim<1>::SymTensor;
"""

    def pyinit00(self,
                 kernel = "const TableKernel<Dim<1>>&",
                 numIntegral = ("const unsigned", "5000u"),
                 numKernel = ("const unsigned", "200u"),
                 useInterpolation = ("const bool", "true")):
        "Construct a tabulated Spherical kernel."

    def pyinit01(self,
                 kernel = "const BSplineKernel<Dim<1>>&",
                 numIntegral = ("const unsigned", "5000u"),
                 numKernel = ("const unsigned", "200u"),
                 useInterpolation = ("const bool", "true")):
        "Construct a tabulated Spherical kernel."

    def pyinit02(self,
                 kernel = "const NBSplineKernel<Dim<1>>&",
                 numIntegral = ("const unsigned", "5000u"),
                 numKernel = ("const unsigned", "200u"),
                 useInterpolation = ("const bool", "true")):
        "Construct a tabulated Spherical kernel."

    def pyinit03(self,
                 kernel = "const W4SplineKernel<Dim<1>>&",
                 numIntegral = ("const unsigned", "5000u"),
                 numKernel = ("const unsigned", "200u"),
                 useInterpolation = ("const bool", "true")):
        "Construct a tabulated Spherical kernel."

    def pyinit04(self,
                 kernel = "const GaussianKernel<Dim<1>>&",
                 numIntegral = ("const unsigned", "5000u"),
                 numKernel = ("const unsigned", "200u"),
                 useInterpolation = ("const bool", "true")):
        "Construct a tabulated Spherical kernel."

    def pyinit05(self,
                 kernel = "const SuperGaussianKernel<Dim<1>>&",
                 numIntegral = ("const unsigned", "5000u"),
                 numKernel = ("const unsigned", "200u"),
                 useInterpolation = ("const bool", "true")):
        "Construct a tabulated Spherical kernel."

    def pyinit06(self,
                 kernel = "const HatKernel<Dim<1>>&",
                 numIntegral = ("const unsigned", "5000u"),
                 numKernel = ("const unsigned", "200u"),
                 useInterpolation = ("const bool", "true")):
        "Construct a tabulated Spherical kernel."

    def pyinit07(self,
                 kernel = "const SincKernel<Dim<1>>&",
                 numIntegral = ("const unsigned", "5000u"),
                 numKernel = ("const unsigned", "200u"),
                 useInterpolation = ("const bool", "true")):
        "Construct a tabulated Spherical kernel."

    def pyinit08(self,
                 kernel = "const NSincPolynomialKernel<Dim<1>>&",
                 numIntegral = ("const unsigned", "5000u"),
                 numKernel = ("const unsigned", "200u"),
                 useInterpolation = ("const bool", "true")):
        "Construct a tabulated Spherical kernel."

    def pyinit09(self,
                 kernel = "const QuarticSplineKernel<Dim<1>>&",
                 numIntegral = ("const unsigned", "5000u"),
                 numKernel = ("const unsigned", "200u"),
                 useInterpolation = ("const bool", "true")):
        "Construct a tabulated Spherical kernel."

    def pyinit10(self,
                 kernel = "const QuinticSplineKernel<Dim<1>>&",
                 numIntegral = ("const unsigned", "5000u"),
                 numKernel = ("const unsigned", "200u"),
                 useInterpolation = ("const bool", "true")):
        "Construct a tabulated Spherical kernel."

    def pyinit11(self,
                 kernel = "const WendlandC2Kernel<Dim<1>>&",
                 numIntegral = ("const unsigned", "5000u"),
                 numKernel = ("const unsigned", "200u"),
                 useInterpolation = ("const bool", "true")):
        "Construct a tabulated Spherical kernel."

    def pyinit12(self,
                 kernel = "const WendlandC4Kernel<Dim<1>>&",
                 numIntegral = ("const unsigned", "5000u"),
                 numKernel = ("const unsigned", "200u"),
                 useInterpolation = ("const bool", "true")):
        "Construct a tabulated Spherical kernel."

    def pyinit13(self,
                 kernel = "const WendlandC6Kernel<Dim<1>>&",
                 numIntegral = ("const unsigned", "5000u"),
                 numKernel = ("const unsigned", "200u"),
                 useInterpolation = ("const bool", "true")):
        "Construct a tabulated Spherical kernel."

    def pyinit_copy(self,
                    rhs = "const SphericalRadialKernel&"):
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
    @PYB11implementation("""[](const SphericalRadialKernel& self, const Vector& etaj, const Vector& etai, const SymTensor& H) -> py::tuple {
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

    @PYB11const
    def volumeNormalization(self,
                            eta = "const Scalar"):
        "Return the normalization 'constant' A as a function of r/h"
        return "Scalar"

    @PYB11const
    def gradAInv(self,
                 eta = "const Scalar"):
        "Return the gradient as a function of r/h of the volume normalization A(r/h)"
        return "Scalar"

    #---------------------------------------------------------------------------
    # Attributes
    Ainterpolator = PYB11property(doc="Interpolator for the kernel normalization value")
    gradAInvInterpolator = PYB11property(doc="Interpolator for the gradient of the kernel normalization")
    baseKernel1d = PYB11property(doc="The base 1D kernel")
    etamax = PYB11property(doc="The maximum kernel extent of the base kernel")
    etacutoff = PYB11property(doc="Interpolations used for eta in [0, etacutoff] -- beyond this point we scale the last values for interpolation")
    useInterpolation = PYB11property("bool", getter="useInterpolation", setter="useInterpolation",
                                     doc="Use interpolation to fit integral or numerically evaluate every time")
