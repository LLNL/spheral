#-------------------------------------------------------------------------------
# QuadraticInterpolator
#-------------------------------------------------------------------------------
from PYB11Generator import *

class QuadraticInterpolator:

    """Encapsulates the algorithm and data for parabolic interpolation in 1D
Assumes the results is interpolated as y_interp = a + b*x + c*x^2"""

    def pyinit(self):
        "Default constuctor -- returns a non-functional interpolator until initialized"
        return

    def pyinit1(self,
                xmin = "const double",
                xmax = "const double",
                yvals = "const std::vector<double>&"):
        "Returns an interpolator for yvals sampled in x in [xmin, xmax]"
        return

    def initialize(self,
                   xmin = "const double",
                   xmax = "const double",
                   yvals = "const std::vector<double>&"):
        "Initializes the interpolator for yvals sampled in x in [xmin, xmax]"

    @PYB11const
    def __call__(self,
                 x = "const double"):
        "Returns the interpolated value <y>(x)"
        return "double"

    @PYB11const
    def prime(self,
              x = "const double"):
        "Interpolator for the first derivative: <dy/dx>(x)"
        return "double"

    @PYB11const
    def prime2(self,
               x = "const double"):
        "Interpolator for the second derivative: <d^2y/dx^2>(x)"
        return "double"

    @PYB11pyname("__call__")
    @PYB11const
    def __call__i0(self,
                   x = "const double",
                   i0 = "const size_t"):
        "Returns the interpolated value <y>(x)"
        return "double"

    @PYB11pycppname("prime")
    @PYB11const
    def prime_i0(self,
                 x = "const double",
                 i0 = "const size_t"):
        "Interpolator for the first derivative: <dy/dx>(x)"
        return "double"

    @PYB11pycppname("prime2")
    @PYB11const
    def prime2_i0(self,
                  x = "const double",
                  i0 = "const size_t"):
        "Interpolator for the second derivative: <d^2y/dx^2>(x)"
        return "double"

    @PYB11const
    def lowerBound(self,
                   x = "const double"):
        "Return the lower bound index in the table for the given x coordinate"
        return "size_t"

    # Attributes
    size = PYB11property(doc="The size of the tabulated coefficient arrays")
    xmin = PYB11property(doc="Minimum x coordinate for table")
    xmax = PYB11property(doc="Maximum x coordinate for table")
    xstep = PYB11property(doc="delta x between tabulated values")
    coeffs = PYB11property(doc="the fitting coefficients")
