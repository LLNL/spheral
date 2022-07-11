#-------------------------------------------------------------------------------
# CubicHermiteInterpolator
#-------------------------------------------------------------------------------
from PYB11Generator import *

class CubicHermiteInterpolator:

    "An (optionally monotonic) form of cubic Hermite interpolation."

    def pyinit(self):
        "Default constuctor -- returns a non-functional interpolator until initialized"
        return

    def pyinit_func(self,
                    xmin = "const double",
                    xmax = "const double",
                    n = "const size_t",
                    F = "const PythonBoundFunctors::SpheralFunctor<double, double>&"):
        "Constructs an interpolator based on the given function"
        return

    def pyinit_gradfunc(self,
                        xmin = "const double",
                        xmax = "const double",
                        n = "const size_t",
                        F = "const PythonBoundFunctors::SpheralFunctor<double, double>&",
                        Fgrad = "const PythonBoundFunctors::SpheralFunctor<double, double>&"):
        "Constructs an interpolator based on the given function and its gradient"
        return

    def initialize(self,
                   xmin = "const double",
                   xmax = "const double",
                   yvals = "const std::vector<double>&"):
        "Initializes the interpolator for yvals sampled in x in [xmin, xmax]"

    def makeMonotonic(self):
        """Force interpolation to be monotonic.  This generally degrades accuracy, and can introduce structure between
tabulated knots (although that structure should still be monotonic between the tabulated function values)."""
        return "void"

    @PYB11const
    def __call__(self,
                 x = "const double"):
        "Returns the interpolated value <y>(x)"
        return "double"

    @PYB11pyname("__call__")
    @PYB11const
    def __call__i0(self,
                   x = "const double",
                   i0 = "const size_t"):
        "Returns the interpolated value <y>(x)"
        return "double"

    @PYB11const
    def lowerBound(self,
                   x = "const double"):
        "Return the lower bound index in the table for the given x coordinate"
        return "size_t"

    @PYB11const
    def h00(self,
            t = "const double"):
        "h00 Hermite basis function"
        return "double"

    @PYB11const
    def h10(self,
            t = "const double"):
        "h10 Hermite basis function"
        return "double"

    @PYB11const
    def h01(self,
            t = "const double"):
        "h01 Hermite basis function"
        return "double"

    @PYB11const
    def h11(self,
            t = "const double"):
        "h11 Hermite basis function"
        return "double"

    # Attributes
    N = PYB11property(doc="The number of the tabulated values used")
    xmin = PYB11property(doc="Minimum x coordinate for table")
    xmax = PYB11property(doc="Maximum x coordinate for table")
    xstep = PYB11property(doc="delta x between tabulated values")
    vals = PYB11property(doc="tabulated values and gradients (sequentially)")
