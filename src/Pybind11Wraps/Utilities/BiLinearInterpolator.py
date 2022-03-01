#-------------------------------------------------------------------------------
# BiLinearInterpolator
#-------------------------------------------------------------------------------
from PYB11Generator import *

class BiLinearInterpolator:

    """Encapsulates the algorithm and data for bi-linear interpolation in 2D
Assumes the results is interpolated as
  <F(x,y)> = c0 + c1*x + c2*y + c3*x*y"""

    def pyinit(self):
        "Default constuctor -- returns a non-functional interpolator until initialized"
        return

    def pyinit1(self,
                xmin = "const double",
                xmax = "const double",
                ymin = "const double",
                ymax = "const double",
                nx = "const size_t",
                ny = "const size_t",
                F = "const Spheral::PythonBoundFunctors::Spheral2ArgFunctor<double, double, double>&"):
        "Returns an interpolator for yvals sampled in x in [xmin, xmax]"
        return

    @PYB11template("Func")
    @PYB11cppname("initialize")
    def initialize_(self,
                    xmin = "const double",
                    xmax = "const double",
                    ymin = "const double",
                    ymax = "const double",
                    nx = "const size_t",
                    ny = "const size_t",
                    F = "const %(Func)s&"):
        "Initializes the interpolator for interpolating the given function"
        return "void"

    initialize = PYB11TemplateMethod(initialize_, "Spheral::PythonBoundFunctors::Spheral2ArgFunctor<double, double, double>")

    def __call__(self,
                 x = "const double",
                 y = "const double"):
        "Returns the interpolated value <F>(x,y)"
        return "double"

    def prime_x(self,
                x = "const double",
                y = "const double"):
        """Returns the interpolated value \\\partial_x <F>(x,y)"""
        return "double"

    def prime_y(self,
                x = "const double",
                y = "const double"):
        """Returns the interpolated value \\\partial_y <F>(x,y)"""
        return "double"

    def lowerBound(self,
                   x = "const double",
                   y = "const double"):
        "Return the index into the coefficient array for the given coordinates"
        return "size_t"

    # Attributes
    size = PYB11property(doc="The size of the tabulated coefficient arrays")
    xmin = PYB11property(doc="Minimum x coordinate for table")
    xmax = PYB11property(doc="Maximum x coordinate for table")
    ymin = PYB11property(doc="Minimum y coordinate for table")
    ymax = PYB11property(doc="Maximum y coordinate for table")
    xstep = PYB11property(doc="delta x between tabulated values")
    ystep = PYB11property(doc="delta y between tabulated values")
    coeffs = PYB11property(doc="the fitting coefficients")
