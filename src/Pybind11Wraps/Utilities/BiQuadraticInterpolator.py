#-------------------------------------------------------------------------------
# BiQuadraticInterpolator
#-------------------------------------------------------------------------------
from PYB11Generator import *

class BiQuadraticInterpolator:

    """Encapsulates the algorithm and data for parabolic interpolation in 2D
Assumes the results is interpolated as
  <F(x,y)> = c0 + c1*x + c2*y + c3*x*y + c4*x^2 + c5*y^2"""

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
                F = "const Spheral::PythonBoundFunctors::Spheral2ArgFunctor<double, double, double>&",
                xlog = ("const bool", "false"),
                ylog = ("const bool", "false")):
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
                    F = "const %(Func)s&",
                    xlog = ("const bool", "false"),
                    ylog = ("const bool", "false")):
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

    def prime2_xx(self,
                  x = "const double",
                  y = "const double"):
        """Returns the interpolated value \\\partial_xx <F>(x,y)"""
        return "double"

    def prime2_xy(self,
                  x = "const double",
                  y = "const double"):
        """Returns the interpolated value \\\partial_xy <F>(x,y)"""
        return "double"

    def prime2_yx(self,
                  x = "const double",
                  y = "const double"):
        """Returns the interpolated value \\\partial_yx <F>(x,y)"""
        return "double"

    def prime2_yy(self,
                  x = "const double",
                  y = "const double"):
        """Returns the interpolated value \\\partial_yy <F>(x,y)"""
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
    xlog = PYB11property(doc="Use logarithmic spacing in x")
    ylog = PYB11property(doc="Use logarithmic spacing in y")
