#-------------------------------------------------------------------------------
# BiLinearInterpolator
#-------------------------------------------------------------------------------
from PYB11Generator import *
from XYInterpolator import *

class BiLinearInterpolator(XYInterpolator):

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
                F = "const Spheral::PythonBoundFunctors::Spheral2ArgFunctor<double, double, double>&",
                xlog = ("const bool", "false"),
                ylog = ("const bool", "false")):
        "Returns an interpolator for ([xmin, xmax], [ymin, ymax])"
        return

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
