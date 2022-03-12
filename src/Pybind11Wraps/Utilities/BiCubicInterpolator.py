#-------------------------------------------------------------------------------
# BiCubicInterpolator
#-------------------------------------------------------------------------------
from PYB11Generator import *

class BiCubicInterpolator:

    """BiCubicInterpolator

Encapsulates the algorithm and data for cubic interpolation in 2D
Assumes the results is interpolated as
  <F(x,y)> = [1 x x^2 x^3][c00 c01 c02 c03][1]
                          [c10 c11 c12 c13][y]
                          [c20 c21 c22 c23][y^2]
                          [c30 c31 c32 c33][y^3]
Assumes we provide functors for F and gradF for fitting, where
     F(x,y) -> double
 gradF(x,y) -> Dim<2>::SymTensor

See https://en.wikipedia.org/wiki/Bicubic_interpolation"""

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
                gradF = "const Spheral::PythonBoundFunctors::Spheral2ArgFunctor<double, double, Dim<2>::SymTensor>&",
                xlog = ("const bool", "false"),
                ylog = ("const bool", "false")):
        "Returns an interpolator for z-vals sampled in (x,y) in ([xmin, xmax], [ymin, ymax])"
        return

    @PYB11template("Func", "GradFunc")
    @PYB11cppname("initialize")
    def initialize_(self,
                    xmin = "const double",
                    xmax = "const double",
                    ymin = "const double",
                    ymax = "const double",
                    nx = "const size_t",
                    ny = "const size_t",
                    F = "const %(Func)s&",
                    gradF = "const %(GradFunc)s&",
                    xlog = ("const bool", "false"),
                    ylog = ("const bool", "false")):
        "Initializes the interpolator for interpolating the given function"
        return "void"

    initialize = PYB11TemplateMethod(initialize_, ("Spheral::PythonBoundFunctors::Spheral2ArgFunctor<double, double, double>",
                                                   "Spheral::PythonBoundFunctors::Spheral2ArgFunctor<double, double, Dim<2>::SymTensor>"))

    @PYB11const
    def __call__(self,
                 x = "const double",
                 y = "const double"):
        "Returns the interpolated value <F>(x,y)"
        return "double"

    @PYB11const
    def prime_x(self,
                x = "const double",
                y = "const double"):
        """Returns the interpolated value \\\partial_x <F>(x,y)"""
        return "double"

    @PYB11const
    def prime_y(self,
                x = "const double",
                y = "const double"):
        """Returns the interpolated value \\\partial_y <F>(x,y)"""
        return "double"

    @PYB11const
    def prime2_xx(self,
                  x = "const double",
                  y = "const double"):
        """Returns the interpolated value \\\partial_xx <F>(x,y)"""
        return "double"

    @PYB11const
    def prime2_xy(self,
                  x = "const double",
                  y = "const double"):
        """Returns the interpolated value \\\partial_xy <F>(x,y)"""
        return "double"

    @PYB11const
    def prime2_yy(self,
                  x = "const double",
                  y = "const double"):
        """Returns the interpolated value \\\partial_yy <F>(x,y)"""
        return "double"

    @PYB11implementation("""[](const BiCubicInterpolator& self, const double x, const double y) -> py::tuple {
                                 size_t ix, iy, i0;
                                 self.lowerBound(x, y, ix, iy, i0);
                                 return py::make_tuple(ix, iy, i0);
                            }""")
    @PYB11const
    def lowerBound(self,
                   x = "const double",
                   y = "const double"):
        "Return the index into the coefficient array for the given coordinates"
        return "py::tuple"

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
