#-------------------------------------------------------------------------------
# BiQuadraticInterpolator
#-------------------------------------------------------------------------------
from PYB11Generator import *

class BiQuadraticInterpolator:

    """Encapsulates the algorithm and data for parabolic interpolation in 2D
Assumes the results is interpolated as
  <F(x,y)> = c0 + c1*x + c2*y + c3*x*y + c4*x^2 + c5*y^2"""

    PYB11typedefs = """
    typedef BiQuadraticInterpolator::Vector Vector;
"""

    def pyinit(self):
        "Default constuctor -- returns a non-functional interpolator until initialized"
        return

    def pyinit1(self,
                xmin = "const Vector&",
                xmax = "const Vector&",
                nx = "const size_t",
                ny = "const size_t",
                logxspace = "const bool",
                logyspace = "const bool",
                F = "const Spheral::PythonBoundFunctors::SpheralFunctor<Vector, double>&"):
        "Returns an interpolator for yvals sampled in x in [xmin, xmax]"
        return

    @PYB11template("Func")
    @PYB11cppname("initialize")
    def initialize_(self,
                    xmin = "const Vector&",
                    xmax = "const Vector&",
                    nx = "const size_t",
                    ny = "const size_t",
                    logxspace = "const bool",
                    logyspace = "const bool",
                    F = "const %(Func)s&"):
        "Initializes the interpolator for interpolating the given function"
        return "void"

    initialize = PYB11TemplateMethod(initialize_, "Spheral::PythonBoundFunctors::SpheralFunctor<Vector, double>")

    def __call__(self,
                 pos = "const Vector&"):
        "Returns the interpolated value <F>(x,y)"
        return "double"

    # Attributes
    size = PYB11property(doc="The size of the tabulated coefficient arrays")
    xmin = PYB11property(doc="Minimum coordinate for table")
    xmax = PYB11property(doc="Maximum coordinate for table")
    xstep = PYB11property(doc="delta x between tabulated values")
    xlog = PYB11property(doc="Fitting in log x space?")
    ylog = PYB11property(doc="Fitting in log y space?")
    coeffs = PYB11property(doc="the fitting coefficients")
