#-------------------------------------------------------------------------------
# XYInterpolator
#-------------------------------------------------------------------------------
from PYB11Generator import *

class XYInterpolator:

    "Base class for 2D interpolator classes in Spheral."

    def pyinit(self):
        "Default constuctor"
        return

    def pyinit1(self,
                order = "const size_t",
                xmin = "const double",
                xmax = "const double",
                ymin = "const double",
                ymax = "const double",
                nx = "const size_t",
                ny = "const size_t",
                xlog = ("const bool", "false"),
                ylog = ("const bool", "false")):
        "Returns an interpolator for yvals sampled in x in [xmin, xmax]"
        return

    @PYB11implementation("""[](const XYInterpolator& self, const double x, const double y) -> py::tuple {
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

    @PYB11implementation("""[](const XYInterpolator& self, const double x, const double y) -> py::tuple {
                                 double etax, etay;
                                 size_t ix, iy, i0;
                                 self.eta_coords(x, y, etax, etay, ix, iy, i0);
                                 return py::make_tuple(etax, etay, ix, iy, i0);
                            }""")
    @PYB11const
    def eta_coords(self,
                   x = "const double",
                   y = "const double"):
        "Return the normalize (etax, etay) coords as well as the result of lowerBound"
        return "py::tuple"

    @PYB11const
    def coord(self):
        "Compute a coordinate value depending on whether we're using log-space"
        return

    @PYB11const
    def xcoord(self):
        "X specialization of coord"
        return

    @PYB11const
    def ycoord(self):
        "Y specialization of coord"
        return

    # Attributes
    xmin = PYB11property(doc="Minimum x coordinate for table")
    xmax = PYB11property(doc="Maximum x coordinate for table")
    ymin = PYB11property(doc="Minimum y coordinate for table")
    ymax = PYB11property(doc="Maximum y coordinate for table")
    xstep = PYB11property(doc="delta x between tabulated values")
    ystep = PYB11property(doc="delta y between tabulated values")
    Ax = PYB11property(doc="A for x log stepping")
    Bx = PYB11property(doc="B for x log stepping")
    Ay = PYB11property(doc="A for y log stepping")
    By = PYB11property(doc="B for y log stepping")
    xlog = PYB11property(doc="Use logarithmic spacing in x")
    ylog = PYB11property(doc="Use logarithmic spacing in y")
    size = PYB11property(doc="The size of the tabulated coefficient arrays")
    coeffs = PYB11property(doc="the fitting coefficients")
