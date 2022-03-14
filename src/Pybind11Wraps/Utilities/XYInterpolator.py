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

    # Attributes
    xmin = PYB11property(doc="Minimum x coordinate for table")
    xmax = PYB11property(doc="Maximum x coordinate for table")
    ymin = PYB11property(doc="Minimum y coordinate for table")
    ymax = PYB11property(doc="Maximum y coordinate for table")
    xstep = PYB11property(doc="delta x between tabulated values")
    ystep = PYB11property(doc="delta y between tabulated values")
    xlog = PYB11property(doc="Use logarithmic spacing in x")
    ylog = PYB11property(doc="Use logarithmic spacing in y")
