#-------------------------------------------------------------------------------
# Facet2d
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11cppname("GeomFacet2d")
class Facet2d:
    """GeomFacet2d -- A facet of a polygon (just two points).

Note a Facet does not maintain it's own copies of it's end points -- the
assumption is that this is a Facet of a GeomPolygon and that polygon owns
the set of vertex positions."""

    PYB11typedefs = """
    typedef GeomFacet2d::Vector Vector;
"""

    #...........................................................................
    # Constructors
    def pyinit(self):
        "Default constructor"

    def pyinit1(self,
                vertices = "const std::vector<Vector>&",
                point1 = "const unsigned",
                point2 = "const unsigned"):
        "Explicit constructor with vertices and edge point indices"

    #...........................................................................
    # Methods
    @PYB11const
    @PYB11pycppname("compare")
    def compare0(self,
                 point = "const Vector&",
                 tol = ("const double tol", "1.0e-8")):
        """Is the given point above, below, or colinear with the facet?
  1 => point above.
  0 => point in plane of facet
 -1 => points below."""
        return "int"

    @PYB11const
    @PYB11pycppname("compare")
    def compare1(self,
                 points = "const std::vector<Vector>&",
                 tol = ("const double tol", "1.0e-8")):
        """Compare a set of points:
  1 => all points above.
  0 => points both above and below (or equal).
 -1 => all points below."""
        return "int"

    @PYB11const
    def distance(self,
                 p = "const Vector&"):
        "Compute the minimum distance from the facet to a point."
        return "double"

    @PYB11const
    def closestPoint(self,
                     p = "const Vector&"):
        "Compute the closest point on the facet to the given point."
        return "Vector"

    #...........................................................................
    # Comparisons
    def __eq__(self):
        return

    def __ne__(self):
        return

    #...........................................................................
    # Properties
    point1 = PYB11property("const Vector&", returnpolicy="reference_internal")
    point2 = PYB11property("const Vector&", returnpolicy="reference_internal")
    ipoint1 = PYB11property("unsigned")
    ipoint2 = PYB11property("unsigned")
    ipoints = PYB11property("const std::vector<unsigned>&", returnpolicy="reference_internal")
    normal = PYB11property("const Vector&", returnpolicy="reference_internal")
    position = PYB11property("Vector")
    area = PYB11property("double")
