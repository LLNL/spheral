from PYB11Generator import *

#-------------------------------------------------------------------------------
# Plane template
#-------------------------------------------------------------------------------
@PYB11template("ndim")
class Plane:
    "The geometric representation for a plane in %(ndim)s dimensions."

    PYB11typedefs = """
    typedef Dim<%(ndim)s>::Vector Vector;
    typedef GeomPlane<Dim<%(ndim)s>> PlaneType;
"""

    # Constructors
    def pyinit0(self):
        "Default constructor"

    def pyinit1(self,
                rhs = "const PlaneType&"):
        "Copy constructor"

    def pyinit2(self,
                point = "const Vector&",
                normal = "const Vector&"):
        "Construct from a (point,normal)."

    def pyinit3(self,
                points = "const std::vector<Vector>&"):
        "Best fit plane for a set of point."

    # Methods
    def signedDistance(self):
        "The signed distance of the given point to the plane."

    def minimumDistance(self):
        "The distance of the given point from the plane."

    def closestPointOnPlane(self):
        "Find the point in the plane closest to the given point."

    def parallel(self):
        "Test if the given plane is parallel to this one."

    def compare(self):
        "Test if the plane is below, in, or above the given point (-1, 0, 1)."

    # Comparisons
    @PYB11pyname("__eq__")
    def __eq__self(self):
        return
    @PYB11pyname("__ne__")
    def __ne__self(self):
        return
    @PYB11pyname("__lt__")
    def __lt__self(self):
        return

    def __eq__(self, rhs="Vector()"):
        return
    def __ne__(self, rhs="Vector()"):
        return
    def __lt__(self, rhs="Vector()"):
        return
    def __gt__(self, rhs="Vector()"):
        return
    def __le__(self, rhs="Vector()"):
        return
    def __ge__(self, rhs="Vector()"):
        return

    # Properties
    point = PYB11property("const Vector&", "point", "point", doc="The point for plane definition.")
    normal = PYB11property("const Vector&", "normal", "normal", doc="The unit normal to the plane.")

    # String representation
    @PYB11implementation("""
[](const GeomPlane<Dim<%(ndim)s>>& self) {
  std::ostringstream ss;
  ss << "Plane%(ndim)sd[point=" << self.point() << " normal=" << self.normal() << "]";
  return ss.str();
}""")
    def __repr__(self):
        return

#-------------------------------------------------------------------------------
# Plane instantiations.
#-------------------------------------------------------------------------------
Plane1d = PYB11TemplateClass(Plane,
                                   template_parameters = ("1"),
                                   cppname = "GeomPlane<Dim<1>>",
                                   pyname = "Plane1d")
Plane2d = PYB11TemplateClass(Plane,
                                   template_parameters = ("2"),
                                   cppname = "GeomPlane<Dim<2>>",
                                   pyname = "Plane2d")
Plane3d = PYB11TemplateClass(Plane,
                                   template_parameters = ("3"),
                                   cppname = "GeomPlane<Dim<3>>",
                                   pyname = "Plane3d")
