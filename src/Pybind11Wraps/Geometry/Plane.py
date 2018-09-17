from PYB11Decorators import *
from PYB11property import *
from PYB11class import *

#-------------------------------------------------------------------------------
# Plane template
#-------------------------------------------------------------------------------
@PYB11template("ndim")
class Plane:
    "The geometric representation for a plane in %(ndim)s dimensions."

    # Constructors
    def pyinit0(self):
        "Default constructor"

    def pyinit1(self,
                rhs = "const GeomPlane<Dim<%(ndim)s>>"):
        "Copy constructor"

    def pyinit2(self,
                point = "const Dim<%(ndim)s>::Vector&",
                normal = "const Dim<%(ndim)s>::Vector&"):
        "Construct from a (point,normal)."

    def pyinit3(self,
                points = "const std::vector<Dim<%(ndim)s>::Vector>&"):
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
    def __eq__(self):
        return
    def __ne__(self):
        return
    def __lt__(self):
        return
    def __eq__(self, rhs="Dim<%(ndim)s>::Vector()"):
        return
    def __ne__(self, rhs="Dim<%(ndim)s>::Vector()"):
        return
    def __lt__(self, rhs="Dim<%(ndim)s>::Vector()"):
        return
    def __gt__(self, rhs="Dim<%(ndim)s>::Vector()"):
        return
    def __le__(self, rhs="Dim<%(ndim)s>::Vector()"):
        return
    def __ge__(self, rhs="Dim<%(ndim)s>::Vector()"):
        return

    # Point
    @PYB11cppname("point")
    @PYB11ignore
    @PYB11const
    def getpoint(self):
        return "const Dim<%(ndim)s>::Vector&"

    @PYB11cppname("point")
    @PYB11ignore
    def setpoint(self,
                 val = "const Dim<%(ndim)s>::Vector&"):
        return "void"

    # Normal
    @PYB11cppname("normal")
    @PYB11ignore
    @PYB11const
    def getnormal(self):
        return "const Dim<%(ndim)s>::Vector&"

    @PYB11cppname("normal")
    @PYB11ignore
    def setnormal(self,
                 val = "const Dim<%(ndim)s>::Vector&"):
        return "void"

    # Properties
    point = property(getpoint, setpoint, doc="The point for plane definition.")
    normal = property(getnormal, setnormal, doc="The unit normal to the plane.")

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
