"""
PolyClipper module.

Binds the PolyClipper geometry operations for clipping polygons & polyhedra.
"""

from PYB11Generator import *
from CXXTypesMOD import *
import types

# Include files.
includes = ['"Geometry/polyclipper.hh"']

namespaces = ["PolyClipper"]

#-------------------------------------------------------------------------------
# Helper to add methods to Planes.
#-------------------------------------------------------------------------------
@PYB11ignore
def addPlaneMethods(cls, ndim):

    # Constructors
    def pyinit0(self):
        "Default constructor"

    def pyinit1(self,
                rhs = "const Plane%id" % ndim):
        "Copy constructor"

    def pyinit2(self,
                dist = "double",
                normal = "const Plane%id::Vector&" % ndim):
        "Construct with a distance and normal."

    def pyinit3(self,
                point = "const Plane%id::Vector&" % ndim,
                normal = "const Plane%id::Vector&" % ndim):
        "Construct with a point and normal."

    # Attributes
    dist = PYB11readwrite(doc="The distance to the origin along the normal.")
    normal = PYB11readwrite(doc="The normal to the plane.")

    for __x in [x for x in dir() if type(eval(x)) == types.FunctionType]: 
        exec("cls.%s = %s" % (__x, __x))

#-------------------------------------------------------------------------------
# Helper to add methods for Vertex.
#-------------------------------------------------------------------------------
@PYB11ignore
class VertexBase:

    # Constructors
    def pyinit0(self):
        "Default constructor"

    def pyinit1(rhs = "const Vertex%(ndim)sd"):
        "Copy constructor"

    def pyinit2(position = "const Plane%(ndim)sd::Vector&"):
        "Construct with a position."

    def pyinit3(position = "const Plane%(ndim)sd::Vector&",
                c = "int"):
        "Construct with a position and initial compare flag."

    # Attributes
    position = PYB11readwrite(doc="The position of the vertex.")
    neighbors = PYB11readwrite(doc="The connectivty of the vertex.")
    comp = PYB11readwrite(doc="The current comparison flag.")
    ID = PYB11readwrite(doc="The ID or index of the vertex.")

    def __eq__(self, other="py::self"):
        return

#-------------------------------------------------------------------------------
# Plane2d
#-------------------------------------------------------------------------------
@PYB11cppname("Plane2d")
class PolyClipperPlane2d:
    """Plane class for polyclipper in 2 dimensions."""

addPlaneMethods(PolyClipperPlane2d, 2)

#-------------------------------------------------------------------------------
# Plane3d
#-------------------------------------------------------------------------------
@PYB11cppname("Plane3d")
class PolyClipperPlane3d:
    """Plane class for polyclipper in 3 dimensions."""

addPlaneMethods(PolyClipperPlane3d, 3)

#-------------------------------------------------------------------------------
# Vertex2d
#-------------------------------------------------------------------------------
@PYB11template_dict({"ndim" : "2"})
class Vertex2d:
    """Vertex class for polyclipper in 2 dimensions."""
PYB11inject(VertexBase, Vertex2d)

#-------------------------------------------------------------------------------
# Vertex3d
#-------------------------------------------------------------------------------
@PYB11template_dict({"ndim" : "3"})
class Vertex3d:
    """Vertex class for polyclipper in 3 dimensions."""
PYB11inject(VertexBase, Vertex3d)

#-------------------------------------------------------------------------------
# Polygon & Polyhedron
#-------------------------------------------------------------------------------
Polygon    = PYB11_bind_vector("PolyClipper::Vertex2d", opaque=True, local=False)
Polyhedron = PYB11_bind_vector("PolyClipper::Vertex3d", opaque=True, local=False)

#-------------------------------------------------------------------------------
# Polygon methods.
#-------------------------------------------------------------------------------
@PYB11namespace("PolyClipper")
def initializePolygon(poly = "Polygon&",
                      positions = "const std::vector<Spheral::Dim<2>::Vector>&",
                      neighbors = "const std::vector<std::vector<int>>&"):
    "Initialize a PolyClipper::Polygon from vertex positions and vertex neighbors."
    return "void"

@PYB11namespace("PolyClipper")
def polygon2string(poly = "Polygon&"):
    "Return a formatted string representation for a PolyClipper::Polygon."
    return "std::string"

@PYB11namespace("PolyClipper")
def convertToPolygon(polygon = "Polygon&",
                     Spheral_polygon = "const Spheral::Dim<2>::FacetedVolume&"):
    "Construct a PolyClipper::Polygon from a Spheral::Polygon."
    return "void"

@PYB11namespace("PolyClipper")
def convertFromPolygon(Spheral_polygon = "Spheral::Dim<2>::FacetedVolume&",
                       polygon = "const Polygon&"):
    "Construct a Spheral::Polygon from a PolyClipper::Polygon."
    return "void"

@PYB11namespace("PolyClipper")
@PYB11implementation("""[](const Polygon& self) {
                                                  double zerothMoment;
                                                  Spheral::Dim<2>::Vector firstMoment;
                                                  moments(zerothMoment, firstMoment, self);
                                                  return py::make_tuple(zerothMoment, firstMoment);
                                                }""")
@PYB11pycppname("moments")
def momentsPolygon(poly = "const Polygon&"):
    "Compute the zeroth and first moment of a PolyClipper::Polygon."
    return "py::tuple"

@PYB11namespace("PolyClipper")
def clipPolygon(poly = "Polygon&",
                planes = "const std::vector<Plane2d>&"):
    "Clip a PolyClipper::Polygon with a collection of planes."
    return "void"

@PYB11namespace("PolyClipper")
@PYB11pycppname("collapseDegenerates")
def collapseDegeneratesPolygon(poly = "Polygon&",
                               tol = "const double"):
    "Collapse edges in a PolyClipper::Polygon below the given tolerance."
    return "void"

@PYB11namespace("PolyClipper")
def splitIntoTriangles(poly = "const Polygon&",
                       tol = ("const double", "0.0")):
    "Split a PolyClipper::Polygon into triangles.\n"
    "The result is returned as a vector<vector<int>>, where each inner vector is a triple of\n"
    "ints representing vertex indices in the input Polygon."
    return "std::vector<std::vector<int>>"

#-------------------------------------------------------------------------------
# Polyhedron methods.
#-------------------------------------------------------------------------------
@PYB11namespace("PolyClipper")
def initializePolyhedron(poly = "Polyhedron&",
                         positions = "const std::vector<Spheral::Dim<3>::Vector>&",
                         neighbors = "const std::vector<std::vector<int>>&"):
    "Initialize a PolyClipper::Polyhedron from vertex positions and vertex neighbors."

@PYB11namespace("PolyClipper")
def polyhedron2string(poly = "Polyhedron&"):
    "Return a formatted string representation for a PolyClipper::Polyhedron."

@PYB11namespace("PolyClipper")
def convertToPolyhedron(polyhedron = "Polyhedron&",
                        Spheral_polyhedron = "const Spheral::Dim<3>::FacetedVolume&"):
    "Construct a PolyClipper::Polyhedron from a Spheral::Polyhedron."

@PYB11namespace("PolyClipper")
def convertFromPolyhedron(Spheral_polyhedron = "Spheral::Dim<3>::FacetedVolume&",
                          polyhedron = "const Polyhedron&"):
    "Construct a Spheral::Polyhedron from a PolyClipper::Polyhedron."

@PYB11namespace("PolyClipper")
def moments(zerothMoment = "double&",
            firstMoment = "Spheral::Dim<3>::Vector&",
            poly = "const Polyhedron&"):
    "Compute the zeroth and first moment of a PolyClipper::Polyhedron."
    return "void"

@PYB11namespace("PolyClipper")
def clipPolyhedron(poly = "Polyhedron&",
                   planes = "const std::vector<Plane3d>&"):
    "Clip a PolyClipper::Polyhedron with a collection of planes."

@PYB11namespace("PolyClipper")
def collapseDegenerates(poly = "Polyhedron&",
                        tol = "const double"):
    "Collapse edges in a PolyClipper::Polyhedron below the given tolerance."
    return "void"

@PYB11namespace("PolyClipper")
def splitIntoTetrahedra(poly = "const Polyhedron&",
                        tol = ("const double", "0.0")):
    "Split a PolyClipper::Polyhedron into tetrahedra.\n"
    "The result is returned as a vector<vector<int>>, where each inner vector is a set of four\n"
    "ints representing vertex indices in the input Polyhedron."
    return "std::vector<std::vector<int>>"
