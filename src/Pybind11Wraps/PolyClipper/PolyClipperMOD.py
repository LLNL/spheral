"""
PolyClipper module.

Binds the PolyClipper geometry operations for clipping polygons & polyhedra.
"""

from PYB11Generator import *
import types

# Include files.
PYB11includes = ['"Geometry/polyclipper.hh"']

PYB11namespaces = ["PolyClipper"]

PYB11opaque = ["std::vector<char>",
               "std::vector<unsigned>",
               "std::vector<uint64_t>",
               "std::vector<int>",
               "std::vector<float>",
               "std::vector<double>",
               "std::vector<std::string>",

               "std::vector<std::vector<char>>",
               "std::vector<std::vector<unsigned>>",
               "std::vector<std::vector<uint64_t>>",
               "std::vector<std::vector<int>>",
               "std::vector<std::vector<float>>",
               "std::vector<std::vector<double>>",
               "std::vector<std::vector<std::string>>"]

#-------------------------------------------------------------------------------
# Planes.
#-------------------------------------------------------------------------------
@PYB11ignore
class PlaneBase:

    # Constructors
    def pyinit0(self):
        "Default constructor"

    def pyinit1(self,
                rhs = "const Plane%(ndim)sd"):
        "Copy constructor"

    def pyinit2(self,
                dist = "double",
                normal = "const Plane%(ndim)sd::Vector&"):
        "Construct with a distance and normal."

    def pyinit3(self,
                point = "const Plane%(ndim)sd::Vector&",
                normal = "const Plane%(ndim)sd::Vector&"):
        "Construct with a point and normal."

    def pyinit4(self,
                point = "const Plane%(ndim)sd::Vector&",
                normal = "const Plane%(ndim)sd::Vector&",
                id = "const int"):
        "Construct with a point, normal, and ID."

    # Attributes
    dist = PYB11readwrite(doc="The distance to the origin along the normal.")
    normal = PYB11readwrite(doc="The normal to the plane.")
    ID = PYB11readwrite(doc="Arbitrary ID number for plane")

@PYB11template_dict({"ndim" : "2"})
class Plane2d:
    """Plane class for polyclipper in 2 dimensions."""
PYB11inject(PlaneBase, Plane2d)

@PYB11template_dict({"ndim" : "3"})
class Plane3d:
    """Plane class for polyclipper in 3 dimensions."""
PYB11inject(PlaneBase, Plane3d)

#-------------------------------------------------------------------------------
# Vertex
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
    clips = PYB11readwrite(doc="The set of plane IDs (if any) responsible for this vertex.")

    def __eq__(self):
        return

@PYB11template_dict({"ndim" : "2"})
class Vertex2d:
    """Vertex class for polyclipper in 2 dimensions."""
PYB11inject(VertexBase, Vertex2d)

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
    "Construct a Spheral::Polygon from a PolyClipper::Polygon.  Returns the set of clip planes responsible for each vertex."
    return "std::vector<std::set<int>>"

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
    """Split a PolyClipper::Polygon into triangles.
The result is returned as a vector<vector<int>>, where each inner vector is a triple of
ints representing vertex indices in the input Polygon."""
    return "std::vector<std::vector<int>>"

#-------------------------------------------------------------------------------
# Polyhedron methods.
#-------------------------------------------------------------------------------
@PYB11namespace("PolyClipper")
def initializePolyhedron(poly = "Polyhedron&",
                         positions = "const std::vector<Spheral::Dim<3>::Vector>&",
                         neighbors = "const std::vector<std::vector<int>>&"):
    "Initialize a PolyClipper::Polyhedron from vertex positions and vertex neighbors."
    return "void"

@PYB11namespace("PolyClipper")
def polyhedron2string(poly = "Polyhedron&"):
    "Return a formatted string representation for a PolyClipper::Polyhedron."
    return "std::string"

@PYB11namespace("PolyClipper")
def convertToPolyhedron(polyhedron = "Polyhedron&",
                        Spheral_polyhedron = "const Spheral::Dim<3>::FacetedVolume&"):
    "Construct a PolyClipper::Polyhedron from a Spheral::Polyhedron."
    return "void"

@PYB11namespace("PolyClipper")
def convertFromPolyhedron(Spheral_polyhedron = "Spheral::Dim<3>::FacetedVolume&",
                          polyhedron = "const Polyhedron&"):
    "Construct a Spheral::Polyhedron from a PolyClipper::Polyhedron.  Returns the set of clip planes responsible for each vertex."
    return "std::vector<std::set<int>>"

@PYB11namespace("PolyClipper")
@PYB11implementation("""[](const Polyhedron& self) {
                                                     double zerothMoment;
                                                     Spheral::Dim<3>::Vector firstMoment;
                                                     moments(zerothMoment, firstMoment, self);
                                                     return py::make_tuple(zerothMoment, firstMoment);
                                                   }""")
@PYB11pycppname("moments")
def momentsPolyhedron(poly = "const Polyhedron&"):
    "Compute the zeroth and first moment of a PolyClipper::Polyhedron."
    return "py::tuple"

@PYB11namespace("PolyClipper")
def clipPolyhedron(poly = "Polyhedron&",
                   planes = "const std::vector<Plane3d>&"):
    "Clip a PolyClipper::Polyhedron with a collection of planes."
    return "void"

@PYB11namespace("PolyClipper")
@PYB11pycppname("collapseDegenerates")
def collapseDegeneratesPolyhedron(poly = "Polyhedron&",
                                  tol = "const double"):
    "Collapse edges in a PolyClipper::Polyhedron below the given tolerance."
    return "void"

@PYB11namespace("PolyClipper")
def splitIntoTetrahedra(poly = "const Polyhedron&",
                        tol = ("const double", "0.0")):
    """Split a PolyClipper::Polyhedron into tetrahedra.
The result is returned as a vector<vector<int>>, where each inner vector is a set of four
ints representing vertex indices in the input Polyhedron."""
    return "std::vector<std::vector<int>>"

