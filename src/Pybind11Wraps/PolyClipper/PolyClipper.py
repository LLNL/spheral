"""
PolyClipper module.

Binds the PolyClipper geometry operations for clipping polygons & polyhedra.
"""

from PYB11Decorators import *
import types

# Include files.
includes = ['"Geometry/polyclipper.hh"']

#-------------------------------------------------------------------------------
# Helper to add the templated methods to Planes.
#-------------------------------------------------------------------------------
@PYB11ignore
def addPlaneMethods(cls, ndim):

    # Constructors
    def pyinit0(self):
        "Default constructor"
        return

    def pyinit1(rhs = "const PolyClipper::Plane%id" % ndim):
        "Copy constructor"
        return

    def pyinit2(dist = "double",
                normal = "const PolyClipper::Plane%id::Vector&" % ndim):
        "Construct with a distance and normal."
        return

    def pyinit3(point = "const PolyClipper::Plane%id::Vector&" % ndim,
                normal = "const PolyClipper::Plane%id::Vector&" % ndim):
        "Construct with a point and normal."
        return

    # Attributes
    @PYB11readwrite
    def dist(self):
        "The distance to the origin along the normal."
        return "double"

    @PYB11readwrite
    def normal(self):
        "The normal to the plane."
        return "Vector"

    for x in [x for x in dir() if type(eval(x)) == types.FunctionType]: 
        exec("cls.%s = %s" % (x, x))

#-------------------------------------------------------------------------------
# Helper to add the templated methods for Vertex.
#-------------------------------------------------------------------------------
@PYB11ignore
def addVertexMethods(cls, ndim):

    # Constructors
    def pyinit0(self):
        "Default constructor"
        return

    def pyinit1(rhs = "const PolyClipper::Vertex%id" % ndim):
        "Copy constructor"
        return

    def pyinit2(position = "const PolyClipper::Plane%id::Vector&" % ndim):
        "Construct with a position."
        return

    def pyinit3(position = "const PolyClipper::Plane%id::Vector&" % ndim,
                c = "int"):
        "Construct with a position and initial compare flag."
        return

    # Attributes
    @PYB11readwrite
    def position(self):
        "The position of the vertex."
        return "None"

    @PYB11readwrite
    def neighbors(self):
        "The connectivty of the vertex."
        return "None"

    @PYB11readwrite
    def comp(self):
        "The current comparison flag."
        return "None"

    @PYB11readwrite
    def ID(self):
        "The ID or index of the vertex."
        return "None"

    def __eq__(self, other="py::self"):
        return "None"

    for x in [x for x in dir() if type(eval(x)) == types.FunctionType]: 
        exec("cls.%s = %s" % (x, x))

#-------------------------------------------------------------------------------
# Plane2d
#-------------------------------------------------------------------------------
@PYB11cppname("PolyClipper::Plane2d")
class PolyClipperPlane2d:
    """Plane class for polyclipper in 2 dimensions."""

addPlaneMethods(PolyClipperPlane2d, 2)

#-------------------------------------------------------------------------------
# Plane3d
#-------------------------------------------------------------------------------
@PYB11cppname("PolyClipper::Plane3d")
class PolyClipperPlane3d:
    """Plane class for polyclipper in 3 dimensions."""

addPlaneMethods(PolyClipperPlane3d, 3)

#-------------------------------------------------------------------------------
# Vertex2d
#-------------------------------------------------------------------------------
@PYB11cppname("PolyClipper::Vertex2d")
class Vertex2d:
    """Vertex class for polyclipper in 2 dimensions."""

addVertexMethods(Vertex2d, 2)

#-------------------------------------------------------------------------------
# Vertex3d
#-------------------------------------------------------------------------------
@PYB11cppname("PolyClipper::Vertex3d")
class Vertex3d:
    """Vertex class for polyclipper in 3 dimensions."""

addVertexMethods(Vertex3d, 3)

