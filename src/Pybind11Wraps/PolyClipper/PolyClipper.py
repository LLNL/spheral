"""
PolyClipper module.

Binds the PolyClipper geometry operations for clipping polygons & polyhedra.
"""

from PYB11Decorators import *

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

    for x in ("pyinit0", "pyinit1", "pyinit2", "pyinit3", "dist", "normal"):
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
