"""
PolyClipper module.

Binds the PolyClipper geometry operations for clipping polygons & polyhedra.
"""

from PYB11Decorators import *

# Include files.
includes = ['"Geometry/polyclipper.hh"']

#-------------------------------------------------------------------------------
# Plane2d
#-------------------------------------------------------------------------------
@PYB11cppname("PolyClipper::Plane2d")
class PolyClipperPlane2d:
    """Plane class for polyclipper in 2 dimensions."""
    ndim = 2    

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
