"""
PolyClipper module.

Binds the PolyClipper geometry operations for clipping polygons & polyhedra.
"""

from PYB11Decorators import *

# Include files.
includes = ['"Geometry/polyclipper.hh"']

#-------------------------------------------------------------------------------
# Plane
#-------------------------------------------------------------------------------
@PYB11cppname("PolyClipper::Plane2d")
class PolyClipperPlane2d:
    "Plane class for polyclipper in %i dimensions." % 2

    # Constructors
    def pyinit0(self):
        "Default constructor"
        return

    def pyinit1(rhs = "const PolyClipper::Plane%id" % 2):
        "Copy constructor"
        return
