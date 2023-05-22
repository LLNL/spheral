#-------------------------------------------------------------------------------
# GridCellPlane template
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11template("Dimension")
class GridCellPlane:

    PYB11typedefs = """
    typedef GridCellIndex<%(Dimension)s> GridCellIndexType;
    //typedef GridCellPlane<%(Dimension)s> GridCellPlaneType;
"""

    #...........................................................................
    # Constructors
    def pyinit0(self):
        "Default constructor"

    def pyinit1(self, point="const GridCellIndexType&", normal="const GridCellIndexType&"):
        "Construct with (point, normal)"

    def pyinit4(self, rhs="const GridCellPlane<%(Dimension)s>&"):
        "Copy constructor"

    def valid(self):
        "Test if the Plane is valid, ready to use"

    #...........................................................................
    # Methods
    def minimumDistance(self):
        "Distance from GridCellIndex to plane"
        return

    def coplanar(self):
        "Test if GridCellIndex is in the plane"
        return

    def parallel(self):
        "Test if the plane is parallel to this one"

    #...........................................................................
    # Operators
    def __eq__(self):
        return
    def __ne__(self):
        return

    def __lt__(self, rhs="GridCellIndexType()"):
        return
    def __gt__(self, rhs="GridCellIndexType()"):
        return
    def __le__(self, rhs="GridCellIndexType()"):
        return
    def __ge__(self, rhs="GridCellIndexType()"):
        return

    #...........................................................................
    # Properties
    point = PYB11property("", "point", "point", "Point in the plane")
    normal = PYB11property("", "normal", "normal", "Normal to plane")
