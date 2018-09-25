#-------------------------------------------------------------------------------
# GridCellPlane template
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11template("Dimension")
class GridCellPlane:

    typedefs="""
    typedef GridCellIndex<%(Dimension)s> GridCellIndexType;
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
    @PYB11pycppname("point")
    @PYB11ignore
    @PYB11const
    def getPoint(self):
        return "GridCellIndexType"

    @PYB11ignore
    def setPoint(self, val="GridCellIndexType"):
        return "void"

    @PYB11pycppname("normal")
    @PYB11ignore
    @PYB11const
    def getNormal(self):
        return "GridCellIndexType"

    @PYB11ignore
    def setNormal(self, val="GridCellIndexType"):
        return "void"

    point = property(getPoint, setPoint, "Point in the plane")
    normal = property(getNormal, setNormal, "Normal to plane")
