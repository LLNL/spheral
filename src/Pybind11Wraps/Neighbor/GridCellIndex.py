#-------------------------------------------------------------------------------
# GridCellIndex template
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11template("Dimension")
class GridCellIndex:

    #...........................................................................
    # Constructors
    def pyinit0(self):
        "Default constructor"

    def pyinit1(self, xIndex="int"):
        "Construct with (x)"

    def pyinit2(self, xIndex="int", yIndex="int"):
        "Construct with (x,y)"

    def pyinit3(self, xIndex="int", yIndex="int", zIndex="int"):
        "Construct with (x,y,z)"

    def pyinit4(self, rhs="const GridCellIndex<%(Dimension)s>&"):
        "Copy constructor"

    #...........................................................................
    # Methods
    @PYB11pycppname("setIndices")
    def setIndices1(self, xIndex="int"):
        "Set the x index"
        return "void"

    @PYB11pycppname("setIndices")
    def setIndices2(self, xIndex="int", yIndex="int"):
        "Set the (x,y) indices"
        return "void"

    @PYB11pycppname("setIndices")
    def setIndices3(self, xIndex="int", yIndex="int", zIndex="int"):
        "Set the (x,y,z) indices"
        return "void"
