#-------------------------------------------------------------------------------
# GridCellIndex template
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11template("Dimension")
class GridCellIndex:

    PYB11typedefs = """
    typedef GridCellPlane<%(Dimension)s> GridCellPlaneType;
"""

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

    def dot(self):
        "Dot product with another GridCellIndex"
        return

    def compare(self):
        "Compare with another gridcell (-1,0,1)"
        return

    def inRange(self):
        "Test if the gridcell is inside the range specified by the arguments"

    def magnitude(self):
        "The magnitude of the Vector"

    def magnitude2(self):
        "The square of the magnitude of the Vector"

    def minElement(self):
        "Minimum int coordinate in Vector"

    def maxElement(self):
        "Maximum int coordinate in Vector"

    def sumElements(self):
        "Sum of the coordinate"

    def productElements(self):
        "Product of the coordinates"

    def indexMin(self):
        "Return the minimum allowed coordinate"

    def indexMax(self):
        "Return the maximum allowed coordinate"

    #...........................................................................
    # Sequence methods
    @PYB11implementation("[](const GridCellIndex<%(Dimension)s>&) { return %(Dimension)s::nDim; }")
    def __len__(self):
        "The size (in number of coordinates) of the GridCellIndex."

    @PYB11implementation("[](const GridCellIndex<%(Dimension)s> &s, int i) { if (i >= %(Dimension)s::nDim) throw py::index_error(); return s(i); }") 
    @PYB11returnpolicy("reference_internal")
    def __getitem__(self):
        "Python indexing to get a coordinate."

    @PYB11implementation("[](GridCellIndex<%(Dimension)s> &s, int i, double v) { if (i >= %(Dimension)s::nDim) throw py::index_error(); s(i) = v; }") 
    def __setitem__(self):
        "Python indexing to set a coordinate."

    @PYB11implementation("[](const GridCellIndex<%(Dimension)s> &s) { return py::make_iterator(std::begin(s), std::end(s)); }, py::keep_alive<0,1>()")
    def __iter__(self):
        "Python iteration through a GridCellIndex."

    @PYB11const
    def __call__(self, i="int"):
        "Index for a coordinate using parens."
        return "int"

    #...........................................................................
    # Operators
    def __neg__(self):
        return
    def __add__(self):
        return
    def __sub__(self):
        return
    def __iadd__(self):
        return
    def __isub__(self):
        return

    def __add__(self, rhs="int()"):
        return
    def __sub__(self, rhs="int()"):
        return
    def __mul__(self, rhs="int()"):
        return
    def __rmul__(self, rhs="int()"):
        return
    def __truediv__(self, rhs="int()"):
        return
    def __add__(self, rhs="int()"):
        return
    def __isub__(self, rhs="int()"):
        return

    def __eq__(self):
        return
    def __ne__(self):
        return
    def __lt__(self):
        return
    def __gt__(self):
        return
    def __le__(self):
        return
    def __ge__(self):
        return

    def __lt__(self, rhs="GridCellPlaneType()"):
        return
    def __gt__(self, rhs="GridCellPlaneType()"):
        return
    def __le__(self, rhs="GridCellPlaneType()"):
        return
    def __ge__(self, rhs="GridCellPlaneType()"):
        return
