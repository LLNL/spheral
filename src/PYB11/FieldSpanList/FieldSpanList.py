from PYB11Generator import *
#import FieldSpanListBase

#-------------------------------------------------------------------------------
# FieldSpanList
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "Value")
@PYB11module("SpheralFieldSpanList")
class FieldSpanList:

    PYB11typedefs = """
    using FieldListType = FieldList<%(Dimension)s, %(Value)s>;
    using FieldSpanListType = FieldSpanList<%(Dimension)s, %(Value)s>;
    using FieldSpanType = FieldSpan<%(Dimension)s, %(Value)s>;
    using Scalar = %(Dimension)s::Scalar;
"""

    def pyinitCopy(self, rhs="FieldSpanListType&"):
        "Copy constructor"

    #...........................................................................
    # Methods
    def assignFields(self, rhs="FieldSpanListType&"):
        "Set all Fields pointed to by this FieldSpanList equal to those of the given FieldSpanList."
        return "void"

    def Zero(self):
        "Set all element values equal to zero"
        return "void"

    @PYB11const
    def size(self):
        "Number of Fields"
        return "size_t"

    #...........................................................................
    # Comparators
    def __eq__(self):
        return

    def __ne__(self):
        return

    def __eq__(self, rhs="%(Value)s()"):
        "Equivalence comparision with a %(Value)s"
        return "bool"

    def __ne__(self, rhs="%(Value)s()"):
        "Not equal comparision with a %(Value)s"
        return "bool"

    #...........................................................................
    # Sequence methods
    @PYB11cppname("size")
    @PYB11const
    def __len__(self):
        return "size_t"

    @PYB11cppname("operator[]")
    @PYB11returnpolicy("reference")
    @PYB11keepalive(0,1)
    def __getitem__(self, index="size_t"):
        return "FieldSpanType&"

    @PYB11returnpolicy("reference")
    @PYB11implementation("[](FieldSpanListType& self) { return py::make_iterator(self.begin(), self.end()); }, py::keep_alive<0,1>()")
    def __iter__(self):
        "Python iteration through a FieldSpanList."

    def __call__(self,
                 fieldIndex = "size_t",
                 nodeIndex = "size_t"):
        "Return the %(Value)s for the given (fieldIndex, nodeIndex)."
        return "%(Value)s&"

    #...........................................................................
    # Properties
    numFields = PYB11property("size_t", doc="Number of Fields")
    numElements = PYB11property("size_t", doc="Number of elements in all the associated Fields")
    numInternalElements = PYB11property("size_t", doc="Number of internal elements in all the associated Fields")
    numGhostElements = PYB11property("size_t", doc="Number of ghost elementes in all the associated Fields")
