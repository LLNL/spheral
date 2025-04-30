import inspect
from PYB11Generator import *
from FieldSpanBase import FieldSpanBase

#-------------------------------------------------------------------------------
# Field
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "Value")
@PYB11module("SpheralField")
class FieldSpan(FieldSpanBase):

    PYB11typedefs = """
  using SelfType = FieldSpan<%(Dimension)s, %(Value)s>;
  using Scalar = %(Dimension)s::Scalar;
"""

    def pyinit(self, field="Field<%(Dimension)s, %(Value)s>&"):
        "Construct from a Field"

    #...........................................................................
    # Comparators
    def __eq__(self):
        return

    def __ne__(self):
        return

    @PYB11pycppname("__eq__")
    def __eq__S__(self, rhs="%(Value)s()"):
        "Equivalence comparision with a %(Value)s"
        return "bool"

    @PYB11pycppname("__ne__")
    def __ne__S__(self, rhs="%(Value)s()"):
        "Not equal comparision with a %(Value)s"
        return "bool"

    #...........................................................................
    # Sequence methods
    @PYB11cppname("size")
    @PYB11const
    def __len__(self):
        return "size_t"

    @PYB11cppname("operator[]")
    @PYB11returnpolicy("reference_internal")
    @PYB11implementation('[](SelfType& self, int i) { const auto n = self.size(); if (size_t(i) >= n) throw py::index_error(); return &self[(i %% n + n) %% n]; }')
    def __getitem__(self):
        #return "%(Value)s&"
        return

    @PYB11implementation("[](SelfType& self, int i, const %(Value)s v) { const auto n = self.size(); if (size_t(i) >= n) throw py::index_error(); self[(i %% n + n) %% n] = v; }")
    def __setitem__(self):
        "Set a value"

    @PYB11implementation("[](SelfType& self) { return py::make_iterator(self.begin(), self.end()); }, py::keep_alive<0,1>()")
    def __iter__(self):
        "Python iteration through a FieldSpan."

    @PYB11returnpolicy("reference_internal")
    @PYB11implementation("[](SelfType& self, int i) { const auto n = self.size(); if (size_t(i) >= n) throw py::index_error(); return &self[(i %% n + n) %% n]; }")
    def __call__(self):
        "Index into a FieldSpan"
        return

    #...........................................................................
    # Methods
    @PYB11virtual
    @PYB11const
    def size(self):
        "Number of elements"
        return "size_t"

    @PYB11virtual
    def Zero(self):
        "Set all element values equal to zero"
        return "void"

    #...........................................................................
    # Properties
    numElements = PYB11property("size_t", doc="Number of elements in span")
    numInternalElements = PYB11property("size_t", doc="Number of internal elements in span")
    numGhostElements = PYB11property("size_t", doc="Number of ghost elements in span")
