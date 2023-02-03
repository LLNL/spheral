import inspect
from PYB11Generator import *
from FieldBase import FieldBase

#-------------------------------------------------------------------------------
# Field
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "Value")
@PYB11module("SpheralField")
class Field(FieldBase):

    PYB11typedefs = """
  typedef Field<%(Dimension)s, %(Value)s> FieldType;
"""

    def pyinit(self, name="std::string"):
        "Construct with a name"

    def pyinit1(self, name="std::string", nodeList="const FieldType&"):
        "Construct as a copy of a Field with a new name"

    def pyinit2(self, name="std::string", nodeList="const NodeList<%(Dimension)s>&"):
        "Construct with a name and NodeList"

    def pyinit3(self,
                name = "std::string",
                nodeList = "const NodeList<%(Dimension)s>&",
                val = "%(Value)s"):
        "Construct with a name, NodeList, and initial value"

    # def pyinit4(self, rhs="const PYB11TrampolineField<%(Dimension)s, %(Value)s>&"):
    #     "Copy constructor"

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
        return "unsigned"

    @PYB11cppname("operator[]")
    @PYB11returnpolicy("reference_internal")
    @PYB11implementation('[](FieldType& self, int i) { const int n = self.size(); if (i >= n) throw py::index_error(); return &self[(i %% n + n) %% n]; }')
    def __getitem__(self):
        #return "%(Value)s&"
        return

    @PYB11implementation("[](FieldType& self, int i, const %(Value)s v) { const int n = self.size(); if (i >= n) throw py::index_error(); self[(i %% n + n) %% n] = v; }")
    def __setitem__(self):
        "Set a value"

    @PYB11implementation("[](const FieldType& self) { return py::make_iterator(self.begin(), self.end()); }, py::keep_alive<0,1>()")
    def __iter__(self):
        "Python iteration through a Field."

    @PYB11returnpolicy("reference_internal")
    @PYB11implementation("[](FieldType& self, int i) { const int n = self.size(); if (i >= n) throw py::index_error(); return &self[(i %% n + n) %% n]; }")
    def __call__(self):
        "Index into a Field"
        #return "%(Value)s&"
        return

    #...........................................................................
    # FieldBase virtual methods
    @PYB11virtual
    @PYB11const
    def size(self):
        "Number of elements"
        return "unsigned"

    @PYB11virtual
    def Zero(self):
        "Set all element values equal to zero"
        return "void"

    @PYB11virtual
    def setNodeList(self, nodeList="const NodeList<%(Dimension)s>&"):
        "Register this Field with the given NodeList"
        return "void"

    @PYB11virtual
    def resizeField(self, size="unsigned"):
        "Set the number of elements"
        return "void"

    @PYB11virtual
    def resizeFieldInternal(self, size="unsigned", oldFirstGhostNode="unsigned"):
        "Set the number of internal elements"
        return "void"

    @PYB11virtual
    def resizeFieldGhost(self, size="unsigned"):
        "Set the number of ghost elements"
        return "void"

    @PYB11virtual
    def deleteElement(self, nodeID="int"):
        "Delete the element at the given index"
        return "void"

    @PYB11virtual
    def deleteElements(self, nodeIDs="const std::vector<int>&"):
        "Delete the elements at the given indices"
        return "void"

    @PYB11virtual
    @PYB11const
    def packValues(self, nodeIDs="const std::vector<int>&"):
        "Serialize the indicated elements into a vector<char>"
        return "std::vector<char>"

    @PYB11virtual
    def unpackValues(self,
                     nodeIDs="const std::vector<int>&",
                     buffer = "const std::vector<char>&"):
        "Deserialize values from the given buffer"
        return "void"

    @PYB11virtual
    def copyElements(self,
                     fromIndices="const std::vector<int>&",
                     toIndices="const std::vector<int>&"):
        "Copy a range of values from/to elements of the Field"
        return "void"

    #...........................................................................
    # Methods
    @PYB11const
    @PYB11implementation("[](const Field<%(Dimension)s, %(Value)s>& self, const int precision) -> py::bytes { return py::bytes(self.string(precision)); }")
    def string(self, precision=("const int", "20")):
        "Serialize Field to a string"
        return "py::bytes"

    @PYB11pycppname("string")
    @PYB11implementation("[](Field<%(Dimension)s, %(Value)s>& self, const py::bytes& buf) -> void { self.string(static_cast<std::string>(buf)); }")
    def string1(self, buf="py::bytes&"):
        "Deserialize from a string"
        return "void"

    @PYB11const
    @PYB11implementation("[](const Field<%(Dimension)s, %(Value)s>& self) -> py::list { const auto vals = self.internalValues(); py::list result; for (const auto& x: vals) result.append(x); return result; }")
    def internalValues(self):
        "Return a python list (as a copy) of just the internal values"
        return "py::list"

    @PYB11const
    @PYB11implementation("[](const Field<%(Dimension)s, %(Value)s>& self) -> py::list { const auto vals = self.ghostValues(); py::list result; for (const auto& x: vals) result.append(x); return result; }")
    def ghostValues(self):
        "Return a python list (as a copy) of just the ghost values"
        return "py::list"

    @PYB11const
    @PYB11implementation("[](const Field<%(Dimension)s, %(Value)s>& self) -> py::list { const auto vals = self.allValues(); py::list result; for (const auto& x: vals) result.append(x); return result; }")
    def allValues(self):
        "Return a python list (as a copy) of all values in the Field"
        return "py::list"

    #...........................................................................
    # Properties
    numElements = PYB11property("unsigned", doc="Number of elements in field")
    numInternalElements = PYB11property("unsigned", doc="Number of internal elements in field")
    numGhostElements = PYB11property("unsigned", doc="Number of ghost elements in field")
