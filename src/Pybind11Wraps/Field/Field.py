import inspect
from PYB11Generator import *
from FieldBase import FieldBase

#-------------------------------------------------------------------------------
# Field
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "Value")
@PYB11module("SpheralField")
class Field(FieldBase):

    typedefs="""
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
    def __getitem__(self, index="unsigned"):
        return "%(Value)s&"

    @PYB11implementation("[](FieldType& self, size_t i, const %(Value)s v) { if (i >= self.size()) throw py::index_error(); self[i] = v; }") 
    def __setitem__(self):
        "Set a value"

    @PYB11implementation("[](const FieldType& self) { return py::make_iterator(self.begin(), self.end()); }, py::keep_alive<0,1>()")
    def __iter__(self):
        "Python iteration through a Field."

    @PYB11returnpolicy("reference_internal")
    def __call__(self, i="int"):
        "Index into a Field"
        return "%(Value)s&"

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
                     numElements = "int",
                     beginInsertionIndex = "int",
                     buffer = "const std::vector<char>&"):
        "Deserialize values from the given buffer"
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
    def internalValues(self):
        "Return a vector<%(Value)s> of just the internal values"
        return "std::vector<%(Value)s>"

    @PYB11const
    def ghostValues(self):
        "Return a vector<%(Value)s> of just the ghost values"
        return "std::vector<%(Value)s>"

    @PYB11const
    def allValues(self):
        "Return a vector<%(Value)s> of all values"
        return "std::vector<%(Value)s>"

    #...........................................................................
    # operators

    #...........................................................................
    # Properties
    numElements = PYB11property("unsigned", doc="Number of elements in field")
    numInternalElements = PYB11property("unsigned", doc="Number of elements in field")

