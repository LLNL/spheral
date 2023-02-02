from PYB11Generator import *

#-------------------------------------------------------------------------------
# FieldView
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "Value")
@PYB11module("SpheralField")
class FieldView():
    "FieldView class"
    PYB11typedefs = """
  typedef Field<%(Dimension)s, %(Value)s> FieldType;
  typedef FieldView<%(Dimension)s, %(Value)s> FieldViewType;
"""


    #...........................................................................
    # Sequence methods
    #@PYB11cppname("size")
    @PYB11const
    @PYB11implementation('[](FieldViewType& self) { return self.get().size(); }')
    def __len__(self):
        return "unsigned"

    @PYB11cppname("operator[]")
    @PYB11returnpolicy("reference_internal")
    @PYB11implementation('[](FieldViewType& self, int i) { const int n = self.get().size(); if (i >= n) throw py::index_error(); return &self[(i %% n + n) %% n]; }')
    def __getitem__(self):
        #return "%(Value)s&"
        return

    @PYB11implementation("[](FieldViewType& self, int i, const %(Value)s v) { const int n = self.get().size(); if (i >= n) throw py::index_error(); self[(i %% n + n) %% n] = v; }")
    def __setitem__(self):
        "Set a value"

    @PYB11implementation("[](const FieldViewType& self) { return py::make_iterator(self.get().begin(), self.get().end()); }, py::keep_alive<0,1>()")
    def __iter__(self):
        "Python iteration through a Field."

    @PYB11returnpolicy("reference_internal")
    @PYB11implementation("[](FieldViewType& self, int i) { const int n = self.get().size(); if (i >= n) throw py::index_error(); return &self[(i %% n + n) %% n]; }")
    def __call__(self):
        "Index into a Field"
        #return "%(Value)s&"
        return

#    #...........................................................................
#    # FieldBase virtual methods
#    @PYB11const
#    def size(self):
#        "Number of elements"
#        return "unsigned"
#
#    def Zero(self):
#        "Set all element values equal to zero"
#        return "void"
#
#    def setNodeList(self, nodeList="const NodeList<%(Dimension)s>&"):
#        "Register this Field with the given NodeList"
#        return "void"
#
#    def resizeField(self, size="unsigned"):
#        "Set the number of elements"
#        return "void"
#
#    def resizeFieldInternal(self, size="unsigned", oldFirstGhostNode="unsigned"):
#        "Set the number of internal elements"
#        return "void"
#
#    def resizeFieldGhost(self, size="unsigned"):
#        "Set the number of ghost elements"
#        return "void"
#
#    def deleteElement(self, nodeID="int"):
#        "Delete the element at the given index"
#        return "void"
#
#    def deleteElements(self, nodeIDs="const std::vector<int>&"):
#        "Delete the elements at the given indices"
#        return "void"
#
#    @PYB11const
#    def packValues(self, nodeIDs="const std::vector<int>&"):
#        "Serialize the indicated elements into a vector<char>"
#        return "std::vector<char>"
#
#    def unpackValues(self,
#                     nodeIDs="const std::vector<int>&",
#                     buffer = "const std::vector<char>&"):
#        "Deserialize values from the given buffer"
#        return "void"
#
#    def copyElements(self,
#                     fromIndices="const std::vector<int>&",
#                     toIndices="const std::vector<int>&"):
#        "Copy a range of values from/to elements of the Field"
#        return "void"
#
#    #...........................................................................
    # Methods
    @PYB11const
    @PYB11implementation("[](const FieldView<%(Dimension)s, %(Value)s>& self, const int precision) -> py::bytes { return py::bytes(self.get().string(precision)); }")
    def string(self, precision=("const int", "20")):
        "Serialize Field to a string"
        return "py::bytes"

    @PYB11pycppname("string")
    @PYB11implementation("[](FieldView<%(Dimension)s, %(Value)s>& self, const py::bytes& buf) -> void { self.get().string(static_cast<std::string>(buf)); }")
    def string1(self, buf="py::bytes&"):
        "Deserialize from a string"
        return "void"

    @PYB11const
    @PYB11implementation("[](const FieldViewType& self) -> py::list { const auto vals = self.get().internalValues(); py::list result; for (const auto& x: vals) result.append(x); return result; }")
    def internalValues(self):
        "Return a python list (as a copy) of just the internal values"
        return "py::list"

    @PYB11const
    @PYB11implementation("[](const FieldView<%(Dimension)s, %(Value)s>& self) -> py::list { const auto vals = self.get().ghostValues(); py::list result; for (const auto& x: vals) result.append(x); return result; }")
    def ghostValues(self):
        "Return a python list (as a copy) of just the ghost values"
        return "py::list"

    @PYB11const
    @PYB11implementation("[](const FieldView<%(Dimension)s, %(Value)s>& self) -> py::list { const auto vals = self.get().allValues(); py::list result; for (const auto& x: vals) result.append(x); return result; }")
    def allValues(self):
        "Return a python list (as a copy) of all values in the Field"
        return "py::list"

    @PYB11const
    @PYB11implementation('[](FieldViewType& self) -> py::bytes { return py::bytes(self.get().name()); }')
    def name(self):
        "Get Field Name"
        return "py::bytes"

    @PYB11implementation('[](FieldViewType& self, std::string s) { self.get().name(s); }')
    def name(self, s="std::string"):
        "Set Field Name"
        return "void"

    #name = PYB11property(getter="__getName", setter="__setName", doc="Name for the FieldView/Field.")
    #name = property(getName, setName, doc="Name for the FieldView/Field.")


#-------------------------------------------------------------------------------
# Inject base field methods
#-------------------------------------------------------------------------------
#PYB11inject(Field, FieldView)
