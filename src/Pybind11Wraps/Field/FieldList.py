from PYB11Generator import *

#-------------------------------------------------------------------------------
# FieldList
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "Value")
@PYB11module("SpheralField")
class FieldList:

    typedefs="""
    typedef FieldList<%(Dimension)s, %(Value)s> FieldListType;
    typedef Field<%(Dimension)s, %(Value)s> FieldType;
    typedef NodeList<%(Dimension)s> NodeListType;
    typedef %(Dimension)s::Vector Vector;
    typedef %(Dimension)s::SymTensor SymTensor;
"""

    def pyinit(self):
        "Default constructor"

    def pyinit1(self, name="FieldStorageType"):
        "Build with the storage type for the fields"

    def pyinitCopy(self, rhs="const FieldListType&"):
        "Copy constructor"

    #...........................................................................
    # Methods
    def storageType(self):
        "The method whereby Fields are stored/referenced by this FieldList"
        return

    def copyFields(self):
        """Force the the FieldList to make copies of all the Fields it currently\\npoints to and point to the internally cached copies instead.\\nAlso sets storageType to CopyFields."""

    @PYB11const
    def haveField(self, field="const FieldType&"):
        "Test if the specified field is pointed to by this FieldList"
        return "bool"

    @PYB11const
    def haveNodeList(self, nodeList="const NodeListType&"):
        "Test if there are any Fields based on the given NodeList currently pointed\\nto by this FieldList"
        return "bool"

    def assignFields(self, fieldList="const FieldListType&"):
        "Set all Fields pointed to by this FieldList equal to those of the given FieldList."
        return "void"

    def appendField(self, field="const FieldType&"):
        "Add the given Field to this FieldList.\\nstorageType==ReferenceFields -> just make a new member Field pointing to the given Field.\\nstorageType==CopyFields -> copy the Field into a locally cached\\nversion, and add a\\n                           reference pointing to that local copy"
        return "void"

    def deleteField(self, field="const FieldType&"):
        "Remove the given Field from this FieldList."
        return "void"

    def appendNewField(self,
                       name = "const std::string",
                       nodeList = "const NodeListType&",
                       value = "const %(Value)s"):
        "Create a new field with the given (name, value) and add it to the FieldList.\\nNote this only makes sense when storageType==CopyFields"
        return "void"

    @PYB11returnpolicy("reference_internal")
    @PYB11implementation("[](const FieldListType& self, const NodeListType& nodeList) { return *self.fieldForNodeList(nodeList); }")
    def fieldForNodeList(self, nodeList="const NodeListType&"):
        "Return the Field associated with the given NodeList"
        return "FieldType*"

    @PYB11const
    def setMasterNodeLists(self,
                           r = "const Vector&",
                           H = "const SymTensor&",
                           masterLists = "std::vector<std::vector<int>>&",
                           coarseNeighbors = "std::vector<std::vector<int>>&"):
        "Set the master neighbor information based on (r,H)"
        return "void"

    @PYB11const
    def setMasterNodeLists(self,
                           r = "const Vector&",
                           masterLists = "std::vector<std::vector<int>>&",
                           coarseNeighbors = "std::vector<std::vector<int>>&"):
        "Set the master neighbor information based on position r (assume zero associated length scale)"
        return "void"

    @PYB11const
    def setRefineNodeLists(self,
                           r = "const Vector&",
                           H = "const SymTensor&",
                           coarseNeighbors = "const std::vector<std::vector<int>>&",
                           refineNeighbors = "std::vector<std::vector<int>>&"):
        "Set the refine neighbor information based on (r,H)"
        return "void"

    @PYB11const
    def setRefineNodeLists(self,
                           r = "const Vector&",
                           coarseNeighbors = "const std::vector<std::vector<int>>&",
                           refineNeighbors = "std::vector<std::vector<int>>&"):
        "Set the refine neighbor information based on r (assume zero associated length scale)"
        return "void"

    def Zero(self):
        "Set all element values equal to zero"
        return "void"

    @PYB11const
    def size(self):
        "Number of Fields"
        return "unsigned"

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def nodeListPtrs(self):
        "The NodeLists for Fields in this FieldList"
        return "const std::vector<NodeListType*>&"

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
    def __getitem__(self, index="const unsigned"):
        return "FieldType*"

    @PYB11implementation("[](const FieldListType& self) { return py::make_iterator(self.begin(), self.end()); }")
    def __iter__(self):
        "Python iteration through a FieldList."

    def __call__(self,
                 fieldIndex = "const unsigned",
                 nodeIndex = "const unsigned"):
        "Return the %(Value)s for the given (fieldIndex, nodeIndex)."
        return "%(Value)s&"

    #...........................................................................
    # Properties (and methods for em)
    @PYB11cppname("storageType")
    @PYB11const
    def getstorageType(self):
        return "FieldStorageType"

    @PYB11cppname("numFields")
    @PYB11const
    def getnumFields(self):
        return "unsigned"

    @PYB11cppname("numNodes")
    @PYB11const
    def getnumNodes(self):
        return "unsigned"

    @PYB11cppname("numInternalNodes")
    @PYB11const
    def getnumInternalNodes(self):
        return "unsigned"

    @PYB11cppname("numGhostNodes")
    @PYB11const
    def getnumGhostNodes(self):
        return "unsigned"

    storageType = property(getstorageType, doc="The method whereby Fields are stored/referenced by this FieldList")
    numFields = property(getnumFields, doc="Number of Fields")
    numNodes = property(getnumNodes, doc="Number of nodes in all the associated Fields")
    numInternalNodes = property(getnumInternalNodes, doc="Number of internal nodes in all the associated Fields")
    numGhostNodes = property(getnumGhostNodes, doc="Number of ghost nodes in all the associated Fields")

#-------------------------------------------------------------------------------
# FieldList with numeric operations
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "Value")
@PYB11pycppname("FieldList")
class ArithmeticFieldList(FieldList):

    # @PYB11const
    # @PYB11cppname("operator()")
    # def valueat(self,
    #             position = "const %(Dimension)s::Vector&",
    #             W = "const TableKernel<%(Dimension)s>&"):
    #     "Return the interpolated value of the FieldList at a position."
    #     return "%(Value)s"

    def __add__(self):
        return

    def __sub__(self):
        return

    def __iadd__(self):
        return

    def __isub__(self):
        return

    def __add__(self, rhs="%(Value)s()"):
        return "FieldType"

    def __sub__(self, rhs="%(Value)s()"):
        return "FieldType"

    def __iadd__(self, rhs="%(Value)s()"):
        return

    def __isub__(self, rhs="%(Value)s()"):
        return

    def __imul__(self, rhs="double()"):
        return

    def __idiv__(self, rhs="double()"):
        return

    @PYB11const
    def sumElements(self):
        "Return the sum of the elements in the Field."
        return

    @PYB11const
    def localSumElements(self):
        "Return the sum of the elements in the Field local to each processor."
        return

    #...........................................................................
    # Comparators
    def __gt__(self):
        return

    def __lt__(self):
        return

    def __ge__(self):
        return "bool"

    def __le__(self):
        return "bool"

    def __gt__(self, rhs="%(Value)s()"):
        "Greater than comparision with a %(Value)s"
        return "bool"

    def __lt__(self, rhs="%(Value)s()"):
        "Less than comparision with a %(Value)s"
        return "bool"

    def __ge__(self, rhs="%(Value)s()"):
        "Greater than or equal comparision with a %(Value)s"
        return "bool"

    def __le__(self, rhs="%(Value)s()"):
        "Less than or equal comparision with a %(Value)s"
        return "bool"

    def applyMin(self):
        "Enforce a floor on the values of the Field."
        return

    def applyMax(self):
        "Enforce a ceiling on the values of the Field."
        return

    @PYB11const
    def min(self):
        "Return the mimimum value in the Field."
        return

    @PYB11const
    def max(self):
        "Return the maximum value in the Field."
        return

    @PYB11const
    def localMin(self):
        "Return the mimimum value in the Field local to each processor."
        return

    @PYB11const
    def localMax(self):
        "Return the maximum value in the Field local to each processor."
        return

#-------------------------------------------------------------------------------
# Add min/max operations to a Field
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "Value")
@PYB11pycppname("FieldList")
class MinMaxFieldList(ArithmeticFieldList):

    def applyScalarMin(self):
        "Enforce a float floor on the values of the Field."
        return

    def applyScalarMax(self):
        "Enforce a float ceiling on the values of the Field."
        return
