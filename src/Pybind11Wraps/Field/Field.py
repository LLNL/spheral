from PYB11Generator import *

#-------------------------------------------------------------------------------
# Field
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "Value")
class Field:

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

    # def pyinit4(self, rhs="const FieldType&"):
    #     "Copy constructor"

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
