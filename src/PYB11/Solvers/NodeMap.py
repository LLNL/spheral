from PYB11Generator import *
from MatrixMap import MatrixMap

@PYB11holder("std::shared_ptr")
@PYB11template("Dimension")
class NodeMap(MatrixMap):
    def pyinit(self,
               flatConnectivity = "FlatConnectivity<%(Dimension)s>&"):
        "Contains the indexing information for local nodes"

    @PYB11virtual
    @PYB11const
    def firstGlobalIndex(self):
        "Global index of first node"
        return "int"

    @PYB11virtual
    @PYB11const
    def lastGlobalIndex(self):
        "Global index of last node"
        return "int"

    @PYB11virtual
    @PYB11const
    def numLocalElements(self):
        "Number of elements on this processor"
        return "int"

    @PYB11virtual
    @PYB11const
    def numGlobalElements(self):
        "Number of global elements"
        return "int"

    @PYB11virtual
    @PYB11const
    def numElementsPerRow(self,
                          localIndex = "int"):
        "Get the number of nonzero columns per row"
        return "int"
    
    @PYB11virtual
    @PYB11const
    def checkClassInvariants(self):
        "Check class invariants"
        return "void"

    @PYB11virtual
    @PYB11const
    def getColumnIndices(self,
                         localRowIndex = "int",
                         globalColumnIndices = "std::vector<int>&"):
        "Get the nonzero column indices for a given row"
        return "void"
    
    @PYB11virtual
    @PYB11const
    def getLocalIndex(self,
                      nodeListIndex = "int",
                      nodeIndex = "int"):
        "Convert node indices to matrix index"
        return "int"
     
    @PYB11virtual
    @PYB11const
    def getNodeIndex(self,
                     localIndex = "int"):
        "Convert matrix index to node indices"
        return "std::pair<int, int>"

    @PYB11virtual
    @PYB11const
    def getGlobalIndex(self,
                       nodeListIndex = "int",
                       nodeIndex = "int"):
        "Convert node indices to global matrix index"
        return "int"

