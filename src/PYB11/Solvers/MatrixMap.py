from PYB11Generator import *

@PYB11holder("std::shared_ptr")
class MatrixMap:
    def pyinit(self):
        "Contains the indexing information for local nodes"

    @PYB11pure_virtual
    @PYB11const
    def firstGlobalIndex(self):
        "Global index of first node"
        return "int"

    @PYB11pure_virtual
    @PYB11const
    def lastGlobalIndex(self):
        "Global index of last node"
        return "int"

    @PYB11pure_virtual
    @PYB11const
    def numLocalElements(self):
        "Number of elements on this processor"
        return "int"

    @PYB11pure_virtual
    @PYB11const
    def numGlobalElements(self):
        "Number of global elements"
        return "int"

    @PYB11pure_virtual
    @PYB11const
    def numElementsPerRow(self,
                          localIndex = "int"):
        "Get the number of nonzero columns per row"
        return "int"
    
    @PYB11virtual
    @PYB11const
    @PYB11pycppname("numElementsPerRow")
    def numElementsPerRow2(self):
        "Get the number of nonzero columns per row"
        return "std::vector<int>"
    
    @PYB11pure_virtual
    @PYB11const
    def getGlobalIndex(self, localIndex = "int"):
        "Get global from local index"
        return "int"
    
    # @PYB11pure_virtual
    # @PYB11const
    # def getColumnIndices(localRowIndex = "int",
    #                      globalColumnIndices = "std::vector<int>&"):
    #     "Get the nonzero column indices for a given row"
    #     return "void"
        
    @PYB11pure_virtual
    @PYB11const
    def checkClassInvariants(self):
        "Check class invariants"
        return "void"
