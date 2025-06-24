from PYB11Generator import *
from MatrixMap import MatrixMap

@PYB11holder("std::shared_ptr")
class SimpleMatrixMap(MatrixMap):
    def pyinit(self,
               numLocalElements = "int"):
        "Form a simple linear system"
        
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
    def getGlobalIndex(self, localIndex = "int"):
        "Get global from local index"
        return "int"
    
    @PYB11virtual
    @PYB11const
    def checkClassInvariants(self):
        "Check class invariants"
        return "void"
