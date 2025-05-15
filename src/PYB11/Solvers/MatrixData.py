from PYB11Generator import *

@PYB11holder("std::shared_ptr")
class MatrixData:
    def pyinit(self):
        "Represents the data needed to create a sparse matrix"
    
    @PYB11pure_virtual
    @PYB11const
    def getMap(self):
        "Get map for this matrix"
        return "std::shared_ptr<MatrixMap>"

    @PYB11pure_virtual
    @PYB11const
    def getRowValues(self,
                     localRowIndex = "int",
                     globalColumnIndices = "std::vector<int>&",
                     globalColumnValues = "std::vector<double>&"):
        "Return the nonzero indices and corresponding values for this row"
        return "void"
    
    @PYB11virtual
    @PYB11const
    def initialGuessAvailable(self):
        "Is initial guess available?"
        return "bool"
    
    @PYB11virtual
    @PYB11const
    def getInitialGuess(self,
                        data = "std::vector<double>&"):
        "Get the initial guess"
        return "void"
    
    @PYB11pure_virtual
    @PYB11const
    def checkClassInvariants(self):
        "Check class invariants"
        return "void"
