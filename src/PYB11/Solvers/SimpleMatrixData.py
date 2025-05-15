from PYB11Generator import *
from MatrixData import MatrixData

@PYB11holder("std::shared_ptr")
class SimpleMatrixData(MatrixData):
    def pyinit(self,
               map = "std::shared_ptr<SimpleMatrixMap>"):
        "Form a simple linear system"

    @PYB11virtual
    def setRowValues(localRowIndex = "int",
                     globalColumnIndices = "std::vector<int> const&",
                     globalColumnValues = "std::vector<double> const&"):
        "Set value of a row"
        return "void"
        
    @PYB11virtual
    @PYB11const
    def getMap(self):
        "Get map for this matrix"
        return "std::shared_ptr<MatrixMap>"

    @PYB11virtual
    @PYB11const
    def getRowValues(self,
                     localRowIndex = "int",
                     globalColumnIndices = "std::vector<int>&",
                     globalColumnValues = "std::vector<double>&"):
        "Return the nonzero indices and corresponding values for this row"
        return "void"
    
    @PYB11virtual
    @PYB11const
    def checkClassInvariants(self):
        "Check class invariants"
        return "void"
