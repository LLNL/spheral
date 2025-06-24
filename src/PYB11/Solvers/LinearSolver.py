from PYB11Generator import *

@PYB11holder("std::shared_ptr")
class LinearSolver:
    def pyinit(self):
        "Pure virtual linear solver class"

    @PYB11virtual
    @PYB11const
    def mapSet(self):
        "Is the map set?"
        return "bool"

    @PYB11virtual
    @PYB11const
    def dataSet(self):
        "Is the data set?"
        return "bool"
    
    @PYB11pure_virtual
    @PYB11const
    def readyToSolve(self):
        "Is the system ready to solve?"
        return "bool"

    @PYB11pure_virtual
    def solve(self,
              input = "const std::vector<double>&",
              output = "std::vector<double>&"):
        "Solve the equations with the current data and map"
        return "void"

    @PYB11pure_virtual
    def multiply(self,
                 input = "const std::vector<double>&",
                 output = "std::vector<double>&"):
        "Multiply the matrix by a vector"
        return "void"
    
    map = PYB11property(returnType="std::shared_ptr<MatrixMap>",
                        getter="map", setter="setMap",
                        getterconst=True,
                        doc="The map, which also computes stuff when set")
    data = PYB11property(returnType="std::shared_ptr<MatrixData>",
                         getter="data", setter="setData",
                         getterconst=True,
                         doc="The data, which also computes stuff when set")

    @PYB11virtual
    @PYB11const
    def statistics(self):
        return "std::vector<std::shared_ptr<IncrementalStatistic<double>>>"
    
    @PYB11protected
    @PYB11pure_virtual
    def initializeGraph(self):
        return "void"
    
    @PYB11protected
    @PYB11pure_virtual
    def initializeMatrix(self):
        return "void"
    
    @PYB11protected
    @PYB11pure_virtual
    def initializeSolver(self):
        return "void"
    
