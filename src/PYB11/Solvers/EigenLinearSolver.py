from PYB11Generator import *
from LinearSolver import LinearSolver

@PYB11holder("std::shared_ptr")
class EigenLinearSolver(LinearSolver):
    def pyinit(self,
               options = "std::shared_ptr<EigenOptions>"):
        "Eigen solver"

    @PYB11virtual
    @PYB11const
    def readyToSolve(self):
        "Is the system ready to solve?"
        return "bool"
    
    @PYB11virtual
    def solve(self,
              input = "const std::vector<double>&",
              output = "std::vector<double>&"):
        "Solve the equations with the current data and map"
        return "void"

    @PYB11virtual
    def multiply(self,
                 input = "const std::vector<double>&",
                 output = "std::vector<double>&"):
        "Multiply the matrix by a vector"
        return "void"

    @PYB11virtual
    @PYB11const
    def statistics(self):
        return "std::vector<std::shared_ptr<IncrementalStatistic<double>>>"
    
    @PYB11protected
    @PYB11virtual
    def initializeGraph(self):
        return "void"
    
    @PYB11protected
    @PYB11virtual
    def initializeMatrix(self):
        return "void"
    
    @PYB11protected
    @PYB11virtual
    def initializeSolver(self):
        return "void"
