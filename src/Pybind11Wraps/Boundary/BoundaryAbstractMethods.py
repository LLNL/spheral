#-------------------------------------------------------------------------------
# Provide the abstact Boundary interface, suitable for injection to our actual
# classes.
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11ignore
class BoundaryAbstractMethods:

    def setGhostNodes(self,
                      nodeList = "NodeList<DIM>&"):
        "Use the given NodeList's neighbor object to select the ghost nodes."
        return "void"

    def updateGhostNodes(self,
                         nodeList = "NodeList<DIM>&"):
        "For the computed set of ghost nodes, set the positions and H's."
        return "void"

    @PYB11pycppname("applyGhostBoundary")
    @PYB11const
    def applyGhostBoundary1(self,
                            field = "Field<DIM, int>&"):
        "Apply the boundary condition to the ghost node values in the given Field."
        return "void"

    @PYB11pycppname("applyGhostBoundary")
    @PYB11const
    def applyGhostBoundary2(self,
                            field = "Field<DIM, Scalar>&"):
        "Apply the boundary condition to the ghost node values in the given Field."
        return "void"

    @PYB11pycppname("applyGhostBoundary")
    @PYB11const
    def applyGhostBoundary3(self,
                            field = "Field<DIM, Vector>&"):
        "Apply the boundary condition to the ghost node values in the given Field."
        return "void"

    @PYB11pycppname("applyGhostBoundary")
    @PYB11const
    def applyGhostBoundary4(self,
                            field = "Field<DIM, Tensor>&"):
        "Apply the boundary condition to the ghost node values in the given Field."
        return "void"

    @PYB11pycppname("applyGhostBoundary")
    @PYB11const
    def applyGhostBoundary5(self,
                            field = "Field<DIM, SymTensor>&"):
        "Apply the boundary condition to the ghost node values in the given Field."
        return "void"

    @PYB11pycppname("applyGhostBoundary")
    @PYB11const
    def applyGhostBoundary6(self,
                            field = "Field<DIM, ThirdRankTensor>&"):
        "Apply the boundary condition to the ghost node values in the given Field."
        return "void"

    def setViolationNodes(self,
                          nodeList = "NodeList<DIM>&"):
        "Find any internal nodes that are in violation of this Boundary."
        return "void"

    def updateViolationNodes(self,
                             nodeList = "NodeList<DIM>&"):
        "For the computed set of nodes in violation of the boundary, bring them back into compliance (for the positions and H's."
        return "void"

    @PYB11pycppname("enforceBoundary")
    @PYB11const
    def enforceBoundary1(self,
                         field = "Field<DIM, int>&"):
        "Apply the boundary condition to the violation node values in the given Field."
        return "void"

    @PYB11pycppname("enforceBoundary")
    @PYB11const
    def enforceBoundary2(self,
                         field = "Field<DIM, Scalar>&"):
        "Apply the boundary condition to the violation node values in the given Field."
        return "void"

    @PYB11pycppname("enforceBoundary")
    @PYB11const
    def enforceBoundary3(self,
                         field = "Field<DIM, Vector>&"):
        "Apply the boundary condition to the violation node values in the given Field."
        return "void"

    @PYB11pycppname("enforceBoundary")
    @PYB11const
    def enforceBoundary4(self,
                         field = "Field<DIM, Tensor>&"):
        "Apply the boundary condition to the violation node values in the given Field."
        return "void"

    @PYB11pycppname("enforceBoundary")
    @PYB11const
    def enforceBoundary5(self,
                         field = "Field<DIM, SymTensor>&"):
        "Apply the boundary condition to the violation node values in the given Field."
        return "void"

    @PYB11pycppname("enforceBoundary")
    @PYB11const
    def enforceBoundary6(self,
                         field = "Field<DIM, ThirdRankTensor>&"):
        "Apply the boundary condition to the violation node values in the given Field."
        return "void"
