#-------------------------------------------------------------------------------
# Provide the abstact Boundary interface, suitable for injection to our actual
# classes.
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11ignore
class BoundaryAbstractMethods:

    def setGhostNodes(self,
                      nodeList = "NodeList<%(Dimension)s>&"):
        "Use the given NodeList's neighbor object to select the ghost nodes."
        return "void"

    def updateGhostNodes(self,
                         nodeList = "NodeList<%(Dimension)s>&"):
        "For the computed set of ghost nodes, set the positions and H's."
        return "void"

    def setViolationNodes(self,
                          nodeList = "NodeList<%(Dimension)s>&"):
        "Find any internal nodes that are in violation of this Boundary."
        return "void"

    def updateViolationNodes(self,
                             nodeList = "NodeList<%(Dimension)s>&"):
        "For the computed set of nodes in violation of the boundary, bring them back into compliance (for the positions and H's."
        return "void"

