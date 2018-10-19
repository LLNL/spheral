#-------------------------------------------------------------------------------
# Boundary base class
#-------------------------------------------------------------------------------
from PYB11Generator import *
from BoundaryAbstractMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralBoundary")
class Boundary:

    typedefs = """
    typedef %(Dimension)s DIM;
    typedef typename DIM::Scalar Scalar;
    typedef typename DIM::Vector Vector;
    typedef typename DIM::Tensor Tensor;
    typedef typename DIM::SymTensor SymTensor;
    typedef typename DIM::ThirdRankTensor ThirdRankTensor;
    typedef typename Boundary<DIM>::BoundaryNodes BoundaryNodes;
"""

    #...........................................................................
    # Constructors
    def pyinit(self):
        "Default constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def setAllGhostNodes(self,
                         dataBase = "DataBase<DIM>&"):
        "Recreate ghost nodes for this Boundary for all NodeLists in the DataBase."
        return "void"

    @PYB11virtual
    def setAllViolationNodes(self,
                             dataBase = "DataBase<DIM>&"):
        "Select any nodes based in the FluidNodeLists in the given DataBase that are in violation of boundary condition."
        return "void"

    @PYB11virtual
    def cullGhostNodes(self,
                       flagSet = "const FieldList<DIM, int>&",
                       old2newIndexMap = "FieldList<DIM, int>&",
                       numNodesRemoved = "std::vector<int>&"):
        "Use a set of flags to cull out inactive ghost nodes."
        return "void"

    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary20(self,
                             field = "Field<DIM, std::vector<Scalar>>&"):
        "We provide default copies for arrays of values, but descendants can override these."
        return "void"

    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary21(self,
                             field = "Field<DIM, std::vector<Vector>>&"):
        "We provide default copies for arrays of values, but descendants can override these."
        return "void"

    @PYB11virtual
    def initializeProblemStartup(self):
        "Some boundaries need to know when a problem is starting up and all the physics packages have been initialized."
        return "void"

    @PYB11virtual
    @PYB11const
    def finalizeGhostBoundary(self):
        "Provide an optional hook that is to be called when all ghost boundaries are to have been set."
        return "void"

    @PYB11virtual
    def reset(self,
              dataBase = "const DataBase<DIM>&"):
        "Overridable hook for clearing out the boundary condition."
        return "void"

    @PYB11virtual
    @PYB11const
    def numGhostNodes(self,):
        "Report the number of ghost nodes in this boundary."
        return "int"

    @PYB11virtual
    @PYB11const
    def clip(self,
             xmin = "Vector&",
             xmax = "Vector&"):
        "Optionally the boundary can clip an input box range. Defaults to no-op."
        return "void"

    @PYB11virtual
    @PYB11const
    def meshGhostNodes(self):
        "Some boundaries have ghosts we should exclude from tessellations. Provide a hook to note such cases."
        return "bool"

    @PYB11virtual
    def addNodeList(self, nodeList="NodeList<DIM>&"):
        return "void"
                    
    #...........................................................................
    # Methods
    @PYB11const
    def haveNodeList(self, nodeList="const NodeList<DIM>&"):
        "Check if we have entries for the given NodeList."
        return "bool"

    @PYB11returnpolicy("reference_internal")
    def accessBoundaryNodes(self, nodeList="NodeList<DIM>&"):
        return "BoundaryNodes&"

    @PYB11returnpolicy("reference_internal")
    @PYB11const
    def controlNodes(self, nodeList="NodeList<DIM>&"):
        "Control nodes for a given NodeList"
        return "const std::vector<int>&"

    @PYB11returnpolicy("reference_internal")
    @PYB11const
    def ghostNodes(self, nodeList="NodeList<DIM>&"):
        "Ghost nodes for a given NodeList"
        return "const std::vector<int>&"

    @PYB11returnpolicy("reference_internal")
    @PYB11const
    def violationNodes(self, nodeList="NodeList<DIM>&"):
        "Violation nodes for a given NodeList"
        return "const std::vector<int>&"

    #...........................................................................
    # applyFieldListGhostBoundary/enforceFieldListBoundary
    @PYB11template("Value")
    @PYB11const
    def applyFieldListGhostBoundary(self, fieldList="FieldList<DIM, %(Value)s>&"):
        "Apply ghost boundary to all Fields in FieldList"
        return "void"

    @PYB11template("Value")
    @PYB11const
    def enforceFieldListBoundary(self, fieldList="FieldList<DIM, %(Value)s>&"):
        "Enforce boundary on all Fields in FieldList"
        return "void"

    for T in ("int", "Scalar", "Vector", "Tensor", "SymTensor", "ThirdRankTensor"):
        exec('''
aflgb%(T)s = PYB11TemplateMethod(applyFieldListGhostBoundary, template_parameters="%(T)s", pyname="applyFieldListGhostBoundary")
eflgb%(T)s = PYB11TemplateMethod(enforceFieldListBoundary, template_parameters="%(T)s", pyname="enforceFieldListBoundary")
''' % {"T" : T})

    #---------------------------------------------------------------------------
    # BoundaryNodes
    #---------------------------------------------------------------------------
    class BoundaryNodes:

        controlNodes = PYB11readwrite()
        ghostNodes = PYB11readwrite()
        violationNodes = PYB11readwrite()

#-------------------------------------------------------------------------------
# Inject abstract interface
#-------------------------------------------------------------------------------
PYB11inject(BoundaryAbstractMethods, Boundary, pure_virtual=True)
