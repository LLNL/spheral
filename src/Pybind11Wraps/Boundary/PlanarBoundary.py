#-------------------------------------------------------------------------------
# PlanarBoundary base class
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Boundary import *
from RestartMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralBoundary")
class PlanarBoundary(Boundary):

    typedefs = """
    typedef %(Dimension)s DIM;
    typedef typename DIM::Scalar Scalar;
    typedef typename DIM::Vector Vector;
    typedef typename DIM::Tensor Tensor;
    typedef typename DIM::SymTensor SymTensor;
    typedef typename DIM::ThirdRankTensor ThirdRankTensor;
    typedef typename Boundary<DIM>::BoundaryNodes BoundaryNodes;
    typedef GeomPlane<DIM> Plane;
"""

    #...........................................................................
    # Constructors
    def pyinit(self):
        "Default constructor"

    def pyinit1(self,
                enterPlane = "const Plane&",
                exitPlane = "const Plane&"):
        "Construct with enter/exit mapping planes"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def setGhostNodes(self,
                      nodeList = "NodeList<DIM>&"):
        "Set ghost nodes for the NodeList"
        return "void"

    @PYB11virtual
    def updateGhostNodes(self,
                         nodeList = "NodeList<DIM>&"):
        "Update position and H for ghost nodes for the NodeList"
        return "void"

    @PYB11virtual
    def setViolationNodes(self,
                      nodeList = "NodeList<DIM>&"):
        "Set violation nodes for the NodeList"
        return "void"

    @PYB11virtual
    def updateViolationNodes(self,
                         nodeList = "NodeList<DIM>&"):
        "Update nodes in violation of this Boundary for the NodeList"
        return "void"
                    
    @PYB11virtual
    @PYB11const
    def enterPlane(self):
        "Get the enter plane"
        return "const Plane&"
        
    @PYB11virtual
    def setEnterPlane(self, plane="const Plane&"):
        "Set the enter plane"
        return "void"

    @PYB11virtual
    @PYB11const
    def exitPlane(self):
        "Get the exit plane"
        return "const Plane&"
        
    @PYB11virtual
    def setExitPlane(self, plane="const Plane&"):
        "Set the exit plane"
        return "void"

    @PYB11virtual
    @PYB11const
    def valid(self):
        return "bool"

    @PYB11virtual
    @PYB11const
    def clip(self,
             xmin = "Vector&",
             xmax = "Vector&"):
        "Override the boundary clip"
        return "void"

    #...........................................................................
    # Methods
    def setGhostNodes(self,
                      nodeList = "NodeList<DIM>&",
                      presetControlNodes = "const std::vector<int>&"):
        "Set the ghost nodes for a predefined set of control nodes"
        return "void"

    @PYB11const
    def mapPosition(self,
                    position = "const Vector&",
                    enterPlane = "const Plane&",
                    exitPlane = "const Plane&"):
        "Map a position through the enter/exit plane transformation"
        return "Vector"

    @PYB11const
    def facesOnPlane(self,
                     mesh = "const Mesh<DIM>&",
                     plane = "const Plane&",
                     tol = "const Scalar"):
        "Provide a method to identify tessellation faces on a plane."
        return "std::vector<unsigned>"

#-------------------------------------------------------------------------------
# Inject restart methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, PlanarBoundary)
