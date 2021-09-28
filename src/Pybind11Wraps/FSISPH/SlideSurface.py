#-------------------------------------------------------------------------------
# FSISPHHydroBase
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Physics import *

@PYB11template("Dimension")
@PYB11module("SpheralFSISPH")
class SlideSurface(Physics):
    "SlideSurface -- basic SPH slide surface feature"

    PYB11typedefs = """
  typedef typename %(Dimension)s::Scalar Scalar;
  typedef typename %(Dimension)s::Vector Vector;
  typedef typename %(Dimension)s::Tensor Tensor;
  typedef typename %(Dimension)s::SymTensor SymTensor;
  typedef typename Physics<%(Dimension)s>::ConstBoundaryIterator ConstBoundaryIterator;
  typedef typename Physics<%(Dimension)s>::TimeStepType TimeStepType;
"""
    
    def pyinit(W = "const TableKernel<%(Dimension)s>&",
               contactTypes = "const vector<int>"):
        "Slide surface constructor"

               
    #...........................................................................
    # Virtual methods
    @PYB11virtual
    @PYB11const
    def evaluateDerivatives(time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<%(Dimension)s>&",
                            state = "const State<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "non-op eval derivs."
        return "void"
    
    @PYB11virtual
    @PYB11const
    def dt(dataBase = "const DataBase<%(Dimension)s>&",
           state = "const State<%(Dimension)s>&",
           derivs = "const StateDerivatives<%(Dimension)s>&",
           currentTime = "const Scalar"):
        "non-op"
        return "TimeStepType"


    @PYB11virtual
    def initializeProblemStartup(dataBase = "DataBase<%(Dimension)s>&"):
        "register the surface normals w/ the database"
        return "void"


    @PYB11virtual
    def initialize(time = "const Scalar",
                   dt = "const Scalar",
                   dataBase = "const DataBase<%(Dimension)s>&",
                   state = "State<%(Dimension)s>&",
                   derivs = "StateDerivatives<%(Dimension)s>&"):
        "calculates surface normals"
        return "void"

    
    @PYB11virtual
    def registerState(dataBase = "DataBase<%(Dimension)s>&",
                      state = "State<%(Dimension)s>&"):
        "register the surface normals"
        return "void"


    @PYB11virtual
    def registerDerivatives(dataBase = "DataBase<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "non-op place filler"
        return "void"

    @PYB11const
    def isSlideSurface(nodeListi = "const int",
                       nodeListj = "const int"):
        "returns true if slide surface between nodelisti and j"
        return "bool"

    @PYB11const
    def slideCorrection(nodeListi = "const int",
                        i = "const int",
                        nodeListj = "const int",
                        j = "const int",
                        vi = "const Vector",
                        vj = "const Vector"):
        "returns true if slide surface between nodelisti and j"
        return "Scalar"

    @PYB11virtual
    @PYB11const
    def label(self):
        return "std::string"
    #...........................................................................
    # Properties
    surfaceNormals = PYB11property("const FieldList<%(Dimension)s, Vector>&", "surfaceNormals", returnpolicy="reference_internal")
    surfaceFraction = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "surfaceFraction", returnpolicy="reference_internal")
    surfaceSmoothness = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "surfaceSmoothness", returnpolicy="reference_internal")
    numNodeLists = PYB11property("int", "numNodeLists", "numNodeLists", doc="number of nodelists.")
    isActive = PYB11property("bool", "isActive", "isActive", doc="switch if we have a slide.")
    

