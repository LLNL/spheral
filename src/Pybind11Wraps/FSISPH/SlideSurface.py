#-------------------------------------------------------------------------------
# FSISPHHydroBase
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11template("Dimension")
@PYB11module("SpheralFSISPH")
class SlideSurface:
    "SlideSurface -- basic SPH slide surface feature"

    PYB11typedefs = """
  typedef typename %(Dimension)s::Scalar Scalar;
  typedef typename %(Dimension)s::Vector Vector;
  typedef typename %(Dimension)s::Tensor Tensor;
  typedef typename %(Dimension)s::SymTensor SymTensor;
  typedef typename Physics<%(Dimension)s>::ConstBoundaryIterator ConstBoundaryIterator;
  typedef typename Physics<%(Dimension)s>::TimeStepType TimeStepType;
"""
    
    def pyinit(dataBase = "DataBase<%(Dimension)s>&",
               contactTypes = "const vector<int>",
               surfaceNormalMethod = "const SurfaceNormalMethod",
               normalsAreSmoothed = "const bool",
               gradientsAreCorrected = "const bool"):
        "Slide surface constructor"

               
    #...........................................................................
    # Virtual methods
    # @PYB11virtual
    # def initializeProblemStartup(dataBase = "DataBase<%(Dimension)s>&"):
    #     "register the surface normals w/ the database"
    #     return "void"
    
    # @PYB11virtual
    # def registerState(dataBase = "DataBase<%(Dimension)s>&",
    #                   state = "State<%(Dimension)s>&"):
    #     "register the surface normals, frac, and smoothness"
    #     return "void"


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

    #...........................................................................
    # Properties
    #surfaceNormals = PYB11property("const FieldList<%(Dimension)s, Vector>&", "surfaceNormals", returnpolicy="reference_internal")
    #surfaceFraction = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "surfaceFraction", returnpolicy="reference_internal")
    #surfaceSmoothness = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "surfaceSmoothness", returnpolicy="reference_internal")
    numNodeLists = PYB11property("int", "numNodeLists", "numNodeLists", doc="number of nodelists.")
    isActive = PYB11property("bool", "isActive", "isActive", doc="switch if we have a slide.")
    normalsAreSmoothed = PYB11property("bool", "normalsAreSmoothed", "normalsAreSmoothed", doc="do we want to smooth our normals.")
    gradientsAreCorrected = PYB11property("bool", "gradientsAreCorrected", "gradientsAreCorrected", doc="switch to correct kernel gradients in interface normal calculation.")
    surfaceNormalMethod = PYB11property("SurfaceNormalMethod", "surfaceNormalMethod", "surfaceNormalMethod",
                                    doc="Flag to select what material we use for normal calculation, same, diff, all")
