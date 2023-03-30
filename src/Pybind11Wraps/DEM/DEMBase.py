#-------------------------------------------------------------------------------
# SPHHydroBase
#-------------------------------------------------------------------------------
from PYB11Generator import *
from RestartMethods import *
from Physics import *

@PYB11template("Dimension")
@PYB11module("SpheralDEM")
class DEMBase(Physics):

    PYB11typedefs = """
  typedef typename %(Dimension)s::Scalar Scalar;
  typedef typename %(Dimension)s::Vector Vector;
  typedef typename DEMDimension<%(Dimension)s>::AngularVector RotationType;
  typedef typename %(Dimension)s::Tensor Tensor;
  typedef typename %(Dimension)s::SymTensor SymTensor;
  typedef typename Physics<%(Dimension)s>::TimeStepType TimeStepType;
"""
    
    def pyinit(dataBase = "const DataBase<%(Dimension)s>&",
               stepsPerCollision = "const double",
               xmin = "const Vector&",
               xmax = "const Vector&"):
        "DEMBase constructor"

    #...........................................................................
    # Virtual methods

    @PYB11virtual
    def initializeProblemStartup(dataBase = "DataBase<%(Dimension)s>&"):
        "Tasks we do once on problem startup."
        return "void"

    @PYB11virtual 
    def registerState(dataBase = "DataBase<%(Dimension)s>&",
                      state = "State<%(Dimension)s>&"):
        "Register the state Hydro expects to use and evolve."
        return "void"

    @PYB11virtual
    def registerDerivatives(dataBase = "DataBase<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Register the derivatives/change fields for updating state."
        return "void"

    @PYB11virtual
    def preStepInitialize(self,
                          dataBase = "const DataBase<%(Dimension)s>&", 
                          state = "State<%(Dimension)s>&",
                          derivs = "StateDerivatives<%(Dimension)s>&"):
        "Optional hook to be called at the beginning of a time step."
        return "void"

    @PYB11virtual
    def initialize(time = "const Scalar",
                   dt = "const Scalar",
                   dataBase = "const DataBase<%(Dimension)s>&",
                   state = "State<%(Dimension)s>&",
                   derivs = "StateDerivatives<%(Dimension)s>&"):
        "Initialize the DEM before we start a derivative evaluation."
        return "void"
                       

    @PYB11virtual
    @PYB11const
    def finalizeDerivatives(time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<%(Dimension)s>&",
                            state = "const State<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Finalize the derivatives."
        return "void"

    @PYB11virtual
    def applyGhostBoundaries(state = "State<%(Dimension)s>&",
                             derivs = "StateDerivatives<%(Dimension)s>&"):
        "Apply boundary conditions to the physics specific fields."
        return "void"

    @PYB11virtual
    def enforceBoundaries(state = "State<%(Dimension)s>&",
                          derivs = "StateDerivatives<%(Dimension)s>&"):
        "Enforce boundary conditions for the physics specific fields."
        return "void"

    @PYB11const
    def momentOfInertia(massi = "const Scalar",
                        partialRadiusi = "const Scalar"):
        return "Scalar"
    #...........................................................................
    # Properties
    xmin = PYB11property("const Vector&", "xmin", "xmin",
                         returnpolicy="reference_internal",
                         doc="Optional minimum coordinate for bounding box for use generating the mesh for the Voronoi mass density update.")
    xmax = PYB11property("const Vector&", "xmax", "xmax",
                         returnpolicy="reference_internal",
                         doc="Optional maximum coordinate for bounding box for use generating the mesh for the Voronoi mass density update.")
    
    stepsPerCollision = PYB11property("Scalar", "stepsPerCollision", "stepsPerCollision", doc="number of steps resolving the collision time scale")
    
    cycle = PYB11property("int", "cycle", "cycle", doc="tracks what cycle we are on.")
    contactRemovalFrequency = PYB11property("int", "contactRemovalFrequency", "contactRemovalFrequency", doc="every n-cycles we prune our contact list to only the active contacts. This is the frequency.")
    firstCycle = PYB11property("bool", "firstCycle", "firstCycle", doc="boolean flag thats true for the first cycle of a run.")
    
    timeStepMask =  PYB11property("const FieldList<%(Dimension)s, int>&", "timeStepMask", returnpolicy="reference_internal")
    DxDt =          PYB11property("const FieldList<%(Dimension)s, Vector>&","DxDt", returnpolicy="reference_internal")
    DvDt =          PYB11property("const FieldList<%(Dimension)s, Vector>&", "DvDt", returnpolicy="reference_internal")
    DomegaDt =      PYB11property("const FieldList<%(Dimension)s, RotationType>&","DomegaDt", returnpolicy="reference_internal")
    omega =         PYB11property("const FieldList<%(Dimension)s, RotationType>&","omega", returnpolicy="reference_internal")
    uniqueIndices = PYB11property("const FieldList<%(Dimension)s, int>&","uniqueIndices", returnpolicy="reference_internal")

    equilibriumOverlap = PYB11property("const FieldList<%(Dimension)s, vector<Scalar>>&","equilibriumOverlap", returnpolicy="reference_internal")
    neighborIndices = PYB11property("const FieldList<%(Dimension)s, vector<int>>&","neighborIndices", returnpolicy="reference_internal")
    shearDisplacement = PYB11property("const FieldList<%(Dimension)s, vector<Vector>>&","shearDisplacement", returnpolicy="reference_internal")
    DDtShearDisplacement = PYB11property("const FieldList<%(Dimension)s, vector<Vector>>&","DDtShearDisplacement", returnpolicy="reference_internal")
    isActiveContact = PYB11property("const FieldList<%(Dimension)s, vector<int>>&","isActiveContact", returnpolicy="reference_internal")
    
#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, DEMBase)

