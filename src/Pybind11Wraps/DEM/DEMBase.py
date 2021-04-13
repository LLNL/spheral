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
  typedef typename %(Dimension)s::Tensor Tensor;
  typedef typename %(Dimension)s::SymTensor SymTensor;
  typedef typename Physics<%(Dimension)s>::TimeStepType TimeStepType;
"""
    
    def pyinit(dataBase = "DataBase<%(Dimension)s>&",
               W = "const TableKernel<%(Dimension)s>&",
               cfl = "const double",
               xmin = "const Vector&",
               xmax = "const Vector&"):
        "DEMBase constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    @PYB11const
    def dt(dataBase = "const DataBase<%(Dimension)s>&", 
           state = "const State<%(Dimension)s>&",
           derivs = "const StateDerivatives<%(Dimension)s>&",
           currentTime = "const Scalar"):
        "Vote on a time step."
        return "TimeStepType"

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
    def evaluateDerivatives(time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<%(Dimension)s>&",
                            state = "const State<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        """Evaluate the derivatives for DEM fields."""
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

    #...........................................................................
    # Properties
    xmin = PYB11property("const Vector&", "xmin", "xmin",
                         returnpolicy="reference_internal",
                         doc="Optional minimum coordinate for bounding box for use generating the mesh for the Voronoi mass density update.")
    xmax = PYB11property("const Vector&", "xmax", "xmax",
                         returnpolicy="reference_internal",
                         doc="Optional maximum coordinate for bounding box for use generating the mesh for the Voronoi mass density update.")
    
    cfl = PYB11property("Scalar", "cfl", "cfl", doc="The Courant-Friedrichs-Lewy timestep limit multiplier")
    
    timeStepMask = PYB11property("const FieldList<%(Dimension)s, int>&",      "timeStepMask", returnpolicy="reference_internal")
    DxDt =         PYB11property("const FieldList<%(Dimension)s, Vector>&",   "DxDt",         returnpolicy="reference_internal")
    DvDt =         PYB11property("const FieldList<%(Dimension)s, Vector>&",   "DvDt",         returnpolicy="reference_internal")
    DomegaDt =     PYB11property("const FieldList<%(Dimension)s, Vector>&",   "DomegaDt",     returnpolicy="reference_internal")
    
#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, DEMBase)
