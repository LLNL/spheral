#-------------------------------------------------------------------------------
# GSPHHydroBase
#-------------------------------------------------------------------------------
from PYB11Generator import *
from GenericRiemannHydro import *
from RestartMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralGSPH")
@PYB11dynamic_attr
class MFVHydroBase(GenericRiemannHydro):

    PYB11typedefs = """
  typedef typename %(Dimension)s::Scalar Scalar;
  typedef typename %(Dimension)s::Vector Vector;
  typedef typename %(Dimension)s::Tensor Tensor;
  typedef typename %(Dimension)s::SymTensor SymTensor;
  typedef typename Physics<%(Dimension)s>::TimeStepType TimeStepType;
"""
    
    def pyinit(smoothingScaleMethod = "const SmoothingScaleBase<%(Dimension)s>&",
               dataBase = "DataBase<%(Dimension)s>&",
               riemannSolver = "RiemannSolverBase<%(Dimension)s>&",
               W = "const TableKernel<%(Dimension)s>&",
               epsDiffusionCoeff = "const Scalar",
               cfl = "const double",
               useVelocityMagnitudeForDt = "const bool",
               compatibleEnergyEvolution = "const bool",
               evolveTotalEnergy = "const bool",
               XSPH = "const bool",
               correctVelocityGradient = "const bool",
               nodeMotionCoefficient = "const double",
               nodeMotionType = "const NodeMotionType",
               gradType = "const GradientType",
               densityUpdate = "const MassDensityType",
               HUpdate = "const HEvolutionType",
               epsTensile = "const double",
               nTensile = "const double",
               xmin = "const Vector&",
               xmax = "const Vector&"):
        "GSPHHydroBase constructor"

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
        "Initialize the Hydro before we start a derivative evaluation."
        return "void"
                       
    @PYB11virtual
    @PYB11const
    def evaluateDerivatives(time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<%(Dimension)s>&",
                            state = "const State<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        """Evaluate the derivatives for the principle hydro 
mass density, velocity, and specific thermal energy."""
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

    DvolumeDt = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "DvolumeDt", returnpolicy="reference_internal")
    nodeMotionCoefficient = PYB11property("double", "nodeMotionCoefficient", "nodeMotionCoefficient",doc="multiplier for XSPH and Fician node motion schemes.")
    nodeMotionType = PYB11property("NodeMotionType","nodeMotionType","nodeMotionType")
    
#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, MFVHydroBase)
