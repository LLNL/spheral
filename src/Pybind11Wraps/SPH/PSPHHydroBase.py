#-------------------------------------------------------------------------------
# PSPHHydroBase
#-------------------------------------------------------------------------------
from PYB11Generator import *
from SPHHydroBase import *
from RestartMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralSPH")
class PSPHHydroBase(SPHHydroBase):

    typedefs = """
  typedef %(Dimension)s DIM;
  typedef typename %(Dimension)s::Scalar Scalar;
  typedef typename %(Dimension)s::Vector Vector;
  typedef typename %(Dimension)s::Tensor Tensor;
  typedef typename %(Dimension)s::SymTensor SymTensor;
  typedef typename Physics<%(Dimension)s>::TimeStepType TimeStepType;
"""
    
    def pyinit(smoothingScaleMethod = "const SmoothingScaleBase<DIM>&",
               Q = "ArtificialViscosity<DIM>&",
               W = "const TableKernel<DIM>&",
               WPi = "const TableKernel<DIM>&",
               filter = "const double",
               cfl = "const double",
               useVelocityMagnitudeForDt = "const bool",
               compatibleEnergyEvolution = "const bool",
               evolveTotalEnergy = "const bool",
               XSPH = "const bool",
               correctVelocityGradient = "const bool",
               HopkinsConductivity = "const bool",
               sumMassDensityOverAllNodeLists = "const bool",
               densityUpdate = "const MassDensityType",
               HUpdate = "const HEvolutionType",
               xmin = "const Vector&",
               xmax = "const Vector&"):
        "PSPHHydroBase constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def initializeProblemStartup(dataBase = "DataBase<DIM>&"):
        "Tasks we do once on problem startup."
        return "void"

    @PYB11virtual 
    def registerState(dataBase = "DataBase<DIM>&",
                      state = "State<DIM>&"):
        "Register the state Hydro expects to use and evolve."
        return "void"

    @PYB11virtual 
    def preStepInitialize(dataBase = "const DataBase<DIM>&",
                          state = "State<DIM>&",
                          derivs = "StateDerivatives<DIM>&"):
        "Pre-step initializations."
        return "void"

    @PYB11virtual
    @PYB11const
    def evaluateDerivatives(time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<DIM>&",
                            state = "const State<DIM>&",
                            derivs = "StateDerivatives<DIM>&"):
        """Evaluate the derivatives for the principle hydro 
mass density, velocity, and specific thermal energy."""
        return "void"

    @PYB11virtual
    @PYB11const
    def finalizeDerivatives(time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<DIM>&",
                            state = "const State<DIM>&",
                            derivs = "StateDerivatives<DIM>&"):
        "Finalize the derivatives."
        return "void"

    @PYB11virtual
    def postStateUpdate(time = "const Scalar",
                        dt = "const Scalar",
                        dataBase = "const DataBase<DIM>&",
                        state = "State<DIM>&",
                        derivs = "StateDerivatives<DIM>&"):
        "Post-state update. For PSPH this is where we recompute the PSPH pressure and corrections."
        return "void"
               
    @PYB11virtual
    def applyGhostBoundaries(state = "State<DIM>&",
                             derivs = "StateDerivatives<DIM>&"):
        "Apply boundary conditions to the physics specific fields."
        return "void"

    @PYB11virtual
    def enforceBoundaries(state = "State<DIM>&",
                          derivs = "StateDerivatives<DIM>&"):
        "Enforce boundary conditions for the physics specific fields."
        return "void"

    #...........................................................................
    # Properties
    HopkinsConductivity = PYB11property("bool", "HopkinsConductivity", "HopkinsConductivity", 
                                        doc="Flag determining if we're applying Hopkins 2014 conductivity.")

    gamma =          PYB11property("const FieldList<DIM, Scalar>&", "gamma",          returnpolicy="reference_internal")
    PSPHcorrection = PYB11property("const FieldList<DIM, Scalar>&", "PSPHcorrection", returnpolicy="reference_internal")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, PSPHHydroBase)
