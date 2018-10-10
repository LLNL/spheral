#-------------------------------------------------------------------------------
# SVPHFactedHydroBase
#-------------------------------------------------------------------------------
from PYB11Generator import *

from GenericHydro import *

@PYB11template("Dimension")
class SVPHFacetedHydroBase(GenericHydro):
    "SVPHFacetedHydroBase -- The fluid SVPH faceted hydro algorithm"

    typedefs = """
    typedef %(Dimension)s DIM;
    typedef typename DIM::Scalar Scalar;
    typedef typename DIM::Vector Vector;
    typedef typename DIM::Tensor Tensor;
    typedef typename DIM::SymTensor SymTensor;
    typedef typename DIM::ThirdRankTensor ThirdRankTensor;
    typedef typename DIM::FourthRankTensor FourthRankTensor;
    typedef typename DIM::FifthRankTensor FifthRankTensor;
    typedef typename Physics<DIM>::TimeStepType TimeStepType;
"""

    def pyinit(self,
               smoothingScaleMethod = "const SmoothingScaleBase<Dimension>&",
               W = "const TableKernel<Dimension>&",
               Q = "ArtificialViscosity<Dimension>&",
               cfl = "const double",
               useVelocityMagnitudeForDt = "const bool",
               compatibleEnergyEvolution = "const bool",
               XSVPH = "const bool",
               linearConsistent = "const bool",
               generateVoid = "const bool",
               densityUpdate = "const MassDensityType",
               HUpdate = "const HEvolutionType",
               fcentroidal = "const Scalar",
               fcellPressure = "const Scalar",
               xmin = "const Vector&",
               xmax = "const Vector&"):
        "SVPH Constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def initializeProblemStartup(self, dataBase = "DataBase<DIM>&"):
        "Tasks we do once on problem startup."
        return "void"

    @PYB11virtual
    @PYB11const
    def dt(self,
           dataBase = "const DataBase<DIM>&",
           state = "const State<DIM>&",
           derivs = "const StateDerivatives<DIM>&",
           currentTime = "const Scalar"):
        "Vote on a time step."
        return "TimeStepType"

    @PYB11virtual 
    def registerState(self,
                      dataBase = "DataBase<DIM>&",
                      state = "State<DIM>&"):
        "Register the state Hydro expects to use and evolve."
        return "void"

    @PYB11virtual
    def registerDerivatives(self,
                            dataBase = "DataBase<DIM>&",
                            derivs = "StateDerivatives<DIM>&"):
        "Register the derivatives/change fields for updating state."
        return "void"

    @PYB11virtual
    def initialize(self,
                   time = "const Scalar",
                   dt = "const Scalar",
                   dataBase = "const DataBase<DIM>&",
                   state = "State<DIM>&",
                   derivs = "StateDerivatives<DIM>&"):
        "Initialize the Hydro before we start a derivative evaluation."
        return "void"
                          
    @PYB11virtual
    @PYB11const
    def evaluateDerivatives(self,
                            time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<DIM>&",
                            state = "const State<DIM>&",
                            derivs = "StateDerivatives<DIM>&"):
        """Evaluate the derivatives for the principle hydro variables:
mass density, velocity, and specific thermal energy."""
        return "void"

    @PYB11virtual
    @PYB11const
    def finalizeDerivatives(self,
                            time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<DIM>&",
                            state = "const State<DIM>&",
                            derivs = "StateDerivatives<DIM>&"):
        "Finalize the derivatives."
        return "void"

    @PYB11virtual
    def finalize(self,
                 time = "const Scalar",
                 dt = "const Scalar",
                 dataBase = "DataBase<DIM>&",
                 state = "State<DIM>&",
                 derivs = "StateDerivatives<DIM>&"):
        "Finalize the hydro at the completion of an integration step."
        return "void"
                  
    @PYB11virtual
    @PYB11const
    def requireConnectivity(self):
        "This algorithm does not use node->node connectivity."
        return "bool"

    @PYB11virtual
    def applyGhostBoundaries(self,
                             state = "State<DIM>&",
                             derivs = "StateDerivatives<DIM>&"):
        "Apply boundary conditions to the physics specific fields."
        return "void"

    @PYB11virtual
    def enforceBoundaries(self,
                          state = "State<DIM>&",
                          derivs = "StateDerivatives<DIM>&"):
        "Enforce boundary conditions for the physics specific fields."
        return "void"

    #...........................................................................
    # Properties
    
