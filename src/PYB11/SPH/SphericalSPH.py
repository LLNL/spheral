#-------------------------------------------------------------------------------
# SphericalSPH
#-------------------------------------------------------------------------------
from PYB11Generator import *
from SPHBase import *

@PYB11template()            # Override the fact SPHBase is templated
@PYB11template_dict({"Dimension" : "Dim<1>"})
@PYB11module("SpheralSPH")
@PYB11dynamic_attr
class SphericalSPH(SPHBase):

    PYB11typedefs = """
  using Scalar = typename %(Dimension)s::Scalar;
  using Vector = typename %(Dimension)s::Vector;
  using Tensor = typename %(Dimension)s::Tensor;
  using SymTensor = typename %(Dimension)s::SymTensor;
  using TimeStepType = typename Physics<%(Dimension)s>::TimeStepType;
  using PairAccelerationsType = typename SphericalSPH::PairAccelerationsType;
"""
    
    def pyinit(dataBase = "DataBase<%(Dimension)s>&",
               Q = "ArtificialViscosityHandle<%(Dimension)s>&",
               W = "const SphericalKernel&",
               WPi = "const SphericalKernel&",
               cfl = "const double",
               useVelocityMagnitudeForDt = "const bool",
               compatibleEnergyEvolution = "const bool",
               evolveTotalEnergy = "const bool",
               gradhCorrection = "const bool",
               XSPH = "const bool",
               correctVelocityGradient = "const bool",
               sumMassDensityOverAllNodeLists = "const bool",
               densityUpdate = "const MassDensityType",
               epsTensile = "const double",
               nTensile = "const double",
               xmin = "const Vector&",
               xmax = "const Vector&"):
        "Spherical SPH constructor"

    #...........................................................................
    # Virtual methods
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
    def applyGhostBoundaries(state = "State<%(Dimension)s>&",
                             derivs = "StateDerivatives<%(Dimension)s>&"):
        "Apply boundary conditions to the physics specific fields."
        return "void"

    @PYB11virtual
    def enforceBoundaries(state = "State<%(Dimension)s>&",
                          derivs = "StateDerivatives<%(Dimension)s>&"):
        "Enforce boundary conditions for the physics specific fields."
        return "void"

    @PYB11virtual
    @PYB11const
    def label(self):
        return "std::string"

    #...........................................................................
    # Properties
    kernel = PYB11property("const SphericalKernel&", "kernel", doc="The interpolation kernel")
    PiKernel = PYB11property("const SphericalKernel&", "PiKernel", doc="The interpolation kernel for the artificial viscosity")
    Qself = PYB11property("double", "Qself", "Qself", doc="Multiplier for Q self-interaction near the origin")
    pairAccelerations = PYB11property("const PairAccelerationsType&", "pairAccelerations", returnpolicy="reference_internal")
    selfAccelerations = PYB11property("const FieldList<%(Dimension)s, Vector>&", "selfAccelerations", returnpolicy="reference_internal")
