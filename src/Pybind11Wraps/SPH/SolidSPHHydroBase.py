#-------------------------------------------------------------------------------
# SolidSPHHydroBase
#-------------------------------------------------------------------------------
from PYB11Generator import *
from SPHHydroBase import *
from RestartMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralSPH")
class SolidSPHHydroBase(SPHHydroBase):
    "SolidSPHHydroBase -- The SPH/ASPH solid material hydrodynamic package for Spheral++."

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
               WGrad = "const TableKernel<DIM>&",
               filter = "const double",
               cfl = "const double",
               useVelocityMagnitudeForDt = "const bool",
               compatibleEnergyEvolution = "const bool",
               evolveTotalEnergy = "const bool",
               gradhCorrection = "const bool",
               XSPH = "const bool",
               correctVelocityGradient = "const bool",
               sumMassDensityOverAllNodeLists = "const bool",
               densityUpdate = "const MassDensityType",
               HUpdate = "const HEvolutionType",
               epsTensile = "const double",
               nTensile = "const double",
               damageRelieveRubble = "const bool",
               xmin = "const Vector&",
               xmax = "const Vector&"):
        "SolidSPHHydroBase constructor"

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
    def registerDerivatives(dataBase = "DataBase<DIM>&",
                            derivs = "StateDerivatives<DIM>&"):
        "Register the derivatives/change fields for updating state."
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
    GradKernel = PYB11property("const TableKernel<DIM>&", "GradKernel",
                               doc="Kernel for estimating velocity gradient")
    damageRelieveRubble = PYB11property("bool", "damageRelieveRubble", "damageRelieveRubble", 
                                        doc="Control whether allow damaged material to have stress relieved.")

    DdeviatoricStressDt =  PYB11property("const FieldList<DIM, SymTensor>&", "DdeviatoricStressDt", returnpolicy="reference_internal")
    bulkModulus =          PYB11property("const FieldList<DIM, Scalar>&",    "bulkModulus",         returnpolicy="reference_internal")
    shearModulus =         PYB11property("const FieldList<DIM, Scalar>&",    "shearModulus",        returnpolicy="reference_internal")
    yieldStrength =        PYB11property("const FieldList<DIM, Scalar>&",    "yieldStrength",       returnpolicy="reference_internal")
    plasticStrain0 =       PYB11property("const FieldList<DIM, Scalar>&",    "plasticStrain0",      returnpolicy="reference_internal")
    Hfield0 =              PYB11property("const FieldList<DIM, SymTensor>&", "Hfield0",             returnpolicy="reference_internal")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, SPHHydroBase)
