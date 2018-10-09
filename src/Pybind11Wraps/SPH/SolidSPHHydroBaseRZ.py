#-------------------------------------------------------------------------------
# SolidSPHHydroBaseRZ
#-------------------------------------------------------------------------------
from PYB11Generator import *
from SolidSPHHydroBase import *
from RestartMethods import *

@PYB11template()            # Override the fact SolidSPHHydroBase is templated
@PYB11template_dict({"Dimension" : "Dim<2>"})
@PYB11module("SpheralSPH")
class SolidSPHHydroBaseRZ(SolidSPHHydroBase):

    typedefs = """
  typedef Dim<2> DIM;
  typedef typename DIM::Scalar Scalar;
  typedef typename DIM::Vector Vector;
  typedef typename DIM::Tensor Tensor;
  typedef typename DIM::SymTensor SymTensor;
  typedef typename Physics<DIM>::TimeStepType TimeStepType;
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
        "SolidSPHHydroBaseRZ constructor"

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
    def finalize(time = "const Scalar",
                 dt = "const Scalar",
                 dataBase = "DataBase<DIM>&",
                 state = "State<DIM>&",
                 derivs = "StateDerivatives<DIM>&"):
        "Finalize the hydro at the completion of an integration step."
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
    deviatoricStressTT =  PYB11property("const FieldList<DIM, Scalar>&", "deviatoricStressTT", returnpolicy="reference_internal")
    DdeviatoricStressTTDt =  PYB11property("const FieldList<DIM, Scalar>&", "DdeviatoricStressTTDt", returnpolicy="reference_internal")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, SolidSPHHydroBaseRZ)
