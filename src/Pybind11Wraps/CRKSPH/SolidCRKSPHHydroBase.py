#-------------------------------------------------------------------------------
# SolidCRKSPHHydroBase
#-------------------------------------------------------------------------------
from PYB11Generator import *
from CRKSPHHydroBase import *
from RestartMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralCRKSPH")
class SolidCRKSPHHydroBase(CRKSPHHydroBase):
    "SolidCRKSPHHydroBase -- The CRKSPH/ACRKSPH solid material hydrodynamic package for Spheral++."

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
               smoothingScaleMethod = "const SmoothingScaleBase<DIM>&",
               Q = "ArtificialViscosity<DIM>&",
               W = "const TableKernel<DIM>&",
               WPi = "const TableKernel<DIM>&",
               filter = "const double",
               cfl = "const double",
               useVelocityMagnitudeForDt = "const bool",
               compatibleEnergyEvolution = "const bool",
               evolveTotalEnergy = "const bool",
               XSPH = "const bool",
               densityUpdate = "const MassDensityType",
               HUpdate = "const HEvolutionType",
               correctionOrder = "const CRKOrder",
               volumeType = "const CRKVolumeType",
               epsTensile = "const double",
               nTensile = "const double",
               damageRelieveRubble = "const bool"):
        "Constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def initializeProblemStartup(self, dataBase = "DataBase<DIM>&"):
        "Tasks we do once on problem startup."
        return "void"

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
    damageRelieveRubble = PYB11property("bool", "damageRelieveRubble", "damageRelieveRubble",
                                        doc="Control whether allow damaged material to have stress relieved.")

    DdeviatoricStressDt = PYB11property("const FieldList<DIM, SymTensor>&", "DdeviatoricStressDt", returnpolicy="reference_internal")
    bulkModulus = PYB11property("const FieldList<DIM, Scalar>&", "bulkModulus", returnpolicy="reference_internal")
    shearModulus = PYB11property("const FieldList<DIM, Scalar>&", "shearModulus", returnpolicy="reference_internal")
    yieldStrength = PYB11property("const FieldList<DIM, Scalar>&", "yieldStrength", returnpolicy="reference_internal")
    plasticStrain0 = PYB11property("const FieldList<DIM, Scalar>&", "plasticStrain0", returnpolicy="reference_internal")
    Hfield0 = PYB11property("const FieldList<DIM, SymTensor>&", "Hfield0", returnpolicy="reference_internal")
    fragIDs = PYB11property("const FieldList<DIM, int>&", "fragIDs", returnpolicy="reference_internal")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, CRKSPHHydroBase)

