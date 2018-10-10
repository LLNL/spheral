#-------------------------------------------------------------------------------
# CRKSPHVariant
#-------------------------------------------------------------------------------
from PYB11Generator import *
from CRKSPHHydroBase import *

@PYB11template("Dimension")
@PYB11module("SpheralCRKSPH")
class CRKSPHVariant(CRKSPHHydroBase):
    "CRKSPHVariant -- A development variant of CRKSPH for experimentation."

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
               nTensile = "const double"):
        "Constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def initializeProblemStartup(self, dataBase = "DataBase<DIM>&"):
        "Tasks we do once on problem startup."
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
