#-------------------------------------------------------------------------------
# TensorDamageModel
#-------------------------------------------------------------------------------
from PYB11Generator import *
from DamageModel import *
from RestartMethods import *

@PYB11template("Dimension")
class TensorDamageModel(DamageModel):
    """TensorDamageModel -- Base class for the tensor damage physics models.
This class does not know how to seed the flaw distribution -- that is 
required of descendant classes."""

    typedefs = """
    typedef %(Dimension)s DIM;
    typedef typename DIM::Scalar Scalar;
    typedef typename DIM::Vector Vector;
    typedef typename DIM::Tensor Tensor;
    typedef typename DIM::SymTensor SymTensor;
    typedef typename Physics<DIM>::TimeStepType TimeStepType;

    typedef Field<DIM, std::vector<double> > FlawStorageType;
"""

    def pyinit(self,
               nodeList = "SolidNodeList<DIM>&",
               strainAlgorithm = "const TensorStrainAlgorithm",
               effDamageAlgorithm = "const EffectiveDamageAlgorithm",
               useDamageGradient = "const bool",
               W = "const TableKernel<DIM>&",
               crackGrowthMultiplier = "const double",
               flawAlgorithm = "const EffectiveFlawAlgorithm",
               criticalDamageThreshold = "const double",
               damageInCompression = "const bool",
               flaws = "const FlawStorageType&"):
        "Constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual 
    @PYB11const
    def evaluateDerivatives(self,
                            time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<DIM>&",
                            state = "const State<DIM>&",
                            derivs = "StateDerivatives<DIM>&"):
        "Compute the derivatives."
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
        "Register our state."
        return "void"

    @PYB11virtual
    def registerDerivatives(self,
                            dataBase = "DataBase<DIM>&",
                            derivs = "StateDerivatives<DIM>&"):
        "Register the derivatives/change fields for updating state."
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
    strain = PYB11property("Field<DIM, SymTensor>&", returnpolicy="reference_internal")
    effectiveStrain = PYB11property("const Field<DIM, SymTensor>&", returnpolicy="reference_internal")
    DdamageDt = PYB11property("const Field<DIM, Scalar>&", returnpolicy="reference_internal")
    newEffectiveDamage = PYB11property("const Field<DIM, SymTensor>&", returnpolicy="reference_internal")
    newDamageGradient = PYB11property("const Field<DIM, Vector>&", returnpolicy="reference_internal")

    strainAlgorithm = PYB11property("TensorStrainAlgorithm")
    effectiveDamageAlgorithm = PYB11property("EffectiveDamageAlgorithm")
    useDamageGradient = PYB11property("bool", "useDamageGradient", "useDamageGradient",
                                      doc="Flag to determine if we compute the gradient of the damage at the start of a timestep.")
    damageInCompression = PYB11property("bool", "damageInCompression", "damageInCompression",
                                        doc="Flag to determine if damage in compression is allowed.")
    criticalDamageThreshold = PYB11property("double", "criticalDamageThreshold", "criticalDamageThreshold",
                                            doc="The critical damage threshold for not setting the time step.")

#-------------------------------------------------------------------------------
# Add the restart methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, TensorDamageModel)
