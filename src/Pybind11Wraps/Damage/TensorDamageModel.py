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

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    typedef typename Physics<%(Dimension)s>::TimeStepType TimeStepType;

    typedef Field<%(Dimension)s, std::vector<double> > FlawStorageType;
"""

    def pyinit(self,
               nodeList = "SolidNodeList<%(Dimension)s>&",
               strainAlgorithm = "const TensorStrainAlgorithm",
               effectiveDamageAlgorithm = "const EffectiveDamageAlgorithm",
               useDamageGradient = "const bool",
               kernel = "const TableKernel<%(Dimension)s>&",
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
                            dataBase = "const DataBase<%(Dimension)s>&",
                            state = "const State<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Compute the derivatives."
        return "void"

    @PYB11virtual
    @PYB11const
    def dt(self,
           dataBase = "const DataBase<%(Dimension)s>&",
           state = "const State<%(Dimension)s>&",
           derivs = "const StateDerivatives<%(Dimension)s>&",
           currentTime = "const Scalar"):
        "Vote on a time step."
        return "TimeStepType"

    @PYB11virtual
    def registerState(self,
                      dataBase = "DataBase<%(Dimension)s>&",
                      state = "State<%(Dimension)s>&"):
        "Register our state."
        return "void"

    @PYB11virtual
    def registerDerivatives(self,
                            dataBase = "DataBase<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Register the derivatives/change fields for updating state."
        return "void"

    @PYB11virtual
    def applyGhostBoundaries(self,
                             state = "State<%(Dimension)s>&",
                             derivs = "StateDerivatives<%(Dimension)s>&"):
        "Apply boundary conditions to the physics specific fields."
        return "void"

    @PYB11virtual
    def enforceBoundaries(self,
                          state = "State<%(Dimension)s>&",
                          derivs = "StateDerivatives<%(Dimension)s>&"):
        "Enforce boundary conditions for the physics specific fields."
        return "void"

    #...........................................................................
    # Properties
    strain = PYB11property("const Field<%(Dimension)s, SymTensor>&", returnpolicy="reference_internal")
    effectiveStrain = PYB11property("const Field<%(Dimension)s, SymTensor>&", returnpolicy="reference_internal")
    DdamageDt = PYB11property("const Field<%(Dimension)s, Scalar>&", returnpolicy="reference_internal")
    newEffectiveDamage = PYB11property("const Field<%(Dimension)s, SymTensor>&", returnpolicy="reference_internal")
    newDamageGradient = PYB11property("const Field<%(Dimension)s, Vector>&", returnpolicy="reference_internal")

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
