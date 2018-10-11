#-------------------------------------------------------------------------------
# JohnsonCookDamage
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Physics import *
from RestartMethods import *

@PYB11template("Dimension")
class JohnsonCookDamage(Physics):
    "JohnsonCookDamage -- an implementation of a Johnson-Cook damage law."

    typedefs = """
    typedef %(Dimension)s DIM;
    typedef typename DIM::Scalar Scalar;
    typedef typename DIM::Vector Vector;
    typedef typename DIM::Tensor Tensor;
    typedef typename DIM::SymTensor SymTensor;
    typedef typename Physics<DIM>::TimeStepType TimeStepType;
"""

    def pyinit(self,
               nodeList = "SolidNodeList<DIM>&",
               D1 = "const Field<DIM, Scalar>&",
               D2 = "const Field<DIM, Scalar>&",
               D3 = "const double",
               D4 = "const double",
               D5 = "const double",
               epsilondot0 = "const double",
               Tcrit = "const double",
               sigmamax = "const double",
               efailmin = "const double"):
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
    nodeList = PYB11property("const SolidNodeList<DIM>&", returnpolicy="reference_internal")
    failureStrain = PYB11property("const Field<DIM, Scalar>&", returnpolicy="reference_internal")
    meltSpecificEnergy = PYB11property("const Field<DIM, Scalar>&", returnpolicy="reference_internal")
    newEffectiveDamage = PYB11property("const Field<DIM, SymTensor>&", returnpolicy="reference_internal")
    D1 = PYB11property("const Field<DIM, Scalar>&", returnpolicy="reference_internal")
    D2 = PYB11property("const Field<DIM, Scalar>&", returnpolicy="reference_internal")
    D3 = PYB11property("double")
    D4 = PYB11property("double")
    D5 = PYB11property("double")
    epsilondot0 = PYB11property("double")
    Tcrit = PYB11property("double")
    sigmamax = PYB11property("double")
    efailmin = PYB11property("double")

#-------------------------------------------------------------------------------
# Add the restart methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, JohnsonCookDamage)
