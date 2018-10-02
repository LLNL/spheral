#-------------------------------------------------------------------------------
# ConstantAcceleration base class
#-------------------------------------------------------------------------------
from PYB11Generator import *
from GenericBodyForce import *

@PYB11template("Dimension")
class ConstantAcceleration(GenericBodyForce):

    typedefs = """
    typedef %(Dimension)s DIM;
    typedef typename DIM::Scalar Scalar;
    typedef typename DIM::Vector Vector;
    typedef typename DIM::Tensor Tensor;
    typedef typename DIM::SymTensor SymTensor;
    typedef typename DIM::ThirdRankTensor ThirdRankTensor;
    typedef typename Physics<DIM>::TimeStepType TimeStepType;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               a0 = "const Vector",
               nodeList = "const NodeList<DIM>&",
               indices = "const std::vector<int>&"):
        "ConstantAcceleration constructor"

    def pyinit1(self,
                a0 = "const Vector",
                nodeList = "const NodeList<DIM>&"):
        "ConstantAcceleration constructor"

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
        "Increment the derivatives."
        return "void"

    @PYB11virtual
    @PYB11const
    def dt(dataBase = "const DataBase<DIM>&", 
           state = "const State<DIM>&",
           derivs = "const StateDerivatives<DIM>&",
           currentTime = "const Scalar"):
        "Vote on a time step."
        return "TimeStepType"

    @PYB11virtual
    @PYB11const
    def label(self):
        "It's useful to have labels for Physics packages.  We'll require this to have the same signature as the restart label."
        return "std::string"

    #...........................................................................
    # Properties
    a0 = PYB11property("Vector", "a0", doc="The fixed acceleration to apply")
    nodeList = PYB11property("const NodeList<DIM>&", "nodeList", returnpolicy="reference_internal", doc="The NodeList this ConstantAcceleration applies to")
    flags = PYB11property("const Field<DIM, int>&", "flags", returnpolicy="reference_internal", doc="The node flags")
