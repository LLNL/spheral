#-------------------------------------------------------------------------------
# ConstantAcceleration base class
#-------------------------------------------------------------------------------
from PYB11Generator import *
from GenericBodyForce import *

@PYB11template("Dimension")
class ConstantAcceleration(GenericBodyForce):

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename Physics<%(Dimension)s>::TimeStepType TimeStepType;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               a0 = "const Vector",
               nodeList = "const NodeList<%(Dimension)s>&",
               indices = "const std::vector<int>&"):
        "ConstantAcceleration constructor"

    def pyinit1(self,
                a0 = "const Vector",
                nodeList = "const NodeList<%(Dimension)s>&"):
        "ConstantAcceleration constructor"

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
        "Increment the derivatives."
        return "void"

    @PYB11virtual
    @PYB11const
    def dt(dataBase = "const DataBase<%(Dimension)s>&", 
           state = "const State<%(Dimension)s>&",
           derivs = "const StateDerivatives<%(Dimension)s>&",
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
    nodeList = PYB11property("const NodeList<%(Dimension)s>&", "nodeList", returnpolicy="reference_internal", doc="The NodeList this ConstantAcceleration applies to")
    flags = PYB11property("const Field<%(Dimension)s, int>&", "flags", returnpolicy="reference_internal", doc="The node flags")
