#-------------------------------------------------------------------------------
# ContactModel base class
#-------------------------------------------------------------------------------
from PYB11Generator import *
from ContactModelBase import *

@PYB11template("Dimension")
@PYB11module("SpheralDEM")
class DampedLinearSpring(ContactModelBase):

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::SymTensor SymTensor;
"""

    #...........................................................................
    # Constructors
    def pyinit(dataBase = "const DataBase<%(Dimension)s>&",
               YoungsModulus = "const Scalar",
               RestitutionCoefficient = "const Scalar"):
        "Damped linear spring contact model constructor"

    #...........................................................................
    # Methods
    @PYB11virtual
    @PYB11const
    def timeStep(dataBase = "const DataBase<%(Dimension)s>&",
                 state =    "const State<%(Dimension)s>&",
                 derivs =   "const StateDerivatives<%(Dimension)s>&",
                 time =     "const Scalar"):
        """Time step for the damped linear spring contact model."""
        return "Scalar"

    @PYB11virtual
    @PYB11const
    def evaluateDerivatives(time =     "const Scalar",
                            dt =       "const Scalar",
                            dataBase = "const DataBase<%(Dimension)s>&",
                            state =    "const State<%(Dimension)s>&",
                            derivs =   "StateDerivatives<%(Dimension)s>&"):
        """get derivs from the damped linear spring contact model."""
        return "void"

    #...........................................................................
    # Attributes
    normalSpringConstant = PYB11property("Scalar", "normalSpringConstant", "normalSpringConstant", doc="linear spring constant for compression")
    resitutionCoefficient = PYB11property("Scalar", "restitutionCoefficient", "restitutionCoefficient", doc="ratio of incoming to outgoing velocity used to calculate damping coefficient")
    beta = PYB11property("Scalar", "beta", "beta", doc="a damping parameter")
    timeStep = PYB11property("Scalar", "timeStep", "timeStep", doc="the time step")
