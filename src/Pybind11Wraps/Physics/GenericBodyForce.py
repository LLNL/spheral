#-------------------------------------------------------------------------------
# GenericBodyForce base class
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Physics import *

@PYB11template("Dimension")
@PYB11module("SpheralPhysics")
class GenericBodyForce(Physics):

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename Physics<%(Dimension)s>::TimeStepType TimeStepType;
"""

    #...........................................................................
    # Constructors
    def pyinit(self):
        "GenericBodyForce constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def registerState(self,
                      dataBase = "DataBase<%(Dimension)s>&", 
                      state = "State<%(Dimension)s>&"):
        "Default state registration for an acceleration source"
        return "void"

    @PYB11virtual
    def registerDerivatives(self,
                            dataBase = "DataBase<%(Dimension)s>&", 
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Default state derivative registration for an acceleration source"
        return "void"

    @PYB11pure_virtual
    @PYB11const
    def dt(dataBase = "const DataBase<%(Dimension)s>&", 
           state = "const State<%(Dimension)s>&",
           derivs = "const StateDerivatives<%(Dimension)s>&",
           currentTime = "const Scalar"):
        "Vote on a time step."
        return "TimeStepType"
    #...........................................................................
    # Attributes
    DxDt = PYB11property("const FieldList<%(Dimension)s, Vector>&", "DxDt", doc="Time derivative for position")
    DvDt = PYB11property("const FieldList<%(Dimension)s, Vector>&", "DvDt", doc="Time derivative for velocity")
