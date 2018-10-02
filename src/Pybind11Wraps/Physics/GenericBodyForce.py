#-------------------------------------------------------------------------------
# GenericBodyForce base class
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Physics import *

@PYB11template("Dimension")
class GenericBodyForce(Physics):

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
    def pyinit(self):
        "GenericBodyForce constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def registerState(self,
                      dataBase = "DataBase<DIM>&", 
                      state = "State<DIM>&"):
        "Default state registration for an acceleration source"
        return "void"

    @PYB11virtual
    def registerDerivatives(self,
                            dataBase = "DataBase<DIM>&", 
                            derivs = "StateDerivatives<DIM>&"):
        "Default state derivative registration for an acceleration source"
        return "void"

    #...........................................................................
    # Attributes
    DxDt = PYB11property("const FieldList<DIM, Vector>&", "DxDt", doc="Time derivative for position")
    DvDt = PYB11property("const FieldList<DIM, Vector>&", "DvDt", doc="Time derivative for velocity")
