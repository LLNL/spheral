#-------------------------------------------------------------------------------
# StateDerivatives
#-------------------------------------------------------------------------------
from PYB11Generator import *
from StateBase import *

@PYB11template("Dimension")
class StateDerivatives(StateBase):

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    typedef typename StateBase<%(Dimension)s>::KeyType KeyType;
    typedef typename StateBase<%(Dimension)s>::FieldName FieldName;
    typedef typename StateBase<%(Dimension)s>::MeshPtr MeshPtr;
    typedef typename StateDerivatives<%(Dimension)s>::PackageList PackageList;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               dataBase = "DataBase<%(Dimension)s>&",
               physicsPackages = "PackageList&"):
        "Construct using the Physics::registerDerivatives methods of the packages"

    #...........................................................................
    # Operators
    def __eq__(self):
        return

    #...........................................................................
    # Methods
    def initializeNodePairInformation(self):
        return "void"

    @PYB11const
    def calculatedNodePairsSymmetric(self):
        return "bool"

    def Zero(self):
        "Set all derivative Fields to zero"
        return "void"
