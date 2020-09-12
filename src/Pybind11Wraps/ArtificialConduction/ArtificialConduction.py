#-------------------------------------------------------------------------------
# ArtificialConduction
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Physics import *

@PYB11template("Dimension")
class ArtificialConduction(Physics):
    "ArtificialConduction -- Artificial smoothing of energy discontinuities"

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    typedef typename %(Dimension)s::ThirdRankTensor ThirdRankTensor;
    typedef typename Physics<%(Dimension)s>::TimeStepType TimeStepType;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               W = "const TableKernel<%(Dimension)s>&",
               alphaArCond = "const Scalar",
               ACcorrectionOrder = ("const RKOrder", "RKOrder::LinearOrder")):
        "Constructor"
            
        #...........................................................................
        # Virtual methods
        @PYB11virtual
        def initializeProblemStartup(self,
                                     dataBase = "DataBase<%(Dimension)s>&"):
            "Do any required one-time initializations on problem start up."
            return "void"
            
        #...........................................................................
        # Properties
        ACcorrectionOrder = PYB11property("RKOrder", "ACcorrectionOrder", "ACcorrectionOrder",
                                          doc="The correction order")

#-------------------------------------------------------------------------------
# Inject physics interface
#-------------------------------------------------------------------------------
PYB11inject(PhysicsAbstractMethods, ArtificialConduction, virtual=True, pure_virtual=False)
