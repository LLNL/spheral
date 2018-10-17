#-------------------------------------------------------------------------------
# ArtificialConduction
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Physics import *

@PYB11template("Dimension")
class ArtificialConduction(Physics):
    "ArtificialConduction -- Artificial smoothing of energy discontinuities"

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
               W = "const TableKernel<DIM>&",
               alphaArCond = "const Scalar",
               ACcorrectionOrder = ("const CRKOrder", "CRKOrder::LinearOrder")):
        "Constructor"
            
        #...........................................................................
        # Virtual methods
        @PYB11virtual
        def initializeProblemStartup(self,
                                     dataBase = "DataBase<DIM>&"):
            "Do any required one-time initializations on problem start up."
            return "void"
            
        #...........................................................................
        # Properties
        ACcorrectionOrder = PYB11property("CRKOrder", "ACcorrectionOrder", "ACcorrectionOrder",
                                          doc="The correction order")

#-------------------------------------------------------------------------------
# Inject physics interface
#-------------------------------------------------------------------------------
PYB11inject(PhysicsAbstractMethods, ArtificialConduction, virtual=True, pure_virtual=False)
