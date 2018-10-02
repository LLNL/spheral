#-------------------------------------------------------------------------------
# CRKSPHMonaghanGingoldViscosity
#-------------------------------------------------------------------------------
from PYB11Generator import *
from MonaghanGingoldViscosity import *
from ArtificialViscosityAbstractMethods import *

@PYB11template("Dimension")
class CRKSPHMonaghanGingoldViscosity(MonaghanGingoldViscosity):

    typedefs = """
    typedef %(Dimension)s DIM;
    typedef typename DIM::Scalar Scalar;
    typedef typename DIM::Vector Vector;
    typedef typename DIM::Tensor Tensor;
    typedef typename DIM::SymTensor SymTensor;
    typedef typename DIM::ThirdRankTensor ThirdRankTensor;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               Clinear = "const Scalar",
               Cquadratic = "const Scalar",
               linearInExpansion = ("bool", "false"),
               quadraticInExpansion = ("bool", "false"),
               etaCritFrac = ("double", "1.0"),
               etaFoldFrac = ("double", "0.2")):
        "CRKSPHMonaghanGingoldViscosity constructor"


    #...........................................................................
    # Virtual methods
    @PYB11virtual
    @PYB11const
    def label(self):
        return "std::string"

    #...........................................................................
    # Properties
    etaCritFrac = PYB11property("double", "etaCritFrac", "etaCritFrac",
                                doc="Critical eta, below which viscosity becomes active regardless of limiter")
    etaFoldFrac = PYB11property("double", "etaFoldFrac", "etaFoldFrac",
                                doc="Smoothness in eta phasing in viscosity around etaCritFrac")
    
#-------------------------------------------------------------------------------
# Inject abstract interface
#-------------------------------------------------------------------------------
PYB11inject(ArtificialViscosityAbstractMethods, CRKSPHMonaghanGingoldViscosity, virtual=True, pure_virtual=False)
