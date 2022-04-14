#-------------------------------------------------------------------------------
# LimitedMonaghanGingoldViscosity
#-------------------------------------------------------------------------------
from PYB11Generator import *
from MonaghanGingoldViscosity import *
from ArtificialViscosityAbstractMethods import *

@PYB11template("Dimension")
class LimitedMonaghanGingoldViscosity(MonaghanGingoldViscosity):

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    typedef typename %(Dimension)s::ThirdRankTensor ThirdRankTensor;
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
        "LimitedMonaghanGingoldViscosity constructor"


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
PYB11inject(ArtificialViscosityAbstractMethods, LimitedMonaghanGingoldViscosity, virtual=True, pure_virtual=False)
