#-------------------------------------------------------------------------------
# LimitedMonaghanGingoldViscosity
#-------------------------------------------------------------------------------
from PYB11Generator import *
from MonaghanGingoldViscosity import *
from ArtificialViscosityAbstractMethods import *

@PYB11template("Dimension")
@PYB11template_dict({"QPiType": "typename %(Dimension)s::Scalar"})
class LimitedMonaghanGingoldViscosity(MonaghanGingoldViscosity):

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using Tensor = typename %(Dimension)s::Tensor;
    using SymTensor = typename %(Dimension)s::SymTensor;
    using ThirdRankTensor = typename %(Dimension)s::ThirdRankTensor;
    using ReturnType = %(QPiType)s;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               Clinear = "const Scalar",
               Cquadratic = "const Scalar",
               linearInExpansion = "bool",
               quadraticInExpansion = "bool",
               etaCritFrac = "double",
               etaFoldFrac = "double",
               kernel = "const TableKernel<%(Dimension)s>&"):
        "LimitedMonaghanGingoldViscosity constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    @PYB11const
    def requireVelocityGradient(self):
        "We need the velocity gradient and set this to true"
        return "bool"

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
