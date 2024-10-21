#-------------------------------------------------------------------------------
# ImplicitIntegrator base class
#-------------------------------------------------------------------------------
from PYB11Generator import *
from IntegratorAbstractMethods import *
from Integrator import *

@PYB11template("Dimension")
@PYB11module("SpheralIntegrator")
class ImplicitIntegrator(Integrator):
    "Base class for all Spheral implicit in time integration algorithms"

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    typedef typename %(Dimension)s::ThirdRankTensor ThirdRankTensor;
"""

    #...........................................................................
    # Constructors
    def pyinit(self):
        "Construct an ImplicitIntegrator"

    def pyinit1(self, dataBase = "DataBase<%(Dimension)s>&"):
        "Construct an ImplicitIntegrator with a DataBase"

    def pyinit2(self,
                dataBase = "DataBase<%(Dimension)s>&",
                physicsPackages = "const std::vector<Physics<%(Dimension)s>*>&"):
        "Construct an ImplicitIntegrator with a DataBase and physics packages"

#-------------------------------------------------------------------------------
# Inject other interfaces
#-------------------------------------------------------------------------------
PYB11inject(IntegratorAbstractMethods, ImplicitIntegrator, virtual=False, pure_virtual=True)
PYB11inject(RestartMethods, ImplicitIntegrator)
