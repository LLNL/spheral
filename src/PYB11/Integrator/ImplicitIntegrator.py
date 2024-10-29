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
    def pyinit(self, dataBase = "DataBase<%(Dimension)s>&"):
        "Construct an ImplicitIntegrator with a DataBase"

    def pyinit1(self,
                dataBase = "DataBase<%(Dimension)s>&",
                physicsPackages = ("const std::vector<Physics<%(Dimension)s>*>&", "std::vector<Physics<%(Dimension)s>*>()"),
                tol = ("Scalar", "1.0e-6")):
        "Construct an ImplicitIntegrator with a DataBase, physics packages, and tolerance"

    #...........................................................................
    # Methods
    @PYB11virtual
    @PYB11const
    def computeResiduals(self,
                         state1 = "const State<%(Dimension)s>&",
                         state0 = "const State<%(Dimension)s>&"):
        "Compute the maximum residual difference across all physics packages between two States"
        return "Scalar"

    #...........................................................................
    # Properties
    convergenceTolerance = PYB11property("Scalar", "convergenceTolerance", "convergenceTolerance", doc="Tolerance for convergence in residuals during integration")

#-------------------------------------------------------------------------------
# Inject other interfaces
#-------------------------------------------------------------------------------
PYB11inject(IntegratorAbstractMethods, ImplicitIntegrator, virtual=False, pure_virtual=True)
