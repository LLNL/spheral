#-------------------------------------------------------------------------------
# VerletIntegrator
#-------------------------------------------------------------------------------
from PYB11Generator import *
from IntegratorAbstractMethods import *
from Integrator import *

@PYB11template("Dimension")
@PYB11cppname("Verlet")
class VerletIntegrator(Integrator):
    """Second-order in time explicit Verlet time integration scheme
This method is symplectic in the absence of dissipation."""

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
    def pyinit(self):
        "Construct an itegrator"

    def pyinit1(self, dataBase = "DataBase<DIM>&"):
        "Construct an integrator with a DataBase"

    def pyinit1(self,
                dataBase = "DataBase<DIM>&",
                physicsPackages = "const std::vector<Physics<DIM>*>&"):
        "Construct an integrator with a DataBase and physics packages"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    @PYB11const
    def label(self):
        return "std::string"

#-------------------------------------------------------------------------------
# Inject other interfaces
#-------------------------------------------------------------------------------
PYB11inject(IntegratorAbstractMethods, VerletIntegrator, pure_virtual=False, virtual=True)
