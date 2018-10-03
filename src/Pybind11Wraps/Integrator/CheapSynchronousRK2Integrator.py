#-------------------------------------------------------------------------------
# CheapSynchronousRK2Integrator
#-------------------------------------------------------------------------------
from PYB11Generator import *
from IntegratorAbstractMethods import *
from Integrator import *

@PYB11template("Dimension")
@PYB11cppname("CheapSynchronousRK2")
class CheapSynchronousRK2Integrator(Integrator):
    """Second-order in time explicit Runge-Kutta time integration scheme
This method cheats and reuses the previous step mid-timestep time derivatives
to predict the state at the middle of the current step, rather than evaluate
the derivatives at the beginning of this step.  This allows us to only do a single
call to evaluate derivatives per timestep, maintaining formal second-order but 
sacrificing some accuracy vs. the the true RK2 algorithm."""

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
PYB11inject(IntegratorAbstractMethods, CheapSynchronousRK2Integrator, pure_virtual=False, virtual=True)
