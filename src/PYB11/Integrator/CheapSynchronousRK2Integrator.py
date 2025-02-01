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

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    typedef typename %(Dimension)s::ThirdRankTensor ThirdRankTensor;
"""

    #...........................................................................
    # Constructors
    def pyinit1(self, dataBase = "DataBase<%(Dimension)s>&"):
        "Construct an integrator with a DataBase"

    def pyinit2(self,
                dataBase = "DataBase<%(Dimension)s>&",
                physicsPackages = "const std::vector<Physics<%(Dimension)s>*>&"):
        "Construct an integrator with a DataBase and physics packages"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    @PYB11pycppname("step")
    def step1(self, maxTime="Scalar"):
        "Take a step"
        return "bool"

    @PYB11virtual
    @PYB11const
    def label(self):
        return "std::string"

#-------------------------------------------------------------------------------
# Inject other interfaces
#-------------------------------------------------------------------------------
PYB11inject(IntegratorAbstractMethods, CheapSynchronousRK2Integrator, pure_virtual=False, virtual=True)
