#-------------------------------------------------------------------------------
# SynchronousRK1Integrator
#-------------------------------------------------------------------------------
from PYB11Generator import *
from IntegratorAbstractMethods import *
from Integrator import *

@PYB11template("Dimension")
@PYB11cppname("SynchronousRK1")
class SynchronousRK1Integrator(Integrator):
    "First-order in time explicit (forward Euler) integration scheme"

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
        "Construct an itegrator"

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
PYB11inject(IntegratorAbstractMethods, SynchronousRK1Integrator, pure_virtual=False, virtual=True)
