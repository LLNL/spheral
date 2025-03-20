#-------------------------------------------------------------------------------
# CrankNicolsonIntegrator
#-------------------------------------------------------------------------------
from PYB11Generator import *
from IntegratorAbstractMethods import *
from ImplicitIntegrator import *

@PYB11template("Dimension")
@PYB11cppname("CrankNicolson")
class CrankNicolsonIntegrator(ImplicitIntegrator):
    "Second-order in time implicit Crank-Nicolson integration scheme"

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using Tensor = typename %(Dimension)s::Tensor;
    using SymTensor = typename %(Dimension)s::SymTensor;
    using ThirdRankTensor = typename %(Dimension)s::ThirdRankTensor;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               dataBase = "DataBase<%(Dimension)s>&",
               packages = ("vector<Physics<%(Dimension)s>*>", "vector<Physics<%(Dimension)s>*>()"),
               alpha = ("Scalar", "0.5"),
               tol = ("Scalar", "1.0e-8"),
               maxIterations = ("size_t", "100u")):
        """Construct a Crank-Nicolson itegrator.
Note: alhpa is the blending of the (n+1) and (n) time solutions in the fixed-point iteration."""
        return

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def step(self,
             maxTime = "Scalar",
             state = "State<%(Dimension)s>&",
             derivs = "StateDerivatives<%(Dimension)s>&"):
        "Take a step"
        return "bool"

    @PYB11virtual
    @PYB11const
    def label(self):
        return "std::string"

    #...........................................................................
    # Properties
    alpha = PYB11property("Scalar", "alpha", "alpha", doc="The blend of (n+1) and (n) state during fixed-point iteration")
    maxIterations = PYB11property("size_t", "maxIterations", "maxIterations", doc="The maximum allowed iterations to try for advancing a step")
    numExplicitSteps = PYB11property("size_t", "numExplicitSteps")
    numImplicitSteps = PYB11property("size_t", "numImplicitSteps")

#-------------------------------------------------------------------------------
# Inject other interfaces
#-------------------------------------------------------------------------------
PYB11inject(IntegratorAbstractMethods, CrankNicolsonIntegrator, pure_virtual=False, virtual=True)
