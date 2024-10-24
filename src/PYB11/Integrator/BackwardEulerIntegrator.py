#-------------------------------------------------------------------------------
# BackwardEulerIntegrator
#-------------------------------------------------------------------------------
from PYB11Generator import *
from IntegratorAbstractMethods import *
from ImplicitIntegrator import *

@PYB11template("Dimension")
@PYB11cppname("BackwardEuler")
class BackwardEulerIntegrator(ImplicitIntegrator):
    "First-order in time implicit (backward Euler) integration scheme"

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
               dataBase = "DataBase<%(Dimension)s>&"):
        """Construct a backward Euler itegrator.
Note: beta is the blending of the (n+1) and (n) time derivatives, so adjusting this parameter
makes this not just backward Euler:
    beta = 1.0 : backward Euler
    beta = 0.0 : forward Euler
    beta = 0.5 : Crank-Nicholson"""
        return

    @PYB11implementation("[](DataBase<%(Dimension)s>& db, py::list pypkgs) { return std::make_unique<BackwardEuler<%(Dimension)s>>(db, PYB11utils::from_list<Physics<%(Dimension)s>*>(pypkgs)); }")
    def pyinit1(self):
        """Construct a backward Euler itegrator.
Note: beta is the blending of the (n+1) and (n) time derivatives, so adjusting this parameter
makes this not just backward Euler:
    beta = 1.0 : backward Euler
    beta = 0.0 : forward Euler
    beta = 0.5 : Crank-Nicholson"""

    @PYB11implementation("[](DataBase<%(Dimension)s>& db, py::list pypkgs, Scalar beta) { return std::make_unique<BackwardEuler<%(Dimension)s>>(db, PYB11utils::from_list<Physics<%(Dimension)s>*>(pypkgs), beta); }")
    def pyinit2(self):
        """Construct a backward Euler itegrator.
Note: beta is the blending of the (n+1) and (n) time derivatives, so adjusting this parameter
makes this not just backward Euler:
    beta = 1.0 : backward Euler
    beta = 0.0 : forward Euler
    beta = 0.5 : Crank-Nicholson"""

    @PYB11implementation("[](DataBase<%(Dimension)s>& db, py::list pypkgs, Scalar beta, Scalar tol) { return std::make_unique<BackwardEuler<%(Dimension)s>>(db, PYB11utils::from_list<Physics<%(Dimension)s>*>(pypkgs), beta, tol); }")
    def pyinit3(self):
        """Construct a backward Euler itegrator.
Note: beta is the blending of the (n+1) and (n) time derivatives, so adjusting this parameter
makes this not just backward Euler:
    beta = 1.0 : backward Euler
    beta = 0.0 : forward Euler
    beta = 0.5 : Crank-Nicholson"""

    @PYB11implementation("[](DataBase<%(Dimension)s>& db, py::list pypkgs, Scalar beta, Scalar tol, size_t maxIterations) { return std::make_unique<BackwardEuler<%(Dimension)s>>(db, PYB11utils::from_list<Physics<%(Dimension)s>*>(pypkgs), beta, tol, maxIterations); }")
    def pyinit4(self):
        """Construct a backward Euler itegrator.
Note: beta is the blending of the (n+1) and (n) time derivatives, so adjusting this parameter
makes this not just backward Euler:
    beta = 1.0 : backward Euler
    beta = 0.0 : forward Euler
    beta = 0.5 : Crank-Nicholson"""

#     @PYB11implementation("""
#         [](DataBase<%(Dimension)s>& db,
#            py::list pypkgs = py::list(),
#            Scalar beta = 1.0,
#            Scalar tol = 1.0e-6,
#            size_t maxIterations = 10u) {
#             std::vector<Physics<%(Dimension)s>*> pkgs = PYB11utils::from_list<Physics<%(Dimension)s>*>(pypkgs);
#             // for (auto x: pypkgs) pkgs.push_back(x.cast<Physics<%(Dimension)s>*>());
#             return std::make_unique<BackwardEuler<%(Dimension)s>>(db, pkgs, beta, tol, maxIterations);
#         }""")
#     def pyinit(self,
#                dataBase = "DataBase<%(Dimension)s>&",
#                physicsPackages = ("py::list", "py::list()"),
#                beta = ("Scalar", "1.0"),
#                tol = ("Scalar", "1.0e-6"),
#                maxIterations = ("size_t", "10u")):
#         """Construct a backward Euler itegrator.
# Note: beta is the blending of the (n+1) and (n) time derivatives, so adjusting this parameter
# makes this not just backward Euler:
#     beta = 1.0 : backward Euler
#     beta = 0.0 : forward Euler
#     beta = 0.5 : Crank-Nicholson"""

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
    beta = PYB11property("Scalar", "beta", "beta", doc="The blend of (n+1) and (n) derivative states for advancement")
    tol = PYB11property("Scalar", "tol", "tol", doc="The tolerance to test for convergence of a step")
    maxIterations = PYB11property("size_t", "maxIterations", "maxIterations", doc="The maximum allowed iterations to try for advancing a step")

#-------------------------------------------------------------------------------
# Inject other interfaces
#-------------------------------------------------------------------------------
PYB11inject(IntegratorAbstractMethods, BackwardEulerIntegrator, pure_virtual=False, virtual=True)
