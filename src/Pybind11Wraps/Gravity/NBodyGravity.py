#-------------------------------------------------------------------------------
# NBodyGravity base class
#-------------------------------------------------------------------------------
from PYB11Generator import *
from GenericBodyForce import *

@PYB11template("Dimension")
class NBodyGravity(GenericBodyForce):

    typedefs = """
    typedef %(Dimension)s DIM;
    typedef typename DIM::Scalar Scalar;
    typedef typename DIM::Vector Vector;
    typedef typename DIM::Tensor Tensor;
    typedef typename DIM::SymTensor SymTensor;
    typedef typename DIM::ThirdRankTensor ThirdRankTensor;
    typedef typename Physics<DIM>::TimeStepType TimeStepType;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               plummerSofteningLength = "const double",
               maxDeltaVelocity = "const double",
               G = "const double",
               compatibleVelocityUpdate = "const bool"):
        "NBodyGravity constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def registerState(self,
                      dataBase = "DataBase<DIM>&",
                      state = "State<DIM>&"):
        "Register the state you want carried around (and potentially evolved), as well as the policies for such evolution."
        return "void"

    @PYB11virtual
    @PYB11const
    def evaluateDerivatives(self,
                            time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<DIM>&",
                            state = "const State<DIM>&",
                            derivs = "StateDerivatives<DIM>&"):
        "Increment the derivatives."
        return "void"

    @PYB11virtual
    @PYB11const
    def dt(dataBase = "const DataBase<DIM>&", 
           state = "const State<DIM>&",
           derivs = "const StateDerivatives<DIM>&",
           currentTime = "const Scalar"):
        "Vote on a time step."
        return "TimeStepType"

    @PYB11virtual
    def initializeProblemStartup(self,
                                 dataBase = "DataBase<DIM>&"):
        "An optional hook to initialize once when the problem is starting up."
        return "void"

    @PYB11virtual
    def preStepInitialize(self,
                          dataBase = "const DataBase<DIM>&", 
                          state = "State<DIM>&",
                          derivs = "StateDerivatives<DIM>&"):
        "Optional hook to be called at the beginning of a time step."
        return "void"

    @PYB11virtual
    def finalize(self,
                 time = "const Scalar", 
                 dt = "const Scalar",
                 dataBase = "DataBase<DIM>&", 
                 state = "State<DIM>&",
                 derivs = "StateDerivatives<DIM>&"):
        "Similarly packages might want a hook to do some post-step finalizations.  Really we should rename this post-step finalize."
        return "void"

    @PYB11virtual
    @PYB11const
    def label(self):
        "It's useful to have labels for Physics packages.  We'll require this to have the same signature as the restart label."
        return "std::string"

    @PYB11virtual
    @PYB11const
    def requireConnectivity(self):
        "Some physics does not require the connectivity be constructed."
        return "bool"

    @PYB11virtual
    @PYB11const
    def extraEnergy(self):
        "Many physics packages will have their own representations of energy in the system (gravitational potential energy, radiative losses, etc.)"
        return "Scalar"

    @PYB11virtual
    @PYB11const
    def valid(self):
        return "bool"

    #...........................................................................
    # Properties
    potential = PYB11property("const FieldList<DIM, Scalar>&", "potential", returnpolicy="reference_internal", doc="The last computed potential")
    G = PYB11property("double", "G", doc="The gravitational constant")
    softeningLength = PYB11property("Scalar", "softeningLength", "softeningLength", doc="The Plummer softening scale")
    compatibleVelocityUpdate = PYB11property("bool", "compatibleVelocityUpdate", "compatibleVelocityUpdate", doc="Experimental compatible update for velocity")
