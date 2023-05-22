#-------------------------------------------------------------------------------
# TreeGravity base class
#-------------------------------------------------------------------------------
from PYB11Generator import *
from RestartMethods import *
from GenericBodyForce import *

@PYB11template("Dimension")
class TreeGravity(GenericBodyForce):

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    typedef typename %(Dimension)s::ThirdRankTensor ThirdRankTensor;
    typedef typename Physics<%(Dimension)s>::TimeStepType TimeStepType;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               G = "const double",
               softeningLength = "const double",
               opening = "const double",
               ftimestep = "const double",
               timeStepChoice = "GravityTimeStepType"):
        "TreeGravity constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def registerState(self,
                      dataBase = "DataBase<%(Dimension)s>&",
                      state = "State<%(Dimension)s>&"):
        "Register the state you want carried around (and potentially evolved), as well as the policies for such evolution."
        return "void"

    @PYB11virtual
    @PYB11const
    def evaluateDerivatives(self,
                            time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<%(Dimension)s>&",
                            state = "const State<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Increment the derivatives."
        return "void"

    @PYB11virtual
    @PYB11const
    def dt(dataBase = "const DataBase<%(Dimension)s>&", 
           state = "const State<%(Dimension)s>&",
           derivs = "const StateDerivatives<%(Dimension)s>&",
           currentTime = "const Scalar"):
        "Vote on a time step."
        return "TimeStepType"

    @PYB11virtual
    def initializeProblemStartup(self,
                                 dataBase = "DataBase<%(Dimension)s>&"):
        "An optional hook to initialize once when the problem is starting up."
        return "void"

    @PYB11virtual
    def initialize(self,
                   time = "const Scalar", 
                   dt = "const Scalar",
                   dataBase = "const DataBase<%(Dimension)s>&", 
                   state = "State<%(Dimension)s>&",
                   derivs = "StateDerivatives<%(Dimension)s>&"):
        "Some packages might want a hook to do some initializations before the evaluateDerivatives() method is called."
        return "void"

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

    @PYB11const
    def dumpTree(self, globalTree="const bool"):
        "Return a dump of the tree structure as a string."
        return "std::string"

    @PYB11const
    def dumpTreeStatistics(self, globalTree="const bool"):
        "Return a string describing the overall statistics of the tree."
        return "std::string"

    #...........................................................................
    # Properties
    potential = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "potential", returnpolicy="reference_internal", doc="The last computed potential")
    G = PYB11property("double", "G", doc="The gravitational constant")
    opening = PYB11property("double", "opening", "opening", doc="The opening angle threshold when we shift to tree cell approximations.")
    softeningLength = PYB11property("Scalar", "softeningLength", "softeningLength", doc="The Plummer softening scale")
    ftimestep = PYB11property("double", "ftimestep", "ftimestep", doc="The current time step scaling factor.")
    timeStepChoice = PYB11property("GravityTimeStepType", "timeStepChoice", "timeStepChoice", doc="The algorithmic choice for setting the time step.")
    xmin = PYB11property("Vector", "xmin", doc="The lower left corner of the computational cube that was last used.")
    xmax = PYB11property("Vector", "xmax", doc="The upper right corner of the computational cube that was last used.")

#-------------------------------------------------------------------------------
# Inject restart methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, TreeGravity)
