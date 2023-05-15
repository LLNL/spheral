#-------------------------------------------------------------------------------
# PointPotential base class
#-------------------------------------------------------------------------------
from PYB11Generator import *
from GenericBodyForce import *
from RestartMethods import *

@PYB11template("Dimension")
class PointPotential(GenericBodyForce):

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
               G = "double",
               mass = "double",
               coreRadius = "double",
               origin = "const Vector",
               metric = ("const Tensor", "Tensor::one")):
        "PointPotential constructor"

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
    @PYB11const
    def extraEnergy(self):
        "Many physics packages will have their own representations of energy in the system (gravitational potential energy, radiative losses, etc.)"
        return "Scalar"

    #...........................................................................
    # Methods
    @PYB11const
    def specificPotential(self, r="const Vector&"):
        "The specific potential at a postion"
        return "Scalar"

    #...........................................................................
    # Properties
    G = PYB11property("Scalar", "G", "G", doc="The gravitational constant")
    mass = PYB11property("Scalar", "mass", "mass", doc="The point mass")
    coreRadius = PYB11property("Scalar", "coreRadius", "coreRadius", doc="The core softening radius")
    origin = PYB11property("const Vector&", "origin", "origin", returnpolicy="reference_internal", doc="The point mass position")
    metric = PYB11property("const Tensor&", "metric", "metric", returnpolicy="reference_internal", doc="metric to map relative positions (with respect to origin)")
    ftimestep = PYB11property("double", "ftimestep", "ftimestep", doc="The current time step scaling factor.")
    potential = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "potential", returnpolicy="reference_internal", doc="The last computed potential")

#-------------------------------------------------------------------------------
# Inject restart methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, PointPotential)
