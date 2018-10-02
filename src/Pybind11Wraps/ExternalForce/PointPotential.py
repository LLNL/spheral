#-------------------------------------------------------------------------------
# PointPotential base class
#-------------------------------------------------------------------------------
from PYB11Generator import *
from GenericBodyForce import *

@PYB11template("Dimension")
class PointPotential(GenericBodyForce):

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
               G = "double",
               mass = "double",
               coreRadius = "double",
               origin = "const Vector&"):
        "PointPotential constructor"

    #...........................................................................
    # Virtual methods
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
    @PYB11const
    def label(self):
        "It's useful to have labels for Physics packages.  We'll require this to have the same signature as the restart label."
        return "std::string"

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
    G = PYB11property("Scalar", "G", "setG", doc="The gravitational constant")
    mass = PYB11property("Scalar", "mass", "setMass", doc="The point mass")
    coreRadius = PYB11property("Scalar", "coreRadius", "setCoreRadius", doc="The core softening radius")
    origin = PYB11property("const Vector&", "origin", "setOrigin", returnpolicy="reference_internal", doc="The point mass position")
    deltaPotentialFraction = PYB11property("Scalar", "deltaPotentialFraction", "setDeltaPotentialFraction", doc="The max allowed fraction potential change, for setting the timestep")
