#-------------------------------------------------------------------------------
# SPHSmoothingScale
#-------------------------------------------------------------------------------
from PYB11Generator import *
from SmoothingScaleBase import *

@PYB11template("Dimension")
class SPHSmoothingScale(SmoothingScaleBase):

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
               HUpdate = "HEvolutionType",
               W = "const TableKernel<%(Dimension)s>&"):
        "SPHSmoothingScale constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def initializeProblemStartup(self,
                                 dataBase = "DataBase<%(Dimension)s>&"):
        """An optional hook to initialize once when the problem is starting up.
Typically this is used to size arrays once all the materials and NodeLists have
been created.  It is assumed after this method has been called it is safe to
call Physics::registerState for instance to create full populated State objects."""
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
    def label(self):
        return "std::string"

    @PYB11virtual
    @PYB11const
    def dumpState(self, file="FileIO&", pathName="const std::string&"):
        "Serialize under the given path in a FileIO object"
        return "void"

    @PYB11virtual
    def restoreState(self, file="const FileIO&", pathName="const std::string&"):
        "Restore state from the given path in a FileIO object"
        return "void"

    #...........................................................................
    # Attributes
    WT = PYB11property("const TableKernel<%(Dimension)s>&", "WT", doc="The interpolation kernel")
    zerothMoment = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "zerothMoment", doc="The zeroth moment storage FieldList")
    firstMoment = PYB11property("const FieldList<%(Dimension)s, Vector>&", "firstMoment", doc="The first moment storage FieldList")
