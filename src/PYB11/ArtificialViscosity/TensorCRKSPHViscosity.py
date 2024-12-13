#-------------------------------------------------------------------------------
# TensorCRKSPHViscosity
#-------------------------------------------------------------------------------
from PYB11Generator import *
from TensorMonaghanViscosity import *

@PYB11template("Dimension")
@PYB11template_dict({"QPiType", "typename %(Dimension)s::Tensor"})
class TensorCRKSPHViscosity(ArtificialViscosity):

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
               Clinear = "const Scalar",
               Cquadratic = "const Scalar",
               kernel = "const TableKernel<%(Dimension)s>&",
               order = "const RKOrder"):
        "TensorCRKSPHViscosity constructor"

    #...........................................................................
    # Methods
    @PYB11virtual
    @PYB11const
    def requireVelocityGradient(self):
        "We need the velocity gradient and set this to true"
        return "bool"

    @PYB11virtual
    def updateVelocityGradient(self,
                               dataBase = "const DataBase<%(Dimension)s>&",
                               state = "const State<%(Dimension)s>&",
                               derivs = "const StateDerivatives<%(Dimension)s>&"):
        "Update the locally stored velocity gradient"
        return "void"

    @PYB11virtual
    @PYB11const
    def requireReproducingKernels(self):
        "Some physics algorithms require reproducing kernels."
        return "std::set<RKOrder>"

    @PYB11virtual
    @PYB11const
    def label(self):
        return "std::string"

    #...........................................................................
    # Properties
    order = PYB11Property("RKOrder", "order", "order",
                          doc="RK order to use to estimate velocity gradient")
