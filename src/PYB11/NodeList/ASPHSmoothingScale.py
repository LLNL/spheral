from PYB11Generator import *
from SmoothingScaleBase import *
from SmoothingScaleAbstractMethods import *

#-------------------------------------------------------------------------------
# ASPHSmoothingScale
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
class ASPHSmoothingScale(SmoothingScaleBase):

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using Tensor = typename %(Dimension)s::Tensor;
    using SymTensor = typename %(Dimension)s::SymTensor;
    using ScalarField = Field<%(Dimension)s, Scalar>;
    using VectorField = Field<%(Dimension)s, Vector>;
    using TensorField = Field<%(Dimension)s, Tensor>;
    using SymTensorField = Field<%(Dimension)s, SymTensor>;
"""

    def pyinit(self,
               W = "const TableKernel<%(Dimension)s>&",
               targetNperh = "const double",
               numPoints = ("const size_t", "0u")):
        "Constructor: setting numPoints == 0 implies create lookup tables with same number of points as TableKernel W"

    @PYB11const
    def equivalentNodesPerSmoothingScale(self,
                                         lambdaPsi = "Scalar"):
        "Compute the nPerh that corresponds to the given eigenvalue of second moment tensor (1/sqrt of the eigenvalue actually)"
        return "Scalar"

    @PYB11const
    def equivalentLambdaPsi(self,
                            nPerh = "Scalar"):
        "Compute the lambda_psi eigenvalue that corresponds to the  nPerh value"
        return "Scalar"

    #...........................................................................
    # Properties
    targetNperh = PYB11property("double", doc="The target nPerh for building the ASPH nperh lookup tables")
    minNperh = PYB11property("double", doc="The lower limit for looking up the effective nPerh")
    maxNperh = PYB11property("double", doc="The upper limit for looking up the effective nPerh")
    nPerhInterpolator = PYB11property(doc = "nperh(x) interpolator")
    WsumInterpolator = PYB11property(doc = "Wsum(x) interpolator")

#-------------------------------------------------------------------------------
# Add the abstract interface
#-------------------------------------------------------------------------------
PYB11inject(SmoothingScaleAbstractMethods, ASPHSmoothingScale, virtual=True)
