#-------------------------------------------------------------------------------
# TensorCRKSPHViscosity
#-------------------------------------------------------------------------------
from PYB11Generator import *
from ArtificialViscosity import *
from ArtificialViscosityAbstractMethods import *

@PYB11template("Dimension")
class TensorCRKSPHViscosity(ArtificialViscosity):

    typedefs = """
    typedef %(Dimension)s DIM;
    typedef typename DIM::Scalar Scalar;
    typedef typename DIM::Vector Vector;
    typedef typename DIM::Tensor Tensor;
    typedef typename DIM::SymTensor SymTensor;
    typedef typename DIM::ThirdRankTensor ThirdRankTensor;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               Clinear = ("const Scalar", "1.0"),
               Cquadratic = ("const Scalar", "1.0")):
        "TensorCRKSPHViscosity constructor"

    #...........................................................................
    # Methods
    @PYB11virtual
    @PYB11const
    def label(self):
        return "std::string"

#-------------------------------------------------------------------------------
# Inject abstract interface
#-------------------------------------------------------------------------------
PYB11inject(ArtificialViscosityAbstractMethods, TensorCRKSPHViscosity, virtual=True, pure_virtual=False)
