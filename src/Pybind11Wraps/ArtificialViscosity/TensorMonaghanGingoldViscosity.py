#-------------------------------------------------------------------------------
# TensorMonaghanGingoldViscosity
#-------------------------------------------------------------------------------
from PYB11Generator import *
from ArtificialViscosity import *
from ArtificialViscosityAbstractMethods import *

@PYB11template("Dimension")
class TensorMonaghanGingoldViscosity(ArtificialViscosity):
    """A modified form of the Monaghan & Gingold viscosity, extended to tensor formalism.
This method is described in
Owen, J Michael (2004), 'A tensor artficial visocity for SPH', Journal of Computational Physics, 201(2), 601-629
"""

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
               Clinear = "const Scalar",
               Cquadratic = "const Scalar"):
        "TensorMonaghanGingoldViscosity constructor"

    #...........................................................................
    # Methods
    @PYB11virtual
    @PYB11const
    def label(self):
        return "std::string"
    
#-------------------------------------------------------------------------------
# Inject abstract interface
#-------------------------------------------------------------------------------
PYB11inject(ArtificialViscosityAbstractMethods, TensorMonaghanGingoldViscosity, virtual=True, pure_virtual=False)
