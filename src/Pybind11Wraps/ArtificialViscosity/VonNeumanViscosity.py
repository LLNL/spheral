#-------------------------------------------------------------------------------
# VonNeumanViscosity
#-------------------------------------------------------------------------------
from PYB11Generator import *
from ArtificialViscosity import *
from ArtificialViscosityAbstractMethods import *
from RestartMethods import *

@PYB11template("Dimension")
class VonNeumanViscosity(ArtificialViscosity):

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
        "VonNeumanViscosity constructor"

    #...........................................................................
    # Methods
    @PYB11virtual
    @PYB11const
    def label(self):
        return "std::string"

    #...........................................................................
    # Properties
    viscousEnergy = PYB11property("const FieldList<DIM, Scalar>&", "viscousEnergy", returnpolicy="reference_internal")
    
#-------------------------------------------------------------------------------
# Inject abstract interface
#-------------------------------------------------------------------------------
PYB11inject(ArtificialViscosityAbstractMethods, VonNeumanViscosity, virtual=True, pure_virtual=False)
PYB11inject(RestartMethods, VonNeumanViscosity, virtual=True, pure_virtual=False)
