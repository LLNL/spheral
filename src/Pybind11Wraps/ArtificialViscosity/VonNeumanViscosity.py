#-------------------------------------------------------------------------------
# VonNeumanViscosity
#-------------------------------------------------------------------------------
from PYB11Generator import *
from ArtificialViscosity import *
from ArtificialViscosityAbstractMethods import *
from RestartMethods import *

@PYB11template("Dimension")
class VonNeumanViscosity(ArtificialViscosity):

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
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
    viscousEnergy = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "viscousEnergy", returnpolicy="reference_internal")
    
#-------------------------------------------------------------------------------
# Inject abstract interface
#-------------------------------------------------------------------------------
PYB11inject(ArtificialViscosityAbstractMethods, VonNeumanViscosity, virtual=True, pure_virtual=False)
PYB11inject(RestartMethods, VonNeumanViscosity, virtual=True, pure_virtual=False)
