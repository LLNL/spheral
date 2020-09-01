#-------------------------------------------------------------------------------
# FiniteVolumeViscosity
#-------------------------------------------------------------------------------
from PYB11Generator import *
from ArtificialViscosity import *
from ArtificialViscosityAbstractMethods import *

@PYB11template("Dimension")
class FiniteVolumeViscosity(ArtificialViscosity):

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
               Cquadratic = ("const Scalar", "1.0"),
               scalar = ("const bool", "false")):
        "FiniteVolumeViscosity constructor"

    #...........................................................................
    # Methods
    @PYB11virtual
    @PYB11const
    def label(self):
        return "std::string"

    #...........................................................................
    # Properties
    scalar = PYB11property("bool", "scalar")
    DvDx = PYB11property("const FieldList<%(Dimension)s, Tensor>&", "DvDx", returnpolicy="reference_internal")
    
#-------------------------------------------------------------------------------
# Inject abstract interface
#-------------------------------------------------------------------------------
PYB11inject(ArtificialViscosityAbstractMethods, FiniteVolumeViscosity, virtual=True, pure_virtual=False)
