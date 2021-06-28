#-------------------------------------------------------------------------------
# ContactModel base class
#-------------------------------------------------------------------------------
from PYB11Generator import *
from ContactModelBaseAbstractMethods import *
#from RestartMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralDEM")
class ContactModelBase:

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::SymTensor SymTensor;
"""

    #...........................................................................
    # Constructors
    def pyinit(self):
        "ContactModel constructor"


#-------------------------------------------------------------------------------
# Inject abstract interface
#-------------------------------------------------------------------------------
PYB11inject(ContactModelBaseAbstractMethods, ContactModelBase, pure_virtual=True)
#PYB11inject(RestartMethods, ContactModel)
