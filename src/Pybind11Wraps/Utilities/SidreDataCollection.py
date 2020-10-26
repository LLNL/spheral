#-------------------------------------------------------------------------------
# SidreDataCollection
#-------------------------------------------------------------------------------
from PYB11Generator import *
#from axom import sidre
#from Field import Field

class SidreDataCollection:

    #...........................................................................
    # Constructors
    def pyinit(self):
        "Default constructor"

    #...........................................................................
    # Methods
    @PYB11template("Dimension", "DataType")
    @PYB11returnpolicy("reference")
    def alloc_view(self, view_name="const std::string&", field="const Spheral::Field<%(Dimension)s, %(DataType)s>&"):
        "Create a sidre view containing the data in a Field"
        return "axom::sidre::View *"

    def printDataStore(self):
        "Print out the data store from root"
        return

    alloc_viewIntField1D = PYB11TemplateMethod(alloc_view, ("Dim<1>", "int"), pyname="alloc_view")
    
    #...........................................................................
    # Attributes
