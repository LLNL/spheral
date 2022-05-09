#-------------------------------------------------------------------------------
# SidreFileIO
#-------------------------------------------------------------------------------
from PYB11Generator import *
from FileIO import *
from FileIOAbstractMethods import *
from FileIOTemplateMethods import *
from spheralDimensions import *
dims = spheralDimensions()

class SidreFileIO(FileIO):
    "Handle FileIO for sidre file system"

    #...........................................................................
    # Constructors
    def pyinit0(self):
        "Default constructor"

    def pyinit1(self,
                filename = "const std::string",
                access = "AccessType"):
        "Save a file name and access to be used when openning the sidre datastore"

    #...........................................................................
    # Override abstract methods
    @PYB11virtual
    def open(self,
             fileName = "const std::string",
             access = "AccessType"):
        "Create a pointer to a sidre datastore for IO"
        return "void"

    @PYB11virtual
    def close(self):
        "Write the data in the sidre datastore to a file and delete the pointer"
        return "void"

#-------------------------------------------------------------------------------
# Override the required virtual interface
#-------------------------------------------------------------------------------
PYB11inject(FileIOAbstractMethods, SidreFileIO, virtual=True, pure_virtual=False)
#PYB11inject(FileIOTemplateMethods, SidreFileIO)
