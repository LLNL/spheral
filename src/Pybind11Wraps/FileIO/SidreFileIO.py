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
        "Construct a sidre datastore with a name and access"

    def pyinit2(self,
                filename = "const std::string",
                access = "AccessType",
                numFiles = "int"):
        "Construct a sidre datastore with a name and access and set number of restart files"

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
