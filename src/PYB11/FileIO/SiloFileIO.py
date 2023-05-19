#-------------------------------------------------------------------------------
# SiloFileIO
#-------------------------------------------------------------------------------
from PYB11Generator import *
from FileIO import *
from FileIOAbstractMethods import *
from FileIOTemplateMethods import *
from spheralDimensions import *
dims = spheralDimensions()

class SiloFileIO(FileIO):
    "Handle FileIO for silo files"

    #...........................................................................
    # Constructors
    def pyinit0(self):
        "Default constructor"

    def pyinit1(self,
                filename = "const std::string",
                access = "AccessType"):
        "Open a silo file with a given file name and access"

    #...........................................................................
    # Override abstract methods
    @PYB11virtual
    def open(self,
             fileName = "const std::string",
             access = "AccessType"):
        "Open a file for IO"
        return "void"

    @PYB11virtual
    def close(self):
        "Close the current file we're pointing at"
        return "void"

#-------------------------------------------------------------------------------
# Override the required virtual interface
#-------------------------------------------------------------------------------
PYB11inject(FileIOAbstractMethods, SiloFileIO, virtual=True, pure_virtual=False)
#PYB11inject(FileIOTemplateMethods, SiloFileIO)
