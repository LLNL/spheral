#-------------------------------------------------------------------------------
# SiloFileIO
#-------------------------------------------------------------------------------
from PYB11Generator import *
from injectFileIOVirtualMethods import *
from spheralDimensions import *
dims = spheralDimensions()

class SiloFileIO:
    "Handle FileIO for silo files"

    #...........................................................................
    # Constructors
    def pyinit0(self):
        "Default constructor"

    def pyinit1(self,
                filename = "const std::string",
                access = "AccessType"):
        "Open a silo file with a given file name and access"

#-------------------------------------------------------------------------------
# Override the required virtual interface
#-------------------------------------------------------------------------------
injectFileIOVirtualMethods(SiloFileIO)
