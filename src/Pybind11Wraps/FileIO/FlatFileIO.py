#-------------------------------------------------------------------------------
# FlatFileIO
#-------------------------------------------------------------------------------
from PYB11Generator import *
from FileIO import *
from FileIOAbstractMethods import *
from FileIOTemplateMethods import *
from spheralDimensions import *
dims = spheralDimensions()

class FlatFileIO(FileIO):
    "FileIO implementation for raw read/writes to a file (ascii or binary)"

    #...........................................................................
    # Constructors
    def pyinit0(self):
        "Default constructor"

    def pyinit1(self,
                filename = "const std::string",
                access = "AccessType",
                format = ("FlatFileFormat", "FlatFileFormat::ascii")):
        "Construct with a given file name, access, and file format (ascii/binary)"

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

    #...........................................................................
    # Methods
    @PYB11const
    def findPathName(self,
                     pathName = "const std::string"):
        "Move the pointer in the file to the specified path"
        return "void"

    @PYB11const
    def beginningOfFile(self):
        "Move the pointer in the file to the beginning"
        return "void"

    #...........................................................................
    # Properties
    precision = PYB11property("int", "precision", "setPrecision", doc="Precision for reading/writing numbers")
    readyToWrite = PYB11property("bool", "readyToWrite", doc="Is the FlatFileIO object valid for writing?")

#-------------------------------------------------------------------------------
# Override the required virtual interface
#-------------------------------------------------------------------------------
PYB11inject(FileIOAbstractMethods, FlatFileIO, virtual=True, pure_virtual=False)
#PYB11inject(FileIOTemplateMethods, FlatFileIO)
