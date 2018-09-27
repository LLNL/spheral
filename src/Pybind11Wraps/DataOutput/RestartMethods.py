#-------------------------------------------------------------------------------
# Helper class to inject standard restart methods
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11ignore
class RestartMethods:

    @PYB11virtual
    @PYB11const
    def label(self):
        "The label used to create a path in the restart file"
        return "std::string"

    @PYB11virtual
    @PYB11const
    def dumpState(self,
                  file = "FileIO&",
                  pathName = "const std::string&"):
        "Write the restart state to the given path in the given file"
        return "void"

    @PYB11virtual
    def restoreState(self,
                     file = "const FileIO&",
                     pathName = "const std::string&"):
        "Restore our state from the given path & file"
        return "void"
