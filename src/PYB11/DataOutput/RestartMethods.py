#-------------------------------------------------------------------------------
# Helper class to inject standard restart methods
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11ignore
class RestartMethods:

    @PYB11virtual
    @PYB11const
    def label(self):
        "Label for restart files"
        return "std::string"

    @PYB11virtual
    @PYB11const
    def dumpState(self, file="FileIO&", pathName="const std::string&"):
        "Serialize under the given path in a FileIO object"
        return "void"

    @PYB11virtual
    def restoreState(self, file="const FileIO&", pathName="const std::string&"):
        "Restore state from the given path in a FileIO object"
        return "void"
