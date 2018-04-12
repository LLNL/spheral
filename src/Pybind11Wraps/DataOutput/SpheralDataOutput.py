"""
Spheral DataOuput module.

Provides the fundamental classes for Restart in Spheral.
"""

import sys
sys.path.append("..")
from PYB11Generator import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
includes = ["<vector>",
            "<string>",
            '"Geometry/Dimension"',
            '"DataOutput/RestartRegistrar.hh"',
            '"DataOutput/RestartableObject.hh"',
            '"FileIO/FileIO.hh"']

#-------------------------------------------------------------------------------
# Namespaces the module is in.
#-------------------------------------------------------------------------------
namespaces = ["Spheral", "DataOutput"]

#-------------------------------------------------------------------------------
# RestartRegistrar
#-------------------------------------------------------------------------------
@PYB11singleton
class RestartRegistrar:
    "A singleton object that holds all objects currently registered for restart."

    def removeExpiredPointers(self):
        "Clear all expired objects from the RestartRegistrar."
        return

    @PYB11const
    def uniqueLabels(self):
        "Return a std::vector of the labels for each known restart handle."
        return

    @PYB11const
    def printLabels(self):
        "Print out the labels for all known restart handles."
        return

    @PYB11const
    def dumpState(self):
        "Dump the state of all restartable handles"
        return

    @PYB11const
    def restoreState(self):
        "Restore the state of all restartable handles"
        return

#-------------------------------------------------------------------------------
# RestartableObject
#-------------------------------------------------------------------------------
class RestartableObject:

    def __init__(self,
                 args = [("PyObject*", "self"),
                         ("const unsigned", "priority")]):
        "Construct with the given object and priority."
        return

    @PYB11virtual
    @PYB11const
    def label(self):
        "Define the label for storing this object in a restart file."
        return "std::string"

    @PYB11virtual
    @PYB11const
    def dumpState(self,
                  args = [("FileIOSpace::FileIO&", "file"),
                          ("const std::string", "pathName")]):
        "Write this objects state to the file under the given path."
        return "void"

    @PYB11virtual
    def restoreState(self,
                     args = [("const FileIOSpace::FileIO&", "file"),
                             ("const std::string", "pathName")]):
        "Read the state for this object from the given file and path."
        return "void"

#-------------------------------------------------------------------------------
# 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# 
#-------------------------------------------------------------------------------
