"""
Spheral DataOuput module.

Provides the fundamental classes for Restart in Spheral.
"""

from PYB11Generator import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
includes = ['"Geometry/Dimension.hh"',
            '"DataOutput/RestartRegistrar.hh"',
            '"RestartableObject.hh"',
            '"FileIO/FileIO.hh"',
            "<vector>",
            "<string>"]
            
#-------------------------------------------------------------------------------
# Namespaces the module is in.
#-------------------------------------------------------------------------------
namespaces = ["Spheral"]

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

    # The instance attribute.  We expose this as a property of the class.
    @PYB11static
    @PYB11cppname("instance")
    @PYB11ignore
    def getinstance(self):
        return "RestartRegistrar&"
    instance = property(getinstance, doc="The static RestartRegistrar instance.")

#-------------------------------------------------------------------------------
# RestartableObject
#-------------------------------------------------------------------------------
class RestartableObject:
    "The base class for building restartable python objects in Spheral."

    def pyinit(self,
               pyself = "py::object&",
               priority = "const unsigned"):
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
                  file = "FileIO&",
                  pathName = "const std::string&"):
        "Write this objects state to the file under the given path."
        return "void"

    @PYB11virtual
    def restoreState(self,
                     file = "const FileIO&",
                     pathName = ("const std::string&", '"pathName"')):
        "Read the state for this object from the given file and path."
        return "void"
