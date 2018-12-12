"""
Spheral DataOuput module.

Provides the fundamental classes for Restart in Spheral.
"""

from PYB11Generator import *
from SpheralCommon import *
from RestartableObject import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"DataOutput/RestartRegistrar.hh"',
                  '"RestartableObject.hh"',
                  '"FileIO/FileIO.hh"']
            
#-------------------------------------------------------------------------------
# Namespaces the module is in.
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

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
    @PYB11cppname("instancePtr")
    @PYB11returnpolicy("take_ownership")
    def instance(self):
        return "RestartRegistrar*"
