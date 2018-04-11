"""
Spheral DataOuput module.

Provides the fundamental classes for Restart in Spheral.
"""

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
includes = ["<vector>",
            "<string>",
            '"Geometry/Dimension"',
            '"DataOutput/RestartRegistrar.hh"',
            '"RestartableObject.hh"',
            '"FileIO/FileIO.hh"']

#-------------------------------------------------------------------------------
# Namespaces the module is in.
#-------------------------------------------------------------------------------
namespaces = ["Spheral", "DataOutput"]

#-------------------------------------------------------------------------------
# RestartRegistrar
#-------------------------------------------------------------------------------
class RestartRegistrar:
    singlegton = True

    def removeExpiredPointers(self):
        "Clear all expired objects from the RestartRegistrar."
        return

    def uniqueLabels(self):
        "Return a std::vector of the labels for each known restart handle."
        return

    def printLabels(self):
        "Print out the labels for all known restart handles."
        return

    def dumpState(self):
        "Dump the state of all restartable handles"
        return

    def restoreState(self):
        "Restore the state of all restartable handles"
        return

#-------------------------------------------------------------------------------
# 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# 
#-------------------------------------------------------------------------------
