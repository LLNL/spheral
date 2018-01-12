from pybindgen import *

from PBGutils import *

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class Gperftools:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir, dims):

        # Includes.
        mod.add_include('"gperftools/profiler.h"')
    
        return

    #---------------------------------------------------------------------------
    # Generate bindings.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        mod.add_function("ProfilerStart", "int", [param("char*", "fname", transfer_ownership=False)], 
                         docstring="Start the Gperftools profile to the given file.")
        mod.add_function("ProfilerStop", None, [],
                         docstring="Stop the Gperftools profile.")
        mod.add_function("ProfilerFlush", None, [],
                         docstring="Flush output to the Gperftools profile file.")
        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return []

