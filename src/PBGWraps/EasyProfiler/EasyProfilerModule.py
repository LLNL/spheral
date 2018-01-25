from pybindgen import *

from PBGutils import *

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class EasyProfiler:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir, dims):

        # Includes.
        mod.add_include('"EasyProfiler/EasyProfiler_helpers.hh"')
    
        return

    #---------------------------------------------------------------------------
    # Generate bindings.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        mod.add_function("EasyProfilerStart", None, [], 
                         docstring = "Fire up the Easy Profiler timers.")
        mod.add_function("EasyProfilerDump", None, [param("std::string", "basename")], 
                         docstring = "Write the Easy Profiler output (per MPI process) to files starting with the given basename.")

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return []

