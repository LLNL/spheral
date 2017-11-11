from pybindgen import *

from PBGutils import *

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class OpenMP:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir, dims):

        # Includes.
        mod.add_include('"OpenMP/OpenMPhelpers.hh"')
    
        return

    #---------------------------------------------------------------------------
    # Generate bindings.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        mod.add_function("wrap_omp_get_thread_num", "int", [], 
                         custom_name = "omp_get_thread_num",
                         docstring = "Get the OpenMP thread ID.")
        mod.add_function("wrap_omp_get_num_threads", "int", [], 
                         custom_name = "omp_get_num_threads",
                         docstring = "Get the number of OpenMP threads.")

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return []

