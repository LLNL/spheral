from pybindgen import *
from PBGutils import *

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class DataOutput:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod):

        # Includes.
        mod.add_include('"DataOutput/DataOutputTypes.hh"')
    
        # Namespace.
        Spheral = mod.add_cpp_namespace("Spheral")
        space = Spheral.add_cpp_namespace("DataOutput")

        # Expose types.
        self.RestartRegistrar = addObject(space, "RestartRegistrar", is_singleton=True)
        self.RestartableObject = addObject(space, "RestartableObject", allow_subclassing=True)
        
        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):
        self.addRestartRegistrarMethods()
        self.addRestartableObjectMethods()

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return ["DataOutput"]

    #---------------------------------------------------------------------------
    # Add RestartRegistrar methods.
    #---------------------------------------------------------------------------
    def addRestartRegistrarMethods(self):

        x = self.RestartRegistrar

        # Get the instance.
        x.add_method("instancePtr", retval("Spheral::DataOutput::RestartRegistrar*", caller_owns_return=True), [], is_static=True, custom_name="instance")

        # Methods.
        x.add_method("removeExpiredPointers", None, [])
        x.add_method("uniqueLabels", "vector_of_string", [], is_const=True)
        x.add_method("printLabels", None, [], is_const=True)
        x.add_method("dumpState", None, [refparam("FileIO", "file")], is_const=True)
        x.add_method("restoreState", None, [constrefparam("FileIO", "file")], is_const=True)

        return

    #---------------------------------------------------------------------------
    # Add RestartableObject methods.
    #---------------------------------------------------------------------------
    def addRestartableObjectMethods(self):

        fileio = "Spheral::FileIOSpace::FileIO"

        x = self.RestartableObject

        # Constructors.
        x.add_constructor([])
        x.add_constructor([param("int", "priority")])

        # Methods.
        x.add_method("label", "std::string", [], is_const=True, is_pure_virtual=True)
        x.add_method("dumpState", None, [refparam(fileio, "file"), param("const std::string", "pathName")], is_const=True, is_pure_virtual=True)
        x.add_method("restoreState", None, [constrefparam(fileio, "file"), param("const std::string", "pathName")], is_pure_virtual=True)

        return
