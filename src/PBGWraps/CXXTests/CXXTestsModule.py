from pybindgen import *

from ref_return_value import *
from MaterialModule import  generateEquationOfStateVirtualBindings
from PhysicsModule import generatePhysicsVirtualBindings
from SolidMaterialModule import generateStrengthModelVirtualBindings

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class CXXTests:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir, dims):

        self.dims = dims

        # Includes.
        mod.add_include('"%s/CXXTests/testNodeIterators.hh"' % topsrcdir)
        mod.add_include('"%s/CXXTests/test_r3d_utils.hh"' % topsrcdir)

        # Namespace.
        self.space = mod.add_cpp_namespace("Spheral")

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        for ndim in self.dims:

            self.space.add_function("test_polygon_to_r2d_poly", "std::string", [])
            self.space.add_function("test_r2d_poly_to_polygon", "std::string", [])
            self.space.add_function("test_polyhedron_to_r3d_poly", "std::string", [])
            self.space.add_function("test_r3d_poly_to_polyhedron", "std::string", [])
            self.space.add_function("test_polygon_to_r2d_poly", "std::string", [])
            self.space.add_function("test_clip_polygon", "std::string", [])
            self.space.add_function("test_orphan_polygon", "std::string", [])
            self.space.add_function("test_clip_polyhedron", "std::string", [])

            exec('''
self.space.add_function("testGlobalAllNodeIterators", "std::string",
                        [constrefparam("DataBase%(ndim)id", "dataBase")],
                        template_parameters = ["%(Dim)s"],
                        custom_name = "testGlobalAllNodeIterators")
self.space.add_function("testGlobalInternalNodeIterators", "std::string",
                        [constrefparam("DataBase%(ndim)id", "dataBase")],
                        template_parameters = ["%(Dim)s"],
                        custom_name = "testGlobalInternalNodeIterators")
self.space.add_function("testGlobalGhostNodeIterators", "std::string",
                        [constrefparam("DataBase%(ndim)id", "dataBase")],
                        template_parameters = ["%(Dim)s"],
                        custom_name = "testGlobalGhostNodeIterators")
self.space.add_function("testGlobalMasterNodeIterators", "std::string",
                        [constrefparam("DataBase%(ndim)id", "dataBase")],
                        template_parameters = ["%(Dim)s"],
                        custom_name = "testGlobalMasterNodeIterators")
self.space.add_function("testGlobalCoarseNodeIterators", "std::string",
                        [constrefparam("DataBase%(ndim)id", "dataBase")],
                        template_parameters = ["%(Dim)s"],
                        custom_name = "testGlobalCoarseNodeIterators")
self.space.add_function("testGlobalRefineNodeIterators", "std::string",
                        [constrefparam("DataBase%(ndim)id", "dataBase")],
                        template_parameters = ["%(Dim)s"],
                        custom_name = "testGlobalRefineNodeIterators")

''' % {"ndim" : ndim,
       "Dim"  : "Dim<%i>" % ndim})



        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return []
