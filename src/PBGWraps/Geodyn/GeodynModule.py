from pybindgen import *

from ref_return_value import *
from MaterialModule import  generateEquationOfStateVirtualBindings

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class Geodyn:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir):

        # Includes.
        mod.add_include('"%s/SolidMaterial/Geodyn.hh"' % topsrcdir)

        # Namespace.
        Spheral = mod.add_cpp_namespace("Spheral")
        Material = Spheral.add_cpp_namespace("Material")
        PhysicsSpace = Spheral.add_cpp_namespace("PhysicsSpace")
        self.space = Spheral.add_cpp_namespace("SolidMaterial")

        self.dimSet = (1, 2, 3)

        for dim in self.dimSet:
            exec('''
EquationOfState%(dim)id = Material.add_class("EquationOfState",
                                             template_parameters=["Spheral::Dim<%(dim)i>"],
                                             custom_name="EquationOfState%(dim)id", 
                                             import_from_module="SpheralModules.Spheral.Material",
                                             allow_subclassing=True)
StrengthModel%(dim)id = findObject(self.space, "StrengthModel%(dim)id")
Physics%(dim)id = PhysicsSpace.add_class("Physics",
                                         template_parameters=["Spheral::Dim<%(dim)i>"],
                                         custom_name="Physics%(dim)id", 
                                         import_from_module="SpheralModules.Spheral.PhysicsSpace",
                                         allow_subclassing=True)
self.Geodyn%(dim)id = self.space.add_class("Geodyn",
                                           template_parameters=["Spheral::Dim<%(dim)i>"],
                                           custom_name="Geodyn%(dim)id", 
                                           parent=[EquationOfState%(dim)id, StrengthModel%(dim)id, Physics%(dim)id],
                                           allow_subclassing=True)
''' % {"dim" : dim})

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        for dim in self.dimSet:
            exec('''
generateGeodynBindings(self.Geodyn%(dim)id, %(dim)i)
''' % {"dim" : dim})

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return []

#---------------------------------------------------------------------------
# Geodyn
#---------------------------------------------------------------------------
def generateGeodynBindings(x, ndim):

    dim = "Spheral::Dim< %i >" % ndim
    me = "Spheral::SolidMaterial::Geodyn%id" % ndim

    # Constructors.
    x.add_constructor([constrefparam("PhysicalConstants", "constants"),
                       param("double", "minimumPressure", default_value="-std::numeric_limits<double>::max()"),
                       param("double", "maximumPressure", default_value="std::numeric_limits<double>::max()"),
                       param("MaterialPressureMinType", "minPressureType", default_value="PressureFloor")])

    # Generic EOS interface.
    generateEquationOfStateVirtualBindings(x, ndim, False)

    # Methods.

    # Attributes.

    return
