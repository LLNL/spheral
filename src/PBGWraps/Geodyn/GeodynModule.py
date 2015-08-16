from pybindgen import *

from ref_return_value import *
from MaterialModule import  generateEquationOfStateVirtualBindings
from PhysicsModule import generatePhysicsVirtualBindings
from SolidMaterialModule import generateStrengthModelVirtualBindings

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class Geodyn:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir, dims):

        self.dims = dims

        # Includes.
        mod.add_include('"%s/SolidMaterial/Geodyn.hh"' % topsrcdir)

        # Namespace.
        Spheral = mod.add_cpp_namespace("Spheral")
        Material = Spheral.add_cpp_namespace("Material")
        PhysicsSpace = Spheral.add_cpp_namespace("PhysicsSpace")
        self.space = Spheral.add_cpp_namespace("SolidMaterial")

        for dim in self.dims:
            exec('''
StrengthModel%(dim)id = findObject(self.space, "StrengthModel%(dim)id")
EquationOfState%(dim)id = findObject(Material, "EquationOfState%(dim)id")
Physics%(dim)id = findObject(PhysicsSpace, "Physics%(dim)id")

self.Geodyn%(dim)id = self.space.add_class("Geodyn",
                                           template_parameters=["Spheral::Dim<%(dim)i>"],
                                           custom_name="Geodyn%(dim)id", 
                                           #parent=[EquationOfState%(dim)id, StrengthModel%(dim)id],
                                           parent=[Physics%(dim)id, EquationOfState%(dim)id, StrengthModel%(dim)id],
                                           allow_subclassing=True)
''' % {"dim" : dim})

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        for dim in self.dims:
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

    # Generic parent interfaces.
    generateEquationOfStateVirtualBindings(x, ndim, False)
    generatePhysicsVirtualBindings(x, ndim, False)
    generateStrengthModelVirtualBindings(x, ndim, False)

    # Methods.

    # Attributes.

    return
