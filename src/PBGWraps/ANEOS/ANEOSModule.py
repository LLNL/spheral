from pybindgen import *

from ref_return_value import *
from MaterialModule import  generateEquationOfStateVirtualBindings

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class ANEOS:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod):

        # Includes.
        mod.add_include('"ANEOS/ANEOSTypes.hh"')

        # Namespace.
        Spheral = mod.add_cpp_namespace("Spheral")
        Material = Spheral.add_cpp_namespace("Material")
        self.space = Spheral.add_cpp_namespace("SolidMaterial")

        self.dimSet = (1, 2, 3)

        for dim in self.dimSet:
            exec('''
EquationOfState%(dim)id = findObject(Material, "EquationOfState%(dim)id")
self.ANEOS%(dim)id = addObject(self.space, "ANEOS%(dim)id", parent=EquationOfState%(dim)id, allow_subclassing=True)
''' % {"dim" : dim})

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        for dim in self.dimSet:
            exec('''
generateANEOSBindings(self.ANEOS%(dim)id, %(dim)i)
''' % {"dim" : dim})

        self.space.add_function("initializeANEOS", None,
                                [param("std::string", "in_filename"),
                                 param("std::string", "out_filename"),
                                 param("vector_of_int", "izetl")],
                                docstring = "Intialize the ANEOS equation of state package.")

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return []

#---------------------------------------------------------------------------
# ANEOS
#---------------------------------------------------------------------------
def generateANEOSBindings(x, ndim):

    dim = "Spheral::Dim< %i >" % ndim
    me = "Spheral::SolidMaterial::ANEOS%id" % ndim

    # Constructors.
    x.add_constructor([param("int", "materialNumber"),
                       param("unsigned int", "numRhoVals"),
                       param("unsigned int", "numTvals"),
                       param("double", "rhoMin"),
                       param("double", "rhoMax"),
                       param("double", "Tmin"),
                       param("double", "Tmax"),
                       constrefparam("PhysicalConstants", "constants"),
                       param("double", "externalPressure", default_value="0.0"),
                       param("double", "minimumPressure", default_value="-std::numeric_limits<double>::max()"),
                       param("double", "maximumPressure", default_value="std::numeric_limits<double>::max()")])

    # Generic EOS interface.
    generateEquationOfStateVirtualBindings(x, ndim, False)

    # Methods.
    x.add_function_as_method("ANEOS_STEvals", "vector_of_vector_of_double", 
                             [param(me, "self")],
                             template_parameters = [dim],
                             custom_name = "specificThermalEnergyVals")

    # Attributes.
    x.add_instance_attribute("materialNumber", "int", getter="materialNumber", is_const=True)
    x.add_instance_attribute("numRhoVals", "unsigned int", getter="numRhoVals", is_const=True)
    x.add_instance_attribute("numTvals", "unsigned int", getter="numTvals", is_const=True)
    x.add_instance_attribute("rhoMin", "double", getter="rhoMin", is_const=True)
    x.add_instance_attribute("rhoMax", "double", getter="rhoMax", is_const=True)
    x.add_instance_attribute("Tmin", "double", getter="Tmin", is_const=True)
    x.add_instance_attribute("Tmax", "double", getter="Tmax", is_const=True)
    x.add_instance_attribute("externalPressure", "double", getter="externalPressure", setter="externalPressure")

    return
