from pybindgen import *

import sys
sys.path.append("..")
sys.path.append("../Material")
from PBGutils import *
from MaterialModule import generateEquationOfStateVirtualBindings

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class Helmholtz:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod):

        # Includes.
        mod.add_include('"Helmholtz/HelmholtzTypes.hh"')
    
        # Namespace.
        Spheral = mod.add_cpp_namespace("Spheral")
        space = Spheral.add_cpp_namespace("Material")

        self.dimSet = (1, 2, 3)

        for dim in self.dimSet:
            exec('''
self.EquationOfState%(dim)id = findObject(space, "EquationOfState%(dim)id")
self.HelmholtzEquationOfState%(dim)id = addObject(space, "HelmholtzEquationOfState%(dim)id", allow_subclassing=True, parent=self.EquationOfState%(dim)id)
''' % {"dim" : dim})

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):
        
        for dim in self.dimSet:
            exec('''
self.generateHelmholtzEquationOfStateBindings(self.HelmholtzEquationOfState%(dim)id, %(dim)i)
                ''' % {"dim" : dim})                                
        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return []


#---------------------------------------------------------------------------
# HelmholtzEquationOfState
#---------------------------------------------------------------------------
    def generateHelmholtzEquationOfStateBindings(self, x, ndim):

        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        nodelist = "Spheral::NodeSpace::NodeList%id" % ndim

        # Constructor.
        x.add_constructor([constrefparam(nodelist, "myNodeList"),
                           constrefparam("PhysicalConstants", "constants"),
                           param("double", "minimumPressure", default_value="-std::numeric_limits<double>::max()"),
                           param("double", "maximumPressure", default_value="std::numeric_limits<double>::max()"),
                           param("double", "minimumTemperature", default_value="-std::numeric_limits<double>::min()"),
                           param("double", "maximumTemperature", default_value="std::numeric_limits<double>::max()"),
                           param("MaterialPressureMinType", "minPressureType", default_value="PressureFloor"),
                           param("double", "abar0", default_value="13.6"),
                           param("double", "zbar0", default_value="6.8")])
                       
        # Attributes.
        x.add_instance_attribute("needUpdate", "bool", getter="getUpdateStatus", setter="setUpdateStatus", is_const=False)

        # Methods

        generateEquationOfStateVirtualBindings(x, ndim, False)
            
        return
