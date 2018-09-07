from pybindgen import *

from PBGutils import *
from MaterialModule import generateEquationOfStateVirtualBindings

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class Helmholtz:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir, dims):

        self.dims = dims

        # Includes.
        mod.add_include('"%s/HelmholtzTypes.hh"' % srcdir)
    
        # Namespace.
        space = mod.add_cpp_namespace("Spheral")

        for dim in self.dims:
            exec('''
self.EquationOfState%(dim)id = findObject(space, "EquationOfState%(dim)id")
self.HelmholtzEquationOfState%(dim)id = addObject(space, "HelmholtzEquationOfState%(dim)id", allow_subclassing=True, parent=self.EquationOfState%(dim)id)
''' % {"dim" : dim})

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):
        
        for dim in self.dims:
            exec('''
self.generateHelmholtzEquationOfStateBindings(self.HelmholtzEquationOfState%(dim)id, %(dim)i)
                ''' % {"dim" : dim})                                
        return

#---------------------------------------------------------------------------
# HelmholtzEquationOfState
#---------------------------------------------------------------------------
    def generateHelmholtzEquationOfStateBindings(self, x, ndim):

        scalarfield = "Spheral::ScalarField%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim

        # Constructor.
        x.add_constructor([constrefparam("PhysicalConstants", "constants"),
                           param("double", "minimumPressure", default_value="-std::numeric_limits<double>::max()"),
                           param("double", "maximumPressure", default_value="std::numeric_limits<double>::max()"),
                           param("double", "minimumTemperature", default_value="-std::numeric_limits<double>::min()"),
                           param("MaterialPressureMinType", "minPressureType", default_value="Spheral::MaterialPressureMinType::PressureFloor"),
                           param("double", "abar0", default_value="13.6"),
                           param("double", "zbar0", default_value="6.8")])
                       
        # Attributes.
        x.add_instance_attribute("needUpdate", "bool", getter="getUpdateStatus", setter="setUpdateStatus", is_const=False)

        # Methods

        generateEquationOfStateVirtualBindings(x, ndim, False)
            
        return
