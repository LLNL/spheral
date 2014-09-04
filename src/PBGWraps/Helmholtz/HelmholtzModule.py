from pybindgen import *

import sys
sys.path.append("..")
from PBGutils import *

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class Material:

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

        # Expose types.
        self.MaterialPressureMinType = space.add_enum("MaterialPressureMinType", ["PressureFloor", "ZeroPressure"])
        self.PhysicalConstants = addObject(space, "PhysicalConstants", allow_subclassing=True)
        for dim in self.dimSet:
            exec('''
self.EquationOfState%(dim)id = findObject(Material, "EquationOfState%(dim)id")
self.HelmholtzEquationOfState%(dim)id = addObject(space, "HelmholtzEquationOfState%(dim)id", allow_subclassing=True, parent=self.EquationOfState%(dim)id)
''' % {"dim" : dim})

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):
        
        for dim in self.dimSet:
            exec('''
                generateHelmholtzEquationOfStateBindings(self.HelmholtzEquationOfState%(dim)id, %(dim)i)
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
    #const_ref_return_value(x, me, "%s::abar" % me, scalarfield, [], "abar")
    #const_ref_return_value(x, me, "%s::zbar" % me, scalarfield, [], "zbar")

    generateEquationOfStateVirtualBindings(x, ndim, False)
            
    return

#-------------------------------------------------------------------------------
# EquationOfState virtual bindings
#-------------------------------------------------------------------------------
def generateEquationOfStateVirtualBindings(x, ndim, pureVirtual):

    scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim

    # Methods.
    x.add_method("setPressure", None, [refparam(scalarfield, "pressure"),
                                       constrefparam(scalarfield, "massDensity"),
                                       constrefparam(scalarfield, "specificThermalEnergy")],
                 is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
    x.add_method("setTemperature", None, [refparam(scalarfield, "temperature"),
                                       constrefparam(scalarfield, "massDensity"),
                                       constrefparam(scalarfield, "specificThermalEnergy")],
                 is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
    x.add_method("setSpecificThermalEnergy", None, [refparam(scalarfield, "specificThermalEnergy"),
                                       constrefparam(scalarfield, "massDensity"),
                                       constrefparam(scalarfield, "temperature")],
                 is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
    x.add_method("setSpecificHeat", None, [refparam(scalarfield, "specificHeat"),
                                       constrefparam(scalarfield, "massDensity"),
                                       constrefparam(scalarfield, "temperature")],
                 is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
    x.add_method("setSoundSpeed", None, [refparam(scalarfield, "soundSpeed"),
                                       constrefparam(scalarfield, "massDensity"),
                                       constrefparam(scalarfield, "specificThermalEnergy")],
                 is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
    x.add_method("setGammaField", None, [refparam(scalarfield, "gammaField"),
                                       constrefparam(scalarfield, "massDensity"),
                                       constrefparam(scalarfield, "specificThermalEnergy")],
                 is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
    x.add_method("setBulkModulus", None, [refparam(scalarfield, "bulkModulus"),
                                       constrefparam(scalarfield, "massDensity"),
                                       constrefparam(scalarfield, "specificThermalEnergy")],
                 is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)

    x.add_method("pressure", "double", [param("double", "massDensity"),
                                        param("double", "specificThermalEnergy")],
                 is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
    x.add_method("temperature", "double", [param("double", "massDensity"),
                                           param("double", "specificThermalEnergy")],
                 is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
    x.add_method("specificThermalEnergy", "double", [param("double", "massDensity"),
                                                     param("double", "temperature")],
                 is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
    x.add_method("specificHeat", "double", [param("double", "massDensity"),
                                            param("double", "temperature")],
                 is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
    x.add_method("soundSpeed", "double", [param("double", "massDensity"),
                                          param("double", "specificThermalEnergy")],
                 is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
    x.add_method("gamma", "double", [param("double", "massDensity"),
                                     param("double", "specificThermalEnergy")],
                 is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
    x.add_method("bulkModulus", "double", [param("double", "massDensity"),
                                           param("double", "specificThermalEnergy")],
                 is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)

    x.add_method("valid", "bool", [], is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)

    return

