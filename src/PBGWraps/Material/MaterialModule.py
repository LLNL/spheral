from pybindgen import *

import sys
sys.path.append("..")
from PBGutils import *

#-------------------------------------------------------------------------------
# EquationOfState virtual bindings
#-------------------------------------------------------------------------------
def generateEquationOfStateVirtualBindings(x, ndim, pureVirtual):

    scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim

    # Constructor.

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

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class Material:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod):

        # Includes.
        mod.add_include('"Material/MaterialTypes.hh"')
    
        # Namespace.
        Spheral = mod.add_cpp_namespace("Spheral")
        space = Spheral.add_cpp_namespace("Material")

        self.unitSet = ("CGS", "MKS", "Cosmological", "Solar")
        self.dimSet = (1, 2, 3)

        # Expose types.
        for dim in self.dimSet:
            exec('''self.EquationOfState%(dim)id = addObject(space, "EquationOfState%(dim)id", allow_subclassing=True)''' % {"dim" : dim})
        for units in self.unitSet:
            exec('''self.%(units)sUnits = addObject(space, "%(units)s")''' % {"units" : units})
            for dim in self.dimSet:
                exec('''
self.GammaLawGas%(units)s%(dim)id = addObject(space, "GammaLawGas%(units)s%(dim)id", allow_subclassing=True, parent=self.EquationOfState%(dim)id)                
self.PolytropicEquationOfState%(units)s%(dim)id = addObject(space, "PolytropicEquationOfState%(units)s%(dim)id", allow_subclassing=True, parent=self.EquationOfState%(dim)id)                
self.IsothermalEquationOfState%(units)s%(dim)id = addObject(space, "IsothermalEquationOfState%(units)s%(dim)id", allow_subclassing=True, parent=self.EquationOfState%(dim)id)                
''' % {"units" : units,
       "dim" : dim})

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        for dim in self.dimSet:
            exec('self.generateEquationOfStateBindings(self.EquationOfState%(dim)id, %(dim)i)' % {"dim" : dim})

        for units in self.unitSet:
            exec('self.generateUnitsBindings(self.%sUnits)' % units)

            for dim in self.dimSet:
                exec('''
self.generateGammaLawGasBindings(self.GammaLawGas%(units)s%(dim)id, %(dim)i)
self.generatePolytropicEquationOfStateBindings(self.PolytropicEquationOfState%(units)s%(dim)id, %(dim)i)
self.generateIsothermalEquationOfStateBindings(self.IsothermalEquationOfState%(units)s%(dim)id, %(dim)i)
''' % {"units" : units,
       "dim"   : dim})

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return ["Material"]

    #---------------------------------------------------------------------------
    # CGSUnits
    #---------------------------------------------------------------------------
    def generateUnitsBindings(self, x):

        # Constructor.
        x.add_constructor([])

        # Attributes
        x.add_instance_attribute("protonMass", "double", getter="protonMass", is_const=True)
        x.add_instance_attribute("electronMass", "double", getter="electronMass", is_const=True)
        x.add_instance_attribute("electronCharge", "double", getter="electronCharge", is_const=True)
        x.add_instance_attribute("G", "double", getter="G", is_const=True)
        x.add_instance_attribute("c", "double", getter="c", is_const=True)
        x.add_instance_attribute("kB", "double", getter="kB", is_const=True)
        x.add_instance_attribute("unitLengthMeters", "double", getter="unitLengthMeters", is_const=True)
        x.add_instance_attribute("unitMassKg", "double", getter="unitMassKg", is_const=True)
        x.add_instance_attribute("unitTimeSec", "double", getter="unitTimeSec", is_const=True)
        x.add_instance_attribute("NAvogadro", "double", is_const=True)
        x.add_instance_attribute("MolarGasConstant", "double", is_const=True)
        x.add_instance_attribute("KelvinsToEnergyPerMole", "double", is_const=True)

        return

    #---------------------------------------------------------------------------
    # EquationOfState
    #---------------------------------------------------------------------------
    def generateEquationOfStateBindings(self, x, ndim):

        # Constructors.
        x.add_constructor([param("double", "minimumPressure", default_value="-std::numeric_limits<double>::max()"),
                           param("double", "maximumPressure", default_value="std::numeric_limits<double>::max()")])

        # Attributes.
        x.add_instance_attribute("minimumPressure", "double", getter="minimumPressure", setter="minimumPressure")
        x.add_instance_attribute("maximumPressure", "double", getter="maximumPressure", setter="maximumPressure")

        generateEquationOfStateVirtualBindings(x, ndim, True)

    #---------------------------------------------------------------------------
    # GammaLawGas
    #---------------------------------------------------------------------------
    def generateGammaLawGasBindings(self, x, ndim):

        # Constructor.
        x.add_constructor([param("double", "gamma"),
                           param("double", "mu"),
                           param("double", "minimumPressure", default_value="-std::numeric_limits<double>::max()"),
                           param("double", "maximumPressure", default_value="std::numeric_limits<double>::max()")])

        # Attributes.
        x.add_instance_attribute("gamma", "double", getter="getGamma", setter="setGamma")
        x.add_instance_attribute("mu", "double", getter="getMolecularWeight", setter="setMolecularWeight")

        generateEquationOfStateVirtualBindings(x, ndim, False)

        return

    #---------------------------------------------------------------------------
    # PolytropicEquationOfState
    #---------------------------------------------------------------------------
    def generatePolytropicEquationOfStateBindings(self, x, ndim):

        # Constructor.
        x.add_constructor([param("double", "K"),
                           param("double", "index"),
                           param("double", "mu"),
                           param("double", "minimumPressure", default_value="-std::numeric_limits<double>::max()"),
                           param("double", "maximumPressure", default_value="std::numeric_limits<double>::max()")])

        # Attributes.
        x.add_instance_attribute("polytropicConstant", "double", getter="polytropicConstant", is_const=True)
        x.add_instance_attribute("polytropicIndex", "double", getter="polytropicIndex", is_const=True)
        x.add_instance_attribute("gamma_", "double", getter="gamma", is_const=True)
        x.add_instance_attribute("molecularWeight", "double", getter="molecularWeight", is_const=True)
        x.add_instance_attribute("externalPressure", "double", getter="externalPressure", setter="setExternalPressure")

        generateEquationOfStateVirtualBindings(x, ndim, False)

        return

    #---------------------------------------------------------------------------
    # IsothermalEquationOfState
    #---------------------------------------------------------------------------
    def generateIsothermalEquationOfStateBindings(self, x, ndim):

        # Constructor.
        x.add_constructor([param("double", "K"),
                           param("double", "mu"),
                           param("double", "minimumPressure", default_value="-std::numeric_limits<double>::max()"),
                           param("double", "maximumPressure", default_value="std::numeric_limits<double>::max()")])

        # Attributes.
        x.add_instance_attribute("K", "double", getter="K", is_const=True)
        x.add_instance_attribute("molecularWeight", "double", getter="molecularWeight", is_const=True)
        x.add_instance_attribute("externalPressure", "double", getter="externalPressure", setter="setExternalPressure")

        generateEquationOfStateVirtualBindings(x, ndim, False)

        return
