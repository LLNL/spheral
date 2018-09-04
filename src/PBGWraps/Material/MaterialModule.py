from pybindgen import *

from PBGutils import *

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class Material:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir, dims):

        self.dims = dims

        # Includes.
        mod.add_include('"%s/MaterialTypes.hh"' % srcdir)
    
        # Namespace.
        space = mod.add_cpp_namespace("Spheral")

        # Expose types.
        self.MaterialPressureMinType = space.add_enum("MaterialPressureMinType", [("PressureFloor", "Spheral::MaterialPressureMinType::PressureFloor"),
                                                                                  ("ZeroPressure", "Spheral::MaterialPressureMinType::ZeroPressure")])
        self.PhysicalConstants = addObject(space, "PhysicalConstants", allow_subclassing=True)
        for dim in self.dims:
            exec('''
self.EquationOfState%(dim)id = addObject(space, "EquationOfState%(dim)id", allow_subclassing=True)
self.GammaLawGas%(dim)id = addObject(space, "GammaLawGas%(dim)id", allow_subclassing=True, parent=self.EquationOfState%(dim)id)
self.PolytropicEquationOfState%(dim)id = addObject(space, "PolytropicEquationOfState%(dim)id", allow_subclassing=True, parent=self.EquationOfState%(dim)id)
self.IsothermalEquationOfState%(dim)id = addObject(space, "IsothermalEquationOfState%(dim)id", allow_subclassing=True, parent=self.EquationOfState%(dim)id)
''' % {"dim" : dim})

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):
        self.generatePhysicalConstantsBindings(self.PhysicalConstants)
        for dim in self.dims:
            exec('''
self.generateEquationOfStateBindings(self.EquationOfState%(dim)id, %(dim)i)
self.generateGammaLawGasBindings(self.GammaLawGas%(dim)id, %(dim)i)
self.generatePolytropicEquationOfStateBindings(self.PolytropicEquationOfState%(dim)id, %(dim)i)
self.generateIsothermalEquationOfStateBindings(self.IsothermalEquationOfState%(dim)id, %(dim)i)
''' % {"dim"   : dim})

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return ["Material"]

    #---------------------------------------------------------------------------
    # PhysicalConstants
    #---------------------------------------------------------------------------
    def generatePhysicalConstantsBindings(self, x):

        me = "PhysicalConstants"

        # Constructor.
        x.add_constructor([param("double", "unitLm"),
                           param("double", "unitMkg"),
                           param("double", "unitTsec")])
        x.add_constructor([constrefparam(me, "rhs")])

        # Attributes
        x.add_instance_attribute("unitLengthMeters", "double", getter="unitLengthMeters", is_const=True)
        x.add_instance_attribute("unitMassKg", "double", getter="unitMassKg", is_const=True)
        x.add_instance_attribute("unitTimeSec", "double", getter="unitTimeSec", is_const=True)
        x.add_instance_attribute("unitMassDensity", "double", getter="unitMassDensity", is_const=True)
        x.add_instance_attribute("protonMass", "double", getter="protonMass", is_const=True)
        x.add_instance_attribute("electronMass", "double", getter="electronMass", is_const=True)
        x.add_instance_attribute("electronCharge", "double", getter="electronCharge", is_const=True)
        x.add_instance_attribute("G", "double", getter="G", is_const=True)
        x.add_instance_attribute("c", "double", getter="c", is_const=True)
        x.add_instance_attribute("kB", "double", getter="kB", is_const=True)
        x.add_instance_attribute("NAvogadro", "double", getter="Navogadro", is_const=True)
        x.add_instance_attribute("molarGasConstant", "double", getter="molarGasConstant", is_const=True)
        x.add_instance_attribute("kelvinsToEnergyPerMole", "double", getter="kelvinsToEnergyPerMole", is_const=True)
        x.add_instance_attribute("stefanBoltzmannConstant", "double", getter="stefanBoltzmannConstant", is_const=True)

        return

    #---------------------------------------------------------------------------
    # EquationOfState
    #---------------------------------------------------------------------------
    def generateEquationOfStateBindings(self, x, ndim):

        # Constructors.
        x.add_constructor([constrefparam("PhysicalConstants", "constants"),
                           param("double", "minimumPressure", default_value="-std::numeric_limits<double>::max()"),
                           param("double", "maximumPressure", default_value="std::numeric_limits<double>::max()"),
                           param("MaterialPressureMinType", "minPressureType", default_value="MaterialPressureMinType::PressureFloor")])

        # Methods.
        x.add_method("specificThermalEnergyForPressure", "double",
                     [param("double", "Ptarget"),
                      param("double", "rho"),
                      param("double", "epsMin"),
                      param("double", "epsMax"),
                      param("double", "epsTol"),
                      param("double", "Ptol"),
                      param("double", "maxIterations", default_value="100")],
                     is_virtual = True,
                     is_const = True,
                     docstring = "Look up a specific thermal energy that gives the requested pressure at the requested density.")
        x.add_method("applyPressureLimits", "double", [param("double", "P")], is_const=True)

        # Attributes.
        x.add_instance_attribute("constants", "PhysicalConstants", getter="constants", is_const=True)
        x.add_instance_attribute("minimumPressure", "double", getter="minimumPressure", setter="minimumPressure")
        x.add_instance_attribute("maximumPressure", "double", getter="maximumPressure", setter="maximumPressure")
        x.add_instance_attribute("minimumPressureType", "MaterialPressureMinType", getter="minimumPressureType", setter="minimumPressureType")

        generateEquationOfStateVirtualBindings(x, ndim, True)

    #---------------------------------------------------------------------------
    # GammaLawGas
    #---------------------------------------------------------------------------
    def generateGammaLawGasBindings(self, x, ndim):

        # Constructor.
        x.add_constructor([param("double", "gamma"),
                           param("double", "mu"),
                           constrefparam("PhysicalConstants", "constants"),
                           param("double", "minimumPressure", default_value="-std::numeric_limits<double>::max()"),
                           param("double", "maximumPressure", default_value="std::numeric_limits<double>::max()"),
                           param("MaterialPressureMinType", "minPressureType", default_value="MaterialPressureMinType::PressureFloor")])

        # Attributes.
        x.add_instance_attribute("gamma", "double", getter="getGamma", setter="setGamma")
        x.add_instance_attribute("mu", "double", getter="getMolecularWeight", setter="setMolecularWeight")

        # Methods.
        x.add_method("pressure", "double", [param("double", "massDensity"),
                                            param("double", "specificThermalEnergy")],
                     is_const=True)
        x.add_method("temperature", "double", [param("double", "massDensity"),
                                               param("double", "specificThermalEnergy")],
                     is_const=True)
        x.add_method("specificThermalEnergy", "double", [param("double", "massDensity"),
                                                         param("double", "temperature")],
                     is_const=True)
        x.add_method("specificHeat", "double", [param("double", "massDensity"),
                                                param("double", "tempernature")],
                     is_const=True)
        x.add_method("soundSpeed", "double", [param("double", "massDensity"),
                                              param("double", "specificThermalEnergy")],
                     is_const=True)
        x.add_method("gamma", "double", [param("double", "massDensity"),
                                         param("double", "specificThermalEnergy")],
                     is_const=True)
        x.add_method("bulkModulus", "double", [param("double", "massDensity"),
                                               param("double", "specificThermalEnergy")],
                     is_const=True)
        x.add_method("entropy", "double", [param("double", "massDensity"),
                                           param("double", "specificThermalEnergy")],
                     is_const=True)

        # Add the EOS virual interface.
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
                           constrefparam("PhysicalConstants", "constants"),
                           param("double", "minimumPressure", default_value="-std::numeric_limits<double>::max()"),
                           param("double", "maximumPressure", default_value="std::numeric_limits<double>::max()"),
                           param("MaterialPressureMinType", "minPressureType", default_value="MaterialPressureMinType::PressureFloor")])

        # Attributes.
        x.add_instance_attribute("polytropicConstant", "double", getter="polytropicConstant", is_const=True)
        x.add_instance_attribute("polytropicIndex", "double", getter="polytropicIndex", is_const=True)
        x.add_instance_attribute("gamma_", "double", getter="gamma", is_const=True)
        x.add_instance_attribute("molecularWeight", "double", getter="molecularWeight", is_const=True)
        x.add_instance_attribute("externalPressure", "double", getter="externalPressure", setter="setExternalPressure")

        # Methods.
        x.add_method("pressure", "double", [param("double", "massDensity"),
                                            param("double", "specificThermalEnergy")],
                     is_const=True)
        x.add_method("temperature", "double", [param("double", "massDensity"),
                                               param("double", "specificThermalEnergy")],
                     is_const=True)
        x.add_method("specificThermalEnergy", "double", [param("double", "massDensity"),
                                                         param("double", "temperature")],
                     is_const=True)
        x.add_method("specificHeat", "double", [param("double", "massDensity"),
                                                param("double", "tempernature")],
                     is_const=True)
        x.add_method("soundSpeed", "double", [param("double", "massDensity"),
                                              param("double", "specificThermalEnergy")],
                     is_const=True)
        x.add_method("gamma", "double", [param("double", "massDensity"),
                                         param("double", "specificThermalEnergy")],
                     is_const=True)
        x.add_method("bulkModulus", "double", [param("double", "massDensity"),
                                               param("double", "specificThermalEnergy")],
                     is_const=True)
        x.add_method("entropy", "double", [param("double", "massDensity"),
                                           param("double", "specificThermalEnergy")],
                     is_const=True)

        # Add the EOS virual interface.
        generateEquationOfStateVirtualBindings(x, ndim, False)

        return

    #---------------------------------------------------------------------------
    # IsothermalEquationOfState
    #---------------------------------------------------------------------------
    def generateIsothermalEquationOfStateBindings(self, x, ndim):

        # Constructor.
        x.add_constructor([param("double", "K"),
                           param("double", "mu"),
                           constrefparam("PhysicalConstants", "constants"),
                           param("double", "minimumPressure", default_value="-std::numeric_limits<double>::max()"),
                           param("double", "maximumPressure", default_value="std::numeric_limits<double>::max()"),
                           param("MaterialPressureMinType", "minPressureType", default_value="MaterialPressureMinType::PressureFloor")])

        # Attributes.
        x.add_instance_attribute("K", "double", getter="K", is_const=True)
        x.add_instance_attribute("molecularWeight", "double", getter="molecularWeight", is_const=True)
        x.add_instance_attribute("externalPressure", "double", getter="externalPressure", setter="setExternalPressure")

        # Methods.
        x.add_method("pressure", "double", [param("double", "massDensity"),
                                            param("double", "specificThermalEnergy")],
                     is_const=True)
        x.add_method("temperature", "double", [param("double", "massDensity"),
                                               param("double", "specificThermalEnergy")],
                     is_const=True)
        x.add_method("specificThermalEnergy", "double", [param("double", "massDensity"),
                                                         param("double", "temperature")],
                     is_const=True)
        x.add_method("specificHeat", "double", [param("double", "massDensity"),
                                                param("double", "tempernature")],
                     is_const=True)
        x.add_method("soundSpeed", "double", [param("double", "massDensity"),
                                              param("double", "specificThermalEnergy")],
                     is_const=True)
        x.add_method("gamma", "double", [param("double", "massDensity"),
                                         param("double", "specificThermalEnergy")],
                     is_const=True)
        x.add_method("bulkModulus", "double", [param("double", "massDensity"),
                                               param("double", "specificThermalEnergy")],
                     is_const=True)
        x.add_method("entropy", "double", [param("double", "massDensity"),
                                           param("double", "specificThermalEnergy")],
                     is_const=True)

        # Add the EOS virual interface.
        generateEquationOfStateVirtualBindings(x, ndim, False)

        return

#-------------------------------------------------------------------------------
# EquationOfState virtual bindings
#-------------------------------------------------------------------------------
def generateEquationOfStateVirtualBindings(x, ndim, pureVirtual):

    scalarfield = "Spheral::ScalarField%id" % ndim

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
    x.add_method("setEntropy", None, [refparam(scalarfield, "entropy"),
                                      constrefparam(scalarfield, "massDensity"),
                                      constrefparam(scalarfield, "specificThermalEnergy")],
                 is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)

    x.add_method("valid", "bool", [], is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)

    return

