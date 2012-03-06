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

        # Expose types.
        self.CGSUnits = addObject(space, "CGS")
        self.MKSUnits = addObject(space, "MKS")
        self.CosmologicalUnits = addObject(space, "Cosmological")

        self.EquationOfState1d = addObject(space, "EquationOfState1d", allow_subclassing=True)
        self.EquationOfState2d = addObject(space, "EquationOfState2d", allow_subclassing=True)
        self.EquationOfState3d = addObject(space, "EquationOfState3d", allow_subclassing=True)

        self.GammaLawGasCGS1d = addObject(space, "GammaLawGasCGS1d", allow_subclassing=True, parent=self.EquationOfState1d)
        self.GammaLawGasCGS2d = addObject(space, "GammaLawGasCGS2d", allow_subclassing=True, parent=self.EquationOfState2d)
        self.GammaLawGasCGS3d = addObject(space, "GammaLawGasCGS3d", allow_subclassing=True, parent=self.EquationOfState3d)

        self.GammaLawGasMKS1d = addObject(space, "GammaLawGasMKS1d", allow_subclassing=True, parent=self.EquationOfState1d)
        self.GammaLawGasMKS2d = addObject(space, "GammaLawGasMKS2d", allow_subclassing=True, parent=self.EquationOfState2d)
        self.GammaLawGasMKS3d = addObject(space, "GammaLawGasMKS3d", allow_subclassing=True, parent=self.EquationOfState3d)

        self.GammaLawGasCosmological1d = addObject(space, "GammaLawGasCosmological1d", allow_subclassing=True, parent=self.EquationOfState1d)
        self.GammaLawGasCosmological2d = addObject(space, "GammaLawGasCosmological2d", allow_subclassing=True, parent=self.EquationOfState2d)
        self.GammaLawGasCosmological3d = addObject(space, "GammaLawGasCosmological3d", allow_subclassing=True, parent=self.EquationOfState3d)

        self.PolytropicEquationOfStateCGS1d = addObject(space, "PolytropicEquationOfStateCGS1d", allow_subclassing=True, parent=self.EquationOfState1d)
        self.PolytropicEquationOfStateCGS2d = addObject(space, "PolytropicEquationOfStateCGS2d", allow_subclassing=True, parent=self.EquationOfState2d)
        self.PolytropicEquationOfStateCGS3d = addObject(space, "PolytropicEquationOfStateCGS3d", allow_subclassing=True, parent=self.EquationOfState3d)

        self.PolytropicEquationOfStateMKS1d = addObject(space, "PolytropicEquationOfStateMKS1d", allow_subclassing=True, parent=self.EquationOfState1d)
        self.PolytropicEquationOfStateMKS2d = addObject(space, "PolytropicEquationOfStateMKS2d", allow_subclassing=True, parent=self.EquationOfState2d)
        self.PolytropicEquationOfStateMKS3d = addObject(space, "PolytropicEquationOfStateMKS3d", allow_subclassing=True, parent=self.EquationOfState3d)

        self.PolytropicEquationOfStateCosmological1d = addObject(space, "PolytropicEquationOfStateCosmological1d", allow_subclassing=True, parent=self.EquationOfState1d)
        self.PolytropicEquationOfStateCosmological2d = addObject(space, "PolytropicEquationOfStateCosmological2d", allow_subclassing=True, parent=self.EquationOfState2d)
        self.PolytropicEquationOfStateCosmological3d = addObject(space, "PolytropicEquationOfStateCosmological3d", allow_subclassing=True, parent=self.EquationOfState3d)

        self.IsothermalEquationOfStateCGS1d = addObject(space, "IsothermalEquationOfStateCGS1d", allow_subclassing=True, parent=self.EquationOfState1d)
        self.IsothermalEquationOfStateCGS2d = addObject(space, "IsothermalEquationOfStateCGS2d", allow_subclassing=True, parent=self.EquationOfState2d)
        self.IsothermalEquationOfStateCGS3d = addObject(space, "IsothermalEquationOfStateCGS3d", allow_subclassing=True, parent=self.EquationOfState3d)

        self.IsothermalEquationOfStateMKS1d = addObject(space, "IsothermalEquationOfStateMKS1d", allow_subclassing=True, parent=self.EquationOfState1d)
        self.IsothermalEquationOfStateMKS2d = addObject(space, "IsothermalEquationOfStateMKS2d", allow_subclassing=True, parent=self.EquationOfState2d)
        self.IsothermalEquationOfStateMKS3d = addObject(space, "IsothermalEquationOfStateMKS3d", allow_subclassing=True, parent=self.EquationOfState3d)

        self.IsothermalEquationOfStateCosmological1d = addObject(space, "IsothermalEquationOfStateCosmological1d", allow_subclassing=True, parent=self.EquationOfState1d)
        self.IsothermalEquationOfStateCosmological2d = addObject(space, "IsothermalEquationOfStateCosmological2d", allow_subclassing=True, parent=self.EquationOfState2d)
        self.IsothermalEquationOfStateCosmological3d = addObject(space, "IsothermalEquationOfStateCosmological3d", allow_subclassing=True, parent=self.EquationOfState3d)

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        self.generateUnitsBindings(self.CGSUnits)
        self.generateUnitsBindings(self.MKSUnits)
        self.generateUnitsBindings(self.CosmologicalUnits)

        self.generateEquationOfStateBindings(self.EquationOfState1d, 1)
        self.generateEquationOfStateBindings(self.EquationOfState2d, 2)
        self.generateEquationOfStateBindings(self.EquationOfState3d, 3)

        self.generateGammaLawGasBindings(self.GammaLawGasCGS1d, 1)
        self.generateGammaLawGasBindings(self.GammaLawGasCGS2d, 2)
        self.generateGammaLawGasBindings(self.GammaLawGasCGS3d, 3)
        self.generateGammaLawGasBindings(self.GammaLawGasMKS1d, 1)
        self.generateGammaLawGasBindings(self.GammaLawGasMKS2d, 2)
        self.generateGammaLawGasBindings(self.GammaLawGasMKS3d, 3)
        self.generateGammaLawGasBindings(self.GammaLawGasCosmological1d, 1)
        self.generateGammaLawGasBindings(self.GammaLawGasCosmological2d, 2)
        self.generateGammaLawGasBindings(self.GammaLawGasCosmological3d, 3)

        self.generatePolytropicEquationOfStateBindings(self.PolytropicEquationOfStateCGS1d, 1)
        self.generatePolytropicEquationOfStateBindings(self.PolytropicEquationOfStateCGS2d, 2)
        self.generatePolytropicEquationOfStateBindings(self.PolytropicEquationOfStateCGS3d, 3)
        self.generatePolytropicEquationOfStateBindings(self.PolytropicEquationOfStateMKS1d, 1)
        self.generatePolytropicEquationOfStateBindings(self.PolytropicEquationOfStateMKS2d, 2)
        self.generatePolytropicEquationOfStateBindings(self.PolytropicEquationOfStateMKS3d, 3)
        self.generatePolytropicEquationOfStateBindings(self.PolytropicEquationOfStateCosmological1d, 1)
        self.generatePolytropicEquationOfStateBindings(self.PolytropicEquationOfStateCosmological2d, 2)
        self.generatePolytropicEquationOfStateBindings(self.PolytropicEquationOfStateCosmological3d, 3)

        self.generateIsothermalEquationOfStateBindings(self.IsothermalEquationOfStateCGS1d, 1)
        self.generateIsothermalEquationOfStateBindings(self.IsothermalEquationOfStateCGS2d, 2)
        self.generateIsothermalEquationOfStateBindings(self.IsothermalEquationOfStateCGS3d, 3)
        self.generateIsothermalEquationOfStateBindings(self.IsothermalEquationOfStateMKS1d, 1)
        self.generateIsothermalEquationOfStateBindings(self.IsothermalEquationOfStateMKS2d, 2)
        self.generateIsothermalEquationOfStateBindings(self.IsothermalEquationOfStateMKS3d, 3)
        self.generateIsothermalEquationOfStateBindings(self.IsothermalEquationOfStateCosmological1d, 1)
        self.generateIsothermalEquationOfStateBindings(self.IsothermalEquationOfStateCosmological2d, 2)
        self.generateIsothermalEquationOfStateBindings(self.IsothermalEquationOfStateCosmological3d, 3)

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
