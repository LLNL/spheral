from pybindgen import *

from ref_return_value import *
from PhysicsModule import generatePhysicsVirtualBindings

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class Damage:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir, dims):

        self.dims = dims

        # Includes.
        mod.add_include('"%s/DamageTypes.hh"' % srcdir)
    
        # Namespaces.
        space = mod.add_cpp_namespace("Spheral")

        self.TensorStrainAlgorithm = space.add_enum("TensorStrainAlgorithm", [("BenzAsphaugStrain", "Spheral::TensorStrainAlgorithm::BenzAsphaugStrain"),
                                                                              ("StrainHistory", "Spheral::TensorStrainAlgorithm::StrainHistory"),
                                                                              ("MeloshRyanAsphaugStrain", "Spheral::TensorStrainAlgorithm::MeloshRyanAsphaugStrain"),
                                                                              ("PlasticStrain", "Spheral::TensorStrainAlgorithm::PlasticStrain"),
                                                                              ("PseudoPlasticStrain", "Spheral::TensorStrainAlgorithm::PseudoPlasticStrain")])
        self.EffectiveDamageAlgorithm = space.add_enum("EffectiveDamageAlgorithm", [("CopyDamage", "Spheral::EffectiveDamageAlgorithm::CopyDamage"),
                                                                                    ("MaxDamage", "Spheral::EffectiveDamageAlgorithm::MaxDamage"),
                                                                                    ("SampledDamage", "Spheral::EffectiveDamageAlgorithm::SampledDamage")])

        self.EffectiveFlawAlgorithm = space.add_enum("EffectiveFlawAlgorithm", [("FullSpectrumFlaws", "Spheral::EffectiveFlawAlgorithm::FullSpectrumFlaws"),
                                                                                ("MinFlaw", "Spheral::EffectiveFlawAlgorithm::MinFlaw"),
                                                                                ("MaxFlaw", "Spheral::EffectiveFlawAlgorithm::MaxFlaw"),
                                                                                ("InverseSumFlaws", "Spheral::EffectiveFlawAlgorithm::InverseSumFlaws"),
                                                                                ("SampledFlaws", "Spheral::EffectiveFlawAlgorithm::SampledFlaws")])

        for dim in self.dims:
            exec('''
Physics%(dim)id = findObject(space, "Physics%(dim)id")
self.DamageModel%(dim)id = addObject(space, "DamageModel%(dim)id", parent=Physics%(dim)id, allow_subclassing=True)
self.TensorDamageModel%(dim)id = addObject(space, "TensorDamageModel%(dim)id", parent=self.DamageModel%(dim)id, allow_subclassing=True)
self.JohnsonCookDamage%(dim)id = addObject(space, "JohnsonCookDamage%(dim)id", parent=Physics%(dim)id, allow_subclassing=True)
self.addWeibullDistributionFunctions(space, %(dim)i)
self.addComputeFragmentField(space, %(dim)i)
''' % {"dim" : dim})

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        for dim in self.dims:
            exec('''
self.generateDamageModelBindings(self.DamageModel%(dim)id, %(dim)i)
self.generateTensorDamageModelBindings(self.TensorDamageModel%(dim)id, %(dim)i)
self.generateJohnsonCookDamageBindings(self.JohnsonCookDamage%(dim)id, %(dim)i)
''' % {"dim" : dim})

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return []

    #---------------------------------------------------------------------------
    # DamageModel
    #---------------------------------------------------------------------------
    def generateDamageModelBindings(self, x, ndim):

        me = "Spheral::DamageModel%id" % ndim
        dim = "Spheral::Dim<%i> " % ndim
        solidnodelist = "Spheral::SolidNodeList%id" % ndim
        vectordoublefield = "Spheral::VectorDoubleField%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        tablekernel = "Spheral::TableKernel%id" % ndim
        fileio = "Spheral::FileIO"

        # Constructors.
        x.add_constructor([refparam(solidnodelist, "nodeList"),
                           constrefparam(tablekernel, "W"),
                           param("double", "crackGrowthMultiplier"),
                           param("EffectiveFlawAlgorithm", "effectiveFlawAlgorithm"),
                           refparam(vectordoublefield, "flaws")])

        # Methods.
        x.add_method("computeScalarDDDt", None, [constrefparam(database, "dataBase"),
                                                 constrefparam(state, "state"),
                                                 param("double", "time"),
                                                 param("double", "dt"),
                                                 refparam(scalarfield, "DDDt")], is_const=True)
        x.add_method("registerState", None, [refparam(database, "dataBase"),
                                             refparam(state, "state")], is_virtual=True)
        x.add_method("preStepInitialize", None, [constrefparam(database, "dataBase"),
                                                 refparam(state, "state"),
                                                 refparam(derivatives, "derivatives")], is_virtual=True)
        x.add_method("postStateUpdate", None, [param("double", "t"),
                                               param("double", "dt"),
                                               constrefparam(database, "dataBase"),
                                               refparam(state, "state"),
                                               refparam(derivatives, "derivatives")], is_virtual=True)
        x.add_method("cullToWeakestFlaws", None, [])
        x.add_method("flawsForNode", "vector_of_double", [param("int", "index")], is_const=True)
        x.add_function_as_method("const_reference_as_pointer",
                                 retval(ptr(solidnodelist), reference_existing_object=True),
                                 [param(me, "self")],
                                 template_parameters = [me, solidnodelist, "&%s::nodeList" % me],
                                 foreign_cpp_namespace = "Spheral",
                                 custom_name = "nodeList")
        x.add_function_as_method("youngsModulusFromDamageModel",
                                 retval(ptr(scalarfield), reference_existing_object=True),
                                 [param(me, "self")],
                                 template_parameters = [dim],
                                 custom_name = "youngsModulus")
        x.add_function_as_method("longitudinalSoundSpeedFromDamageModel",
                                 retval(ptr(scalarfield), reference_existing_object=True),
                                 [param(me, "self")],
                                 template_parameters = [dim],
                                 custom_name = "longitudinalSoundSpeed")
        x.add_function_as_method("const_reference_as_pointer",
                                 retval(ptr(vectordoublefield), reference_existing_object=True),
                                 [param(me, "self")],
                                 template_parameters = [me, vectordoublefield, "&%s::flaws" % me],
                                 foreign_cpp_namespace = "Spheral",
                                 custom_name = "flaws")
        x.add_function_as_method("const_reference_as_pointer",
                                 retval(ptr(scalarfield), reference_existing_object=True),
                                 [param(me, "self")],
                                 template_parameters = [me, scalarfield, "&%s::effectiveFlaws" % me],
                                 foreign_cpp_namespace = "Spheral",
                                 custom_name = "effectiveFlaws")
        const_ref_return_value(x, me, "%s::kernel" % me, tablekernel, [], "kernel")

        x.add_method("sumActivationEnergiesPerNode", scalarfield, [], is_const=True)
        x.add_method("numFlawsPerNode", scalarfield, [], is_const=True)

        x.add_method("label", "std::string", [], is_const=True, is_virtual=True)
        x.add_method("dumpState", None, [refparam(fileio, "FileIO"), refparam("std::string", "pathName")], is_const=True, is_virtual=True)
        x.add_method("restoreState", None, [constrefparam(fileio, "FileIO"), refparam("std::string", "pathName")], is_virtual=True)

        # Attributes.
        x.add_instance_attribute("crackGrowthMultiplier", "double", getter="crackGrowthMultiplier", is_const=True)
        x.add_instance_attribute("effectiveFlawAlgorithm", "EffectiveFlawAlgorithm", getter="effectiveFlawAlgorithm", is_const=True)
        x.add_instance_attribute("criticalNodesPerSmoothingScale", "double", getter="criticalNodesPerSmoothingScale", setter="criticalNodesPerSmoothingScale")
        x.add_instance_attribute("excludeNodes", "vector_of_int", getter="excludeNodes", setter="excludeNodes")

        return

    #---------------------------------------------------------------------------
    # TensorDamageModel
    #---------------------------------------------------------------------------
    def generateTensorDamageModelBindings(self, x, ndim):

        me = "Spheral::TensorDamageModel%id" % ndim
        dim = "Spheral::Dim<%i> " % ndim
        solidnodelist = "Spheral::SolidNodeList%id" % ndim
        vectordoublefield = "Spheral::VectorDoubleField%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        tablekernel = "Spheral::TableKernel%id" % ndim
        fileio = "Spheral::FileIO"

        # Constructors.
        x.add_constructor([refparam(solidnodelist, "nodeList"),
                           param("TensorStrainAlgorithm", "strainAlgorithm"),
                           param("EffectiveDamageAlgorithm", "effectiveDamageAlgorithm"),
                           param("bool", "useDamageGradient"),
                           constrefparam(tablekernel, "kernel"),
                           param("double", "crackGrowthMultiplier"),
                           param("EffectiveFlawAlgorithm", "flawAlgorithm"),
                           param("double", "criticalDamageThreshold"),
                           param("bool", "damageInCompression"),
                           refparam(vectordoublefield, "flaws")])

        # Physics interface.
        generatePhysicsVirtualBindings(x, ndim, False)

        # Methods.
        x.add_method("applyGhostBoundaries", None, [refparam(state, "state"), refparam(derivatives, "derivatives")], is_virtual=True)
        x.add_method("enforceBoundaries", None, [refparam(state, "state"), refparam(derivatives, "derivatives")], is_virtual=True)
        x.add_function_as_method("const_reference_as_pointer",
                                 retval(ptr(symtensorfield), reference_existing_object=True),
                                 [param(me, "self")],
                                 template_parameters = [me, symtensorfield, "&%s::strain" % me],
                                 foreign_cpp_namespace = "Spheral",
                                 custom_name = "strain")
        x.add_function_as_method("const_reference_as_pointer",
                                 retval(ptr(symtensorfield), reference_existing_object=True),
                                 [param(me, "self")],
                                 template_parameters = [me, symtensorfield, "&%s::effectiveStrain" % me],
                                 foreign_cpp_namespace = "Spheral",
                                 custom_name = "effectiveStrain")
        x.add_function_as_method("const_reference_as_pointer",
                                 retval(ptr(scalarfield), reference_existing_object=True),
                                 [param(me, "self")],
                                 template_parameters = [me, scalarfield, "&%s::DdamageDt" % me],
                                 foreign_cpp_namespace = "Spheral",
                                 custom_name = "DdamageDt")
        x.add_function_as_method("const_reference_as_pointer",
                                 retval(ptr(symtensorfield), reference_existing_object=True),
                                 [param(me, "self")],
                                 template_parameters = [me, symtensorfield, "&%s::newEffectiveDamage" % me],
                                 foreign_cpp_namespace = "Spheral",
                                 custom_name = "newEffectiveDamage")
        x.add_function_as_method("const_reference_as_pointer",
                                 retval(ptr(vectorfield), reference_existing_object=True),
                                 [param(me, "self")],
                                 template_parameters = [me, vectorfield, "&%s::newDamageGradient" % me],
                                 foreign_cpp_namespace = "Spheral",
                                 custom_name = "newDamageGradient")

        x.add_method("label", "std::string", [], is_const=True, is_virtual=True)
        x.add_method("dumpState", None, [refparam(fileio, "FileIO"), refparam("std::string", "pathName")], is_const=True, is_virtual=True)
        x.add_method("restoreState", None, [constrefparam(fileio, "FileIO"), refparam("std::string", "pathName")], is_virtual=True)

        # Attributes.
        x.add_instance_attribute("strainAlgorithm", "TensorStrainAlgorithm", getter="strainAlgorithm", is_const=True)
        x.add_instance_attribute("effectiveDamageAlgorithm", "EffectiveDamageAlgorithm", getter="effectiveDamageAlgorithm", is_const=True)
        x.add_instance_attribute("useDamageGradient", "bool", getter="useDamageGradient", setter="useDamageGradient")
        x.add_instance_attribute("damageInCompression", "bool", getter="damageInCompression", setter="damageInCompression")
        x.add_instance_attribute("criticalDamageThreshold", "double", getter="criticalDamageThreshold", setter="criticalDamageThreshold")

        return

    #---------------------------------------------------------------------------
    # JohnsonCookDamage
    #---------------------------------------------------------------------------
    def generateJohnsonCookDamageBindings(self, x, ndim):

        me = "Spheral::JohnsonCookDamage%id" % ndim
        dim = "Spheral::Dim<%i> " % ndim
        solidnodelist = "Spheral::SolidNodeList%id" % ndim
        vectordoublefield = "Spheral::VectorDoubleField%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        tablekernel = "Spheral::TableKernel%id" % ndim
        fileio = "Spheral::FileIO"

        # Constructors.
        x.add_constructor([refparam(solidnodelist, "nodeList"),
                           constrefparam(scalarfield, "D1"),
                           constrefparam(scalarfield, "D2"),
                           param("double", "D3"),
                           param("double", "D4"),
                           param("double", "D5"),
                           param("double", "epsilondot0"),
                           param("double", "Tcrit"),
                           param("double", "sigmamax"),
                           param("double", "efailmin")])

        # Physics interface.
        generatePhysicsVirtualBindings(x, ndim, False)

        # Methods.
        x.add_method("applyGhostBoundaries", None, [refparam(state, "state"), refparam(derivatives, "derivatives")], is_virtual=True)
        x.add_method("enforceBoundaries", None, [refparam(state, "state"), refparam(derivatives, "derivatives")], is_virtual=True)
        x.add_function_as_method("const_reference_as_pointer",
                                 retval(ptr(solidnodelist), reference_existing_object=True),
                                 [param(me, "self")],
                                 template_parameters = [me, solidnodelist, "&%s::nodeList" % me],
                                 foreign_cpp_namespace = "Spheral",
                                 custom_name = "nodeList")
        x.add_function_as_method("const_reference_as_pointer",
                                 retval(ptr(scalarfield), reference_existing_object=True),
                                 [param(me, "self")],
                                 template_parameters = [me, scalarfield, "&%s::failureStrain" % me],
                                 foreign_cpp_namespace = "Spheral",
                                 custom_name = "failureStrain")
        x.add_function_as_method("const_reference_as_pointer",
                                 retval(ptr(scalarfield), reference_existing_object=True),
                                 [param(me, "self")],
                                 template_parameters = [me, scalarfield, "&%s::meltSpecificEnergy" % me],
                                 foreign_cpp_namespace = "Spheral",
                                 custom_name = "meltSpecificEnergy")
        x.add_function_as_method("const_reference_as_pointer",
                                 retval(ptr(symtensorfield), reference_existing_object=True),
                                 [param(me, "self")],
                                 template_parameters = [me, symtensorfield, "&%s::newEffectiveDamage" % me],
                                 foreign_cpp_namespace = "Spheral",
                                 custom_name = "newEffectiveDamage")
        x.add_function_as_method("const_reference_as_pointer",
                                 retval(ptr(scalarfield), reference_existing_object=True),
                                 [param(me, "self")],
                                 template_parameters = [me, scalarfield, "&%s::D1" % me],
                                 foreign_cpp_namespace = "Spheral",
                                 custom_name = "D1")
        x.add_function_as_method("const_reference_as_pointer",
                                 retval(ptr(scalarfield), reference_existing_object=True),
                                 [param(me, "self")],
                                 template_parameters = [me, scalarfield, "&%s::D2" % me],
                                 foreign_cpp_namespace = "Spheral",
                                 custom_name = "D2")
        x.add_method("label", "std::string", [], is_const=True, is_virtual=True)
        x.add_method("dumpState", None, [refparam(fileio, "FileIO"), refparam("std::string", "pathName")], is_const=True, is_virtual=True)
        x.add_method("restoreState", None, [constrefparam(fileio, "FileIO"), refparam("std::string", "pathName")], is_virtual=True)

        # Attributes.
        x.add_instance_attribute("D3", "double", getter="D3", is_const=True)
        x.add_instance_attribute("D4", "double", getter="D4", is_const=True)
        x.add_instance_attribute("D5", "double", getter="D5", is_const=True)
        x.add_instance_attribute("epsilondot0", "double", getter="epsilondot0", is_const=True)
        x.add_instance_attribute("Tcrit", "double", getter="Tcrit", is_const=True)
        x.add_instance_attribute("sigmamax", "double", getter="sigmamax", is_const=True)
        x.add_instance_attribute("efailmin", "double", getter="efailmin", is_const=True)

        return

    #---------------------------------------------------------------------------
    # weibull flaw distributions
    #---------------------------------------------------------------------------
    def addWeibullDistributionFunctions(self, space, ndim):

        dim = "Spheral::Dim<%i> " % ndim
        fluidnodelist = "Spheral::FluidNodeList%id" % ndim
        vectordoublefield = "Spheral::VectorDoubleField%id" % ndim

        space.add_function("weibullFlawDistributionBenzAsphaug",
                           vectordoublefield,
                           [param("double", "volume"),
                            param("double", "volumeStretchFactor"),
                            param("int", "seed"),
                            param("double", "kWeibull"),
                            param("double", "mWeibull"),
                            constrefparam(fluidnodelist, "nodeList"),
                            param("int", "minFlawsPerNode"),
                            param("int", "minTotalFlaws")],
                           template_parameters = [dim],
                           custom_name = "weibullFlawDistributionBenzAsphaug%id" % ndim)

        space.add_function("weibullFlawDistributionOwen",
                           vectordoublefield,
                           [param("int", "seed"),
                            param("double", "kWeibull"),
                            param("double", "mWeibull"),
                            constrefparam(fluidnodelist, "nodeList"),
                            param("int", "minFlawsPerNode"),
                            param("double", "volumeMultiplier", default_value="1.0")],
                           template_parameters = [dim],
                           custom_name = "weibullFlawDistributionOwen%id" % ndim)

        return

    #---------------------------------------------------------------------------
    # computeFragmentField
    #---------------------------------------------------------------------------
    def addComputeFragmentField(self, space, ndim):

        dim = "Spheral::Dim<%i> " % ndim
        intfield = "Spheral::IntField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim

        space.add_function("computeFragmentField", intfield,
                           [constrefparam(nodelist, "nodeList"),
                            param("double", "linkRadius"),
                            constrefparam(symtensorfield, "damage"),
                            param("double", "damageThreshold"),
                            param("bool", "assignDustToFragments")],
                           template_parameters = [dim],
                           custom_name = "computeFragmentField%id" % ndim)
        return
