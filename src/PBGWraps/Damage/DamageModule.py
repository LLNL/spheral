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
    def __init__(self, mod, srcdir, topsrcdir):

        # Includes.
        mod.add_include('"%s/DamageTypes.hh"' % srcdir)
    
        # Namespaces.
        SolidSpheral = mod.add_cpp_namespace("Spheral")
        space = SolidSpheral.add_cpp_namespace("PhysicsSpace")

        Spheral = mod.add_cpp_namespace("Spheral")
        PhysicsSpace = Spheral.add_cpp_namespace("PhysicsSpace")

        Physics1d = findObject(PhysicsSpace, "Physics1d")
        Physics2d = findObject(PhysicsSpace, "Physics2d")
        Physics3d = findObject(PhysicsSpace, "Physics3d")

        self.TensorStrainAlgorithm = space.add_enum("TensorStrainAlgorithm", ["BenzAsphaug",
                                                                              "StrainHistory",
                                                                              "MeloshRyanAsphaug",
                                                                              "PlasticStrain",
                                                                              "PseudoPlasticStrain"])
        self.EffectiveDamageAlgorithm = space.add_enum("EffectiveDamageAlgorithm", ["Copy", "Max", "Sampled"])

        self.EffectiveFlawAlgorithm = space.add_enum("EffectiveFlawAlgorithm", ["FullSpectrumFlaws",
                                                                                "MinFlaw",
                                                                                "MaxFlaw",
                                                                                "InverseSumFlaws",
                                                                                "SampledFlaws"])

        self.DamageModel1d = addObject(space, "DamageModel1d", parent=Physics1d, allow_subclassing=True)
        self.DamageModel2d = addObject(space, "DamageModel2d", parent=Physics2d, allow_subclassing=True)
        self.DamageModel3d = addObject(space, "DamageModel3d", parent=Physics3d, allow_subclassing=True)

        self.TensorDamageModel1d = addObject(space, "TensorDamageModel1d", parent=self.DamageModel1d, allow_subclassing=True)
        self.TensorDamageModel2d = addObject(space, "TensorDamageModel2d", parent=self.DamageModel2d, allow_subclassing=True)
        self.TensorDamageModel3d = addObject(space, "TensorDamageModel3d", parent=self.DamageModel3d, allow_subclassing=True)

        self.addWeibullDistributionFunctions(space, 1)
        self.addWeibullDistributionFunctions(space, 2)
        self.addWeibullDistributionFunctions(space, 3)

        self.addComputeFragmentField(SolidSpheral, 1)
        self.addComputeFragmentField(SolidSpheral, 2)
        self.addComputeFragmentField(SolidSpheral, 3)

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        self.generateDamageModelBindings(self.DamageModel1d, 1)
        self.generateDamageModelBindings(self.DamageModel2d, 2)
        self.generateDamageModelBindings(self.DamageModel3d, 3)

        self.generateTensorDamageModelBindings(self.TensorDamageModel1d, 1)
        self.generateTensorDamageModelBindings(self.TensorDamageModel2d, 2)
        self.generateTensorDamageModelBindings(self.TensorDamageModel3d, 3)

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return ["SolidMaterial"]

    #---------------------------------------------------------------------------
    # DamageModel
    #---------------------------------------------------------------------------
    def generateDamageModelBindings(self, x, ndim):

        me = "Spheral::PhysicsSpace::DamageModel%id" % ndim
        dim = "Spheral::Dim<%i> " % ndim
        solidnodelist = "Spheral::SolidMaterial::SolidNodeList%id" % ndim
        vectordoublefield = "Spheral::FieldSpace::VectorDoubleField%id" % ndim
        database = "Spheral::DataBaseSpace::DataBase%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        scalarfieldlist = "Spheral::FieldSpace::ScalarFieldList%id" % ndim
        tablekernel = "Spheral::KernelSpace::TableKernel%id" % ndim
        fileio = "Spheral::FileIOSpace::FileIO"

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
        x.add_method("initialize", None, [param("const double", "time"),
                                          param("const double", "dt"),
                                          constrefparam(database, "dataBase"),
                                          refparam(state, "state"),
                                          refparam(derivatives, "derivatives")], is_virtual=True)
        x.add_method("postStateUpdate", None, [constrefparam(database, "dataBase"),
                                               refparam(state, "state"),
                                               constrefparam(derivatives, "derivatives")], is_const=True, is_virtual=True)
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

        me = "Spheral::PhysicsSpace::TensorDamageModel%id" % ndim
        dim = "Spheral::Dim<%i> " % ndim
        solidnodelist = "Spheral::SolidMaterial::SolidNodeList%id" % ndim
        vectordoublefield = "Spheral::FieldSpace::VectorDoubleField%id" % ndim
        database = "Spheral::DataBaseSpace::DataBase%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        vectorfield = "Spheral::FieldSpace::VectorField%id" % ndim
        symtensorfield = "Spheral::FieldSpace::SymTensorField%id" % ndim
        scalarfieldlist = "Spheral::FieldSpace::ScalarFieldList%id" % ndim
        tablekernel = "Spheral::KernelSpace::TableKernel%id" % ndim
        fileio = "Spheral::FileIOSpace::FileIO"

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
    # weibull flaw distributions
    #---------------------------------------------------------------------------
    def addWeibullDistributionFunctions(self, space, ndim):

        dim = "Spheral::Dim<%i> " % ndim
        fluidnodelist = "Spheral::NodeSpace::FluidNodeList%id" % ndim
        vectordoublefield = "Spheral::FieldSpace::VectorDoubleField%id" % ndim

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
        intfield = "Spheral::FieldSpace::IntField%id" % ndim
        symtensorfield = "Spheral::FieldSpace::SymTensorField%id" % ndim
        nodelist = "Spheral::NodeSpace::NodeList%id" % ndim

        space.add_function("computeFragmentField", intfield,
                           [constrefparam(nodelist, "nodeList"),
                            param("double", "linkRadius"),
                            constrefparam(symtensorfield, "damage"),
                            param("double", "damageThreshold"),
                            param("bool", "assignDustToFragments")],
                           template_parameters = [dim],
                           custom_name = "computeFragmentField%id" % ndim)
        return
