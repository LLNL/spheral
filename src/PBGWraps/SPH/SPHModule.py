from pybindgen import *

from PBGutils import *
from ref_return_value import *

from PhysicsModule import generatePhysicsVirtualBindings

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class SPH:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir, dims):

        self.dims = dims

        # Includes.
        mod.add_include('"%s/SPHTypes.hh"' % srcdir)
    
        # Namespace.
        self.space = mod.add_cpp_namespace("Spheral")
        
        # Expose types.
        self.NodeCoupling = addObject(self.space, "NodeCoupling", allow_subclassing=True)
        self.PerNodeListNodeCoupling = addObject(self.space, "PerNodeListNodeCoupling", allow_subclassing=True, parent=self.NodeCoupling)

        for dim in self.dims:
            exec('''
generichydro%(dim)id = findObject(self.space, "GenericHydro%(dim)id")

# Expose types.
self.SPHHydroBase%(dim)id = addObject(self.space, "SPHHydroBase%(dim)id", allow_subclassing=True, parent=generichydro%(dim)id)
self.PSPHHydroBase%(dim)id = addObject(self.space, "PSPHHydroBase%(dim)id", allow_subclassing=True, parent=self.SPHHydroBase%(dim)id)
self.DamagedNodeCoupling%(dim)id = addObject(self.space, "DamagedNodeCoupling%(dim)id", allow_subclassing=True, parent=self.NodeCoupling)
self.SolidSPHHydroBase%(dim)id = addObject(self.space, "SolidSPHHydroBase%(dim)id", allow_subclassing=True, parent=self.SPHHydroBase%(dim)id)
''' % {"dim" : dim})

        if 2 in dims:
            self.SPHHydroBaseRZ = addObject(self.space, "SPHHydroBaseRZ", allow_subclassing=True, parent=self.SPHHydroBase2d)
            self.SPHHydroBaseGSRZ = addObject(self.space, "SPHHydroBaseGSRZ", allow_subclassing=True, parent=self.SPHHydroBase2d)
            self.SolidSPHHydroBaseRZ = addObject(self.space, "SolidSPHHydroBaseRZ", allow_subclassing=True, parent=self.SolidSPHHydroBase2d)

        return

    #---------------------------------------------------------------------------
    # Generate bindings.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        self.generateNodeCouplingBindings(self.NodeCoupling)

        for dim in self.dims:
            exec('''
self.generateSPHHydroBaseBindings(self.SPHHydroBase%(dim)id, %(dim)i)
self.generatePSPHHydroBaseBindings(self.PSPHHydroBase%(dim)id, %(dim)i)
self.generateDamagedNodeCouplingBindings(self.DamagedNodeCoupling%(dim)id, %(dim)i)
self.generateSolidSPHHydroBaseBindings(self.SolidSPHHydroBase%(dim)id, %(dim)i)
''' % {"dim" : dim})
            self.generateDimBindings(dim)

        if 2 in self.dims:
            self.generateSPHHydroBaseRZBindings()
            self.generateSPHHydroBaseGSRZBindings()
            self.generateSolidSPHHydroBaseRZBindings()

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return ["SPHSpace"]

    #---------------------------------------------------------------------------
    # Add the types per dimension.
    #---------------------------------------------------------------------------
    def generateDimBindings(self, ndim):

        dim = "Spheral::Dim<%i>" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::SymTensorFieldList%id" % ndim
        connectivitymap = "Spheral::ConnectivityMap%id" % ndim
        tablekernel = "Spheral::TableKernel%id" % ndim

        # SPH sum density.
        self.space.add_function("computeSPHSumMassDensity", None,
                                [constrefparam(connectivitymap, "connectivityMap"),
                                 constrefparam(tablekernel, "W"),
                                 param("bool", "sumMassDensityOverAllNodeLists"),
                                 constrefparam(vectorfieldlist, "position"),
                                 constrefparam(scalarfieldlist, "mass"),
                                 constrefparam(symtensorfieldlist, "H"),
                                 refparam(scalarfieldlist, "massDensity")],
                                template_parameters = [dim],
                                custom_name = "computeSPHSumMassDensity%id" % ndim)

        # SPH grad h correction factor (omega).
        self.space.add_function("computeSPHOmegaGradhCorrection", None,
                                [constrefparam(connectivitymap, "connectivityMap"),
                                 constrefparam(tablekernel, "W"),
                                 constrefparam(vectorfieldlist, "position"),
                                 constrefparam(symtensorfieldlist, "H"),
                                 refparam(scalarfieldlist, "omegaGradh")],
                                template_parameters = [dim],
                                custom_name = "computeSPHOmegaGradhCorrection%id" % ndim)

        return

    #---------------------------------------------------------------------------
    # Bindings (SPHHydroBase).
    #---------------------------------------------------------------------------
    def generateSPHHydroBaseBindings(self, x, ndim):

        # Object names.
        me = "Spheral::SPHHydroBase%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        fieldbase = "Spheral::FieldBase%id" % ndim
        intfield = "Spheral::IntField%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        vector3dfield = "Spheral::Vector3dField%id" % ndim
        tensorfield = "Spheral::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::VectorDoubleField%id" % ndim
        vectorvectorfield = "Spheral::VectorVectorField%id" % ndim
        vectorsymtensorfield = "Spheral::VectorSymTensorField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        vector3dfieldlist = "Spheral::Vector3dFieldList%id" % ndim
        tensorfieldlist = "Spheral::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::ThirdRankTensorFieldList%id" % ndim
        vectordoublefieldlist = "Spheral::VectorDoubleFieldList%id" % ndim
        vectorvectorfieldlist = "Spheral::VectorVectorFieldList%id" % ndim
        vectorsymtensorfieldlist = "Spheral::VectorSymTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        connectivitymap = "Spheral::ConnectivityMap%id" % ndim
        key = "pair_NodeList%id_string" % ndim
        vectorkeys = "vector_of_pair_NodeList%id_string" % ndim
        tablekernel = "Spheral::TableKernel%id" % ndim
        artificialviscosity = "Spheral::ArtificialViscosity%id" % ndim
        fileio = "Spheral::FileIO"
        smoothingscalebase = "Spheral::SmoothingScaleBase%id" % ndim

        # Constructors.
        x.add_constructor([constrefparam(smoothingscalebase, "smoothingScaleMethod"),
                           refparam(artificialviscosity, "Q"),
                           constrefparam(tablekernel, "W"),
                           constrefparam(tablekernel, "WPi"),
                           param("double", "filter", default_value="0.0"),
                           param("double", "cfl", default_value="0.25"),
                           param("int", "useVelocityMagnitudeForDt", default_value="false"),
                           param("int", "compatibleEnergyEvolution", default_value="true"),
                           param("int", "evolveTotalEnergy", default_value="false"),
                           param("int", "gradhCorrection", default_value="false"),
                           param("int", "XSPH", default_value="true"),
                           param("int", "correctVelocityGradient", default_value="false"),
                           param("int", "sumMassDensityOverAllNodeLists", default_value="true"),
                           param("MassDensityType", "densityUpdate", default_value="Spheral::MassDensityType::RigorousSumDensity"),
                           param("HEvolutionType", "HUpdate", default_value="Spheral::HEvolutionType::IdealH"),
                           param("double", "epsTensile", default_value="0.0"),
                           param("double", "nTensile", default_value="4.0"),
                           param(vector, "xmin", default_value="%s(-1e10, -1e10, -1e10)" % vector),
                           param(vector, "xmax", default_value="%s( 1e10,  1e10,  1e10)" % vector)])

        # Methods.
        x.add_method("registerState", None, [refparam(database, "dataBase"),
                                             refparam(state, "state")],
                     is_virtual=True)
        x.add_method("registerDerivatives", None, [refparam(database, "dataBase"),
                                                   refparam(derivatives, "derivatives")],
                     is_virtual=True)
        x.add_method("initializeProblemStartup", None, [refparam(database, "dataBase")], is_virtual=True)
        x.add_method("initialize", None, [param("const double", "time"),
                                          param("const double", "dt"),
                                          constrefparam(database, "dataBase"),
                                          refparam(state, "state"),
                                          refparam(derivatives, "derivatives")], is_virtual=True)
        x.add_method("evaluateDerivatives", None, [param("const double", "time"),
                                                   param("const double", "dt"),
                                                   constrefparam(database, "dataBase"),
                                                   constrefparam(state, "state"),
                                                   refparam(derivatives, "derivatives")],
                     is_const=True, is_virtual=True)
        x.add_method("finalizeDerivatives", None, [param("const double", "time"),
                                                   param("const double", "dt"),
                                                   constrefparam(database, "dataBase"),
                                                   constrefparam(state, "state"),
                                                   refparam(derivatives, "derivatives")], is_const=True, is_virtual=True)
        x.add_method("finalize", None, [param("const double", "time"),
                                        param("const double", "dt"),
                                        refparam(database, "dataBase"),
                                        refparam(state, "state"),
                                        refparam(derivatives, "derivatives")], is_virtual=True)
        x.add_method("applyGhostBoundaries", None, [refparam(state, "state"), refparam(derivatives, "derivatives")], is_virtual=True)
        x.add_method("enforceBoundaries", None, [refparam(state, "state"), refparam(derivatives, "derivatives")], is_virtual=True)
        x.add_method("updateVolume", None, [refparam(state, "state"), param("bool", "boundaries")], is_const=True)
        x.add_method("label", "std::string", [], is_const=True, is_virtual=True)
        x.add_method("dumpState", None, [refparam(fileio, "fileIO"),
                                         param("std::string", "pathName")],
                     is_const = True,
                     is_virtual = True)
        x.add_method("restoreState", None, [constrefparam(fileio, "fileIO"),
                                            param("std::string", "pathName")],
                     is_virtual = True)
        
        # Attributes.
        x.add_instance_attribute("densityUpdate", "MassDensityType", getter="densityUpdate", setter="densityUpdate")
        x.add_instance_attribute("sumForMassDensity", "MassDensityType", getter="densityUpdate", setter="densityUpdate")  # Provided for backwards compatibility
        x.add_instance_attribute("HEvolution", "HEvolutionType", getter="HEvolution", setter="HEvolution")
        x.add_instance_attribute("compatibleEnergyEvolution", "bool", getter="compatibleEnergyEvolution", setter="compatibleEnergyEvolution")
        x.add_instance_attribute("evolveTotalEnergy", "bool", getter="evolveTotalEnergy", setter="evolveTotalEnergy")
        x.add_instance_attribute("gradhCorrection", "bool", getter="gradhCorrection", setter="gradhCorrection")
        x.add_instance_attribute("correctVelocityGradient", "bool", getter="correctVelocityGradient", setter="correctVelocityGradient")
        x.add_instance_attribute("XSPH", "bool", getter="XSPH", setter="XSPH")
        x.add_instance_attribute("sumMassDensityOverAllNodeLists", "bool", getter="sumMassDensityOverAllNodeLists", setter="sumMassDensityOverAllNodeLists")
        x.add_instance_attribute("filter", "double", getter="filter", setter="filter")
        x.add_instance_attribute("epsilonTensile", "double", getter="epsilonTensile", setter="epsilonTensile")
        x.add_instance_attribute("nTensile", "double", getter="nTensile", setter="nTensile")
        x.add_instance_attribute("xmin", vector, getter="xmin", setter="xmin")
        x.add_instance_attribute("xmax", vector, getter="xmax", setter="xmax")

        const_ref_return_value(x, me, "%s::smoothingScaleMethod" % me, smoothingscalebase, [], "smoothingScaleMethod")
        const_ref_return_value(x, me, "%s::timeStepMask" % me, intfieldlist, [], "timeStepMask")
        const_ref_return_value(x, me, "%s::pressure" % me, scalarfieldlist, [], "pressure")
        const_ref_return_value(x, me, "%s::soundSpeed" % me, scalarfieldlist, [], "soundSpeed")
        const_ref_return_value(x, me, "%s::omegaGradh" % me, scalarfieldlist, [], "omegaGradh")
        const_ref_return_value(x, me, "%s::specificThermalEnergy0" % me, scalarfieldlist, [], "specificThermalEnergy0")
        const_ref_return_value(x, me, "%s::entropy" % me, scalarfieldlist, [], "entropy")
        const_ref_return_value(x, me, "%s::Hideal" % me, symtensorfieldlist, [], "Hideal")
        const_ref_return_value(x, me, "%s::maxViscousPressure" % me, scalarfieldlist, [], "maxViscousPressure")
        const_ref_return_value(x, me, "%s::effectiveViscousPressure" % me, scalarfieldlist, [], "effectiveViscousPressure")
        const_ref_return_value(x, me, "%s::viscousWork" % me, scalarfieldlist, [], "viscousWork")
        const_ref_return_value(x, me, "%s::massDensitySum" % me, scalarfieldlist, [], "massDensitySum")
        const_ref_return_value(x, me, "%s::weightedNeighborSum" % me, scalarfieldlist, [], "weightedNeighborSum")
        const_ref_return_value(x, me, "%s::massSecondMoment" % me, symtensorfieldlist, [], "massSecondMoment")
        const_ref_return_value(x, me, "%s::XSPHWeightSum" % me, scalarfieldlist, [], "XSPHWeightSum")
        const_ref_return_value(x, me, "%s::XSPHDeltaV" % me, vectorfieldlist, [], "XSPHDeltaV")
        const_ref_return_value(x, me, "%s::DxDt" % me, vectorfieldlist, [], "DxDt")
        const_ref_return_value(x, me, "%s::DvDt" % me, vectorfieldlist, [], "DvDt")
        const_ref_return_value(x, me, "%s::DmassDensityDt" % me, scalarfieldlist, [], "DmassDensityDt")
        const_ref_return_value(x, me, "%s::DspecificThermalEnergyDt" % me, scalarfieldlist, [], "DspecificThermalEnergyDt")
        const_ref_return_value(x, me, "%s::DHDt" % me, symtensorfieldlist, [], "DHDt")
        const_ref_return_value(x, me, "%s::DvDx" % me, tensorfieldlist, [], "DvDx")
        const_ref_return_value(x, me, "%s::internalDvDx" % me, tensorfieldlist, [], "internalDvDx")
        const_ref_return_value(x, me, "%s::M" % me, tensorfieldlist, [], "M")
        const_ref_return_value(x, me, "%s::localM" % me, tensorfieldlist, [], "localM")
        const_ref_return_value(x, me, "%s::pairAccelerations" % me, vectorvectorfieldlist, [], "pairAccelerations")

        return

    #---------------------------------------------------------------------------
    # Bindings (PSPHHydroBase).
    #---------------------------------------------------------------------------
    def generatePSPHHydroBaseBindings(self, x, ndim):

        # Object names.
        me = "Spheral::PSPHHydroBase%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        fieldbase = "Spheral::FieldBase%id" % ndim
        intfield = "Spheral::IntField%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        vector3dfield = "Spheral::Vector3dField%id" % ndim
        tensorfield = "Spheral::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::VectorDoubleField%id" % ndim
        vectorvectorfield = "Spheral::VectorVectorField%id" % ndim
        vectorsymtensorfield = "Spheral::VectorSymTensorField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        vector3dfieldlist = "Spheral::Vector3dFieldList%id" % ndim
        tensorfieldlist = "Spheral::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::ThirdRankTensorFieldList%id" % ndim
        vectordoublefieldlist = "Spheral::VectorDoubleFieldList%id" % ndim
        vectorvectorfieldlist = "Spheral::VectorVectorFieldList%id" % ndim
        vectorsymtensorfieldlist = "Spheral::VectorSymTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        connectivitymap = "Spheral::ConnectivityMap%id" % ndim
        key = "pair_NodeList%id_string" % ndim
        vectorkeys = "vector_of_pair_NodeList%id_string" % ndim
        tablekernel = "Spheral::TableKernel%id" % ndim
        artificialviscosity = "Spheral::ArtificialViscosity%id" % ndim
        fileio = "Spheral::FileIO"
        smoothingscalebase = "Spheral::SmoothingScaleBase%id" % ndim

        # Constructors.
        x.add_constructor([constrefparam(smoothingscalebase, "smoothingScaleMethod"),
                           refparam(artificialviscosity, "Q"),
                           constrefparam(tablekernel, "W"),
                           constrefparam(tablekernel, "WPi"),
                           param("double", "filter", default_value="0.0"),
                           param("double", "cfl", default_value="0.25"),
                           param("int", "useVelocityMagnitudeForDt", default_value="false"),
                           param("int", "compatibleEnergyEvolution", default_value="true"),
                           param("int", "evolveTotalEnergy", default_value="false"),
                           param("int", "XSPH", default_value="true"),
                           param("int", "correctVelocityGradient", default_value="false"),
                           param("int", "HopkinsConductivity", default_value="false"),
                           param("int", "sumMassDensityOverAllNodeLists", default_value="true"),
                           param("MassDensityType", "densityUpdate", default_value="Spheral::MassDensityType::RigorousSumDensity"),
                           param("HEvolutionType", "HUpdate", default_value="Spheral::HEvolutionType::IdealH"),
                           param(vector, "xmin", default_value="%s(-1e10, -1e10, -1e10)" % vector),
                           param(vector, "xmax", default_value="%s( 1e10,  1e10,  1e10)" % vector)])

        # Attributes.
        x.add_instance_attribute("HopkinsConductivity", "bool", getter="HopkinsConductivity", setter="HopkinsConductivity")

        const_ref_return_value(x, me, "%s::gamma" % me, scalarfieldlist, [], "gamma")
        const_ref_return_value(x, me, "%s::PSPHcorrection" % me, scalarfieldlist, [], "PSPHcorrection")

        return

    #---------------------------------------------------------------------------
    # Bindings (SPHHydroBaseRZ).
    #---------------------------------------------------------------------------
    def generateSPHHydroBaseRZBindings(self):

        # Object names.
        x = self.SPHHydroBaseRZ
        ndim = 2
        me = "Spheral::SPHHydroBase%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        fieldbase = "Spheral::FieldBase%id" % ndim
        intfield = "Spheral::IntField%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        vector3dfield = "Spheral::Vector3dField%id" % ndim
        tensorfield = "Spheral::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::VectorDoubleField%id" % ndim
        vectorvectorfield = "Spheral::VectorVectorField%id" % ndim
        vectorsymtensorfield = "Spheral::VectorSymTensorField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        vector3dfieldlist = "Spheral::Vector3dFieldList%id" % ndim
        tensorfieldlist = "Spheral::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::ThirdRankTensorFieldList%id" % ndim
        vectordoublefieldlist = "Spheral::VectorDoubleFieldList%id" % ndim
        vectorvectorfieldlist = "Spheral::VectorVectorFieldList%id" % ndim
        vectorsymtensorfieldlist = "Spheral::VectorSymTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        connectivitymap = "Spheral::ConnectivityMap%id" % ndim
        key = "pair_NodeList%id_string" % ndim
        vectorkeys = "vector_of_pair_NodeList%id_string" % ndim
        tablekernel = "Spheral::TableKernel%id" % ndim
        artificialviscosity = "Spheral::ArtificialViscosity%id" % ndim
        fileio = "Spheral::FileIO"
        smoothingscalebase = "Spheral::SmoothingScaleBase%id" % ndim

        # Constructors.
        x.add_constructor([constrefparam(smoothingscalebase, "smoothingScaleMethod"),
                           refparam(artificialviscosity, "Q"),
                           constrefparam(tablekernel, "W"),
                           constrefparam(tablekernel, "WPi"),
                           param("double", "filter", default_value="0.0"),
                           param("double", "cfl", default_value="0.25"),
                           param("int", "useVelocityMagnitudeForDt", default_value="false"),
                           param("int", "compatibleEnergyEvolution", default_value="true"),
                           param("int", "evolveTotalEnergy", default_value="false"),
                           param("int", "gradhCorrection", default_value="false"),
                           param("int", "XSPH", default_value="true"),
                           param("int", "correctVelocityGradient", default_value="false"),
                           param("int", "sumMassDensityOverAllNodeLists", default_value="true"),
                           param("MassDensityType", "densityUpdate", default_value="Spheral::MassDensityType::RigorousSumDensity"),
                           param("HEvolutionType", "HUpdate", default_value="Spheral::HEvolutionType::IdealH"),
                           param("double", "epsTensile", default_value="0.0"),
                           param("double", "nTensile", default_value="4.0"),
                           param(vector, "xmin", default_value="%s(-1e10, -1e10, -1e10)" % vector),
                           param(vector, "xmax", default_value="%s( 1e10,  1e10,  1e10)" % vector)])

        return

    #---------------------------------------------------------------------------
    # Bindings (SPHHydroBaseGSRZ).
    #---------------------------------------------------------------------------
    def generateSPHHydroBaseGSRZBindings(self):

        # Object names.
        x = self.SPHHydroBaseGSRZ
        ndim = 2
        me = "Spheral::SPHHydroBase%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        fieldbase = "Spheral::FieldBase%id" % ndim
        intfield = "Spheral::IntField%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        vector3dfield = "Spheral::Vector3dField%id" % ndim
        tensorfield = "Spheral::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::VectorDoubleField%id" % ndim
        vectorvectorfield = "Spheral::VectorVectorField%id" % ndim
        vectorsymtensorfield = "Spheral::VectorSymTensorField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        vector3dfieldlist = "Spheral::Vector3dFieldList%id" % ndim
        tensorfieldlist = "Spheral::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::ThirdRankTensorFieldList%id" % ndim
        vectordoublefieldlist = "Spheral::VectorDoubleFieldList%id" % ndim
        vectorvectorfieldlist = "Spheral::VectorVectorFieldList%id" % ndim
        vectorsymtensorfieldlist = "Spheral::VectorSymTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        connectivitymap = "Spheral::ConnectivityMap%id" % ndim
        key = "pair_NodeList%id_string" % ndim
        vectorkeys = "vector_of_pair_NodeList%id_string" % ndim
        tablekernel = "Spheral::TableKernel%id" % ndim
        artificialviscosity = "Spheral::ArtificialViscosity%id" % ndim
        fileio = "Spheral::FileIO"
        smoothingscalebase = "Spheral::SmoothingScaleBase%id" % ndim

        # Constructors.
        x.add_constructor([constrefparam(smoothingscalebase, "smoothingScaleMethod"),
                           refparam(artificialviscosity, "Q"),
                           constrefparam(tablekernel, "W"),
                           constrefparam(tablekernel, "WPi"),
                           param("double", "filter", default_value="0.0"),
                           param("double", "cfl", default_value="0.25"),
                           param("int", "useVelocityMagnitudeForDt", default_value="false"),
                           param("int", "compatibleEnergyEvolution", default_value="true"),
                           param("int", "evolveTotalEnergy", default_value="false"),
                           param("int", "gradhCorrection", default_value="false"),
                           param("int", "XSPH", default_value="true"),
                           param("int", "correctVelocityGradient", default_value="false"),
                           param("int", "sumMassDensityOverAllNodeLists", default_value="true"),
                           param("MassDensityType", "densityUpdate", default_value="Spheral::MassDensityType::RigorousSumDensity"),
                           param("HEvolutionType", "HUpdate", default_value="Spheral::HEvolutionType::IdealH"),
                           param("double", "epsTensile", default_value="0.0"),
                           param("double", "nTensile", default_value="4.0"),
                           param(vector, "xmin", default_value="%s(-1e10, -1e10, -1e10)" % vector),
                           param(vector, "xmax", default_value="%s( 1e10,  1e10,  1e10)" % vector)])

        return

    #---------------------------------------------------------------------------
    # Bindings (NodeCoupling).
    #---------------------------------------------------------------------------
    def generateNodeCouplingBindings(self, x):

        # Object names.
        me = "Spheral::NodeCoupling"

        # Constructors.
        x.add_constructor([])

        # Methods.
        x.add_method("operator()", "double", [param("unsigned int", "nodeListi"),
                                              param("unsigned int", "i"),
                                              param("unsigned int", "nodeListj"),
                                              param("unsigned int", "j")],
                     custom_name = "__call__",
                     is_const = True,
                     is_virtual = True)
        return

    #---------------------------------------------------------------------------
    # Bindings (DamagedNodeCoupling).
    #---------------------------------------------------------------------------
    def generateDamagedNodeCouplingBindings(self, x, ndim):

        # Object names.
        me = "Spheral::DamagedNodeCoupling%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        vector3dfieldlist = "Spheral::Vector3dFieldList%id" % ndim
        tensorfieldlist = "Spheral::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::SymTensorFieldList%id" % ndim

        # Constructors.
        x.add_constructor([constrefparam(symtensorfieldlist, "damage"),
                           constrefparam(vectorfieldlist, "damageGradient"),
                           constrefparam(symtensorfieldlist, "H")])

        x.add_method("operator()", "double",
                     [param("unsigned int", "nodeListi"),
                      param("unsigned int", "i"),
                      param("unsigned int", "nodeListj"),
                      param("unsigned int", "j")],
                     custom_name = "__call__")

        return

    #---------------------------------------------------------------------------
    # Bindings (SolidSPHHydroBase).
    #---------------------------------------------------------------------------
    def generateSolidSPHHydroBaseBindings(self, x, ndim):

        # Object names.
        me = "Spheral::SolidSPHHydroBase%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        fieldbase = "Spheral::FieldBase%id" % ndim
        intfield = "Spheral::IntField%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        vector3dfield = "Spheral::Vector3dField%id" % ndim
        tensorfield = "Spheral::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::VectorDoubleField%id" % ndim
        vectorvectorfield = "Spheral::VectorVectorField%id" % ndim
        vectorsymtensorfield = "Spheral::VectorSymTensorField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        vector3dfieldlist = "Spheral::Vector3dFieldList%id" % ndim
        tensorfieldlist = "Spheral::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::ThirdRankTensorFieldList%id" % ndim
        vectordoublefieldlist = "Spheral::VectorDoubleFieldList%id" % ndim
        vectorvectorfieldlist = "Spheral::VectorVectorFieldList%id" % ndim
        vectorsymtensorfieldlist = "Spheral::VectorSymTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        connectivitymap = "Spheral::ConnectivityMap%id" % ndim
        key = "pair_NodeList%id_string" % ndim
        vectorkeys = "vector_of_pair_NodeList%id_string" % ndim
        tablekernel = "Spheral::TableKernel%id" % ndim
        artificialviscosity = "Spheral::ArtificialViscosity%id" % ndim
        fileio = "Spheral::FileIO"
        smoothingscalebase = "Spheral::SmoothingScaleBase%id" % ndim

        # Constructors.
        x.add_constructor([constrefparam(smoothingscalebase, "smoothingScaleMethod"),
                           refparam(artificialviscosity, "Q"),
                           constrefparam(tablekernel, "W"),
                           constrefparam(tablekernel, "WPi"),
                           constrefparam(tablekernel, "WGrad"),
                           param("double", "filter", default_value="0.0"),
                           param("double", "cfl", default_value="0.25"),
                           param("int", "useVelocityMagnitudeForDt", default_value="false"),
                           param("int", "compatibleEnergyEvolution", default_value="true"),
                           param("int", "evolveTotalEnergy", default_value="false"),
                           param("int", "gradhCorrection", default_value="false"),
                           param("int", "XSPH", default_value="true"),
                           param("int", "correctVelocityGradient", default_value="true"),
                           param("int", "sumMassDensityOverAllNodeLists", default_value="false"),
                           param("MassDensityType", "densityUpdate", default_value="Spheral::MassDensityType::RigorousSumDensity"),
                           param("HEvolutionType", "HUpdate", default_value="Spheral::HEvolutionType::IdealH"),
                           param("double", "epsTensile", default_value="0.0"),
                           param("double", "nTensile", default_value="4.0"),
                           param("bool", "damageRelieveRubble", default_value="true"),
                           param(vector, "xmin", default_value="%s(-1e10, -1e10, -1e10)" % vector),
                           param(vector, "xmax", default_value="%s( 1e10,  1e10,  1e10)" % vector)])

        # Methods.
        x.add_method("initializeProblemStartup", None, [refparam(database, "dataBase")], is_virtual=True)
        x.add_method("registerState", None, [refparam(database, "dataBase"),
                                             refparam(state, "state")],
                     is_virtual=True)
        x.add_method("registerDerivatives", None, [refparam(database, "dataBase"),
                                                   refparam(derivatives, "derivatives")],
                     is_virtual=True)
        x.add_method("evaluateDerivatives", None, [param("const double", "time"),
                                                   param("const double", "dt"),
                                                   constrefparam(database, "dataBase"),
                                                   constrefparam(state, "state"),
                                                   refparam(derivatives, "derivatives")],
                     is_const=True, is_virtual=True)
        x.add_method("applyGhostBoundaries", None, [refparam(state, "state"), refparam(derivatives, "derivatives")], is_virtual=True)
        x.add_method("enforceBoundaries", None, [refparam(state, "state"), refparam(derivatives, "derivatives")], is_virtual=True)
        x.add_method("label", "std::string", [], is_const=True, is_virtual=True)
        x.add_method("dumpState", None, [refparam(fileio, "fileIO"),
                                         refparam("std::string", "pathName")],
                     is_const = True,
                     is_virtual = True)
        x.add_method("restoreState", None, [refparam(fileio, "fileIO"),
                                            refparam("std::string", "pathName")],
                     is_virtual = True)
        
        # Attributes.
        const_ref_return_value(x, me, "%s::DdeviatoricStressDt" % me, symtensorfieldlist, [], "DdeviatoricStressDt")
        const_ref_return_value(x, me, "%s::bulkModulus" % me, scalarfieldlist, [], "bulkModulus")
        const_ref_return_value(x, me, "%s::shearModulus" % me, scalarfieldlist, [], "shearModulus")
        const_ref_return_value(x, me, "%s::yieldStrength" % me, scalarfieldlist, [], "yieldStrength")
        const_ref_return_value(x, me, "%s::plasticStrain0" % me, scalarfieldlist, [], "plasticStrain0")
        const_ref_return_value(x, me, "%s::Hfield0" % me, symtensorfieldlist, [], "Hfield0")
        x.add_instance_attribute("damageRelieveRubble", "bool", getter="damageRelieveRubble", setter="damageRelieveRubble")

        return

    #---------------------------------------------------------------------------
    # Bindings (SolidSPHHydroBaseRZ).
    #---------------------------------------------------------------------------
    def generateSolidSPHHydroBaseRZBindings(self):

        # Object names.
        x = self.SolidSPHHydroBaseRZ
        ndim = 2
        me = "Spheral::SolidSPHHydroBaseRZ"
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        fieldbase = "Spheral::FieldBase%id" % ndim
        intfield = "Spheral::IntField%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        vector3dfield = "Spheral::Vector3dField%id" % ndim
        tensorfield = "Spheral::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::VectorDoubleField%id" % ndim
        vectorvectorfield = "Spheral::VectorVectorField%id" % ndim
        vectorsymtensorfield = "Spheral::VectorSymTensorField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        vector3dfieldlist = "Spheral::Vector3dFieldList%id" % ndim
        tensorfieldlist = "Spheral::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::ThirdRankTensorFieldList%id" % ndim
        vectordoublefieldlist = "Spheral::VectorDoubleFieldList%id" % ndim
        vectorvectorfieldlist = "Spheral::VectorVectorFieldList%id" % ndim
        vectorsymtensorfieldlist = "Spheral::VectorSymTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        connectivitymap = "Spheral::ConnectivityMap%id" % ndim
        key = "pair_NodeList%id_string" % ndim
        vectorkeys = "vector_of_pair_NodeList%id_string" % ndim
        tablekernel = "Spheral::TableKernel%id" % ndim
        artificialviscosity = "Spheral::ArtificialViscosity%id" % ndim
        fileio = "Spheral::FileIO"
        smoothingscalebase = "Spheral::SmoothingScaleBase%id" % ndim

        # Constructors.
        x.add_constructor([constrefparam(smoothingscalebase, "smoothingScaleMethod"),
                           refparam(artificialviscosity, "Q"),
                           constrefparam(tablekernel, "W"),
                           constrefparam(tablekernel, "WPi"),
                           constrefparam(tablekernel, "WGrad"),
                           param("double", "filter", default_value="0.0"),
                           param("double", "cfl", default_value="0.25"),
                           param("int", "useVelocityMagnitudeForDt", default_value="false"),
                           param("int", "compatibleEnergyEvolution", default_value="true"),
                           param("int", "evolveTotalEnergy", default_value="false"),
                           param("int", "gradhCorrection", default_value="false"),
                           param("int", "XSPH", default_value="true"),
                           param("int", "correctVelocityGradient", default_value="true"),
                           param("int", "sumMassDensityOverAllNodeLists", default_value="false"),
                           param("MassDensityType", "densityUpdate", default_value="Spheral::MassDensityType::RigorousSumDensity"),
                           param("HEvolutionType", "HUpdate", default_value="Spheral::HEvolutionType::IdealH"),
                           param("double", "epsTensile", default_value="0.0"),
                           param("double", "nTensile", default_value="4.0"),
                           param("bool", "damageRelieveRubble", default_value="false"),
                           param(vector, "xmin", default_value="%s(-1e10, -1e10, -1e10)" % vector),
                           param(vector, "xmax", default_value="%s( 1e10,  1e10,  1e10)" % vector)])

        # Attributes.
        const_ref_return_value(x, me, "%s::deviatoricStressTT" % me, scalarfieldlist, [], "deviatoricStressTT")
        const_ref_return_value(x, me, "%s::DdeviatoricStressTTDt" % me, scalarfieldlist, [], "DdeviatoricStressTTDt")

        return
