from pybindgen import *

import sys
sys.path.append("..")
from PBGutils import *
from ref_return_value import *

sys.path.append("../Physics")
from PhysicsModule import generatePhysicsVirtualBindings

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class SPH:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod):

        # Includes.
        mod.add_include('"SPH/SPHTypes.hh"')
    
        # Namespace.
        Spheral = mod.add_cpp_namespace("Spheral")
        PhysicsSpace = Spheral.add_cpp_namespace("PhysicsSpace")
        self.space = Spheral.add_cpp_namespace("SPHSpace")
        
        generichydro1d = findObject(PhysicsSpace, "GenericHydro1d")
        generichydro2d = findObject(PhysicsSpace, "GenericHydro2d")
        generichydro3d = findObject(PhysicsSpace, "GenericHydro3d")

        # Expose types.
        self.SPHHydroBase1d = addObject(self.space, "SPHHydroBase1d", allow_subclassing=True, parent=generichydro1d)
        self.SPHHydroBase2d = addObject(self.space, "SPHHydroBase2d", allow_subclassing=True, parent=generichydro2d)
        self.SPHHydroBase3d = addObject(self.space, "SPHHydroBase3d", allow_subclassing=True, parent=generichydro3d)

        return

    #---------------------------------------------------------------------------
    # Generate bindings.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        self.generateSPHHydroBaseBindings(self.SPHHydroBase1d, 1)
        self.generateSPHHydroBaseBindings(self.SPHHydroBase2d, 2)
        self.generateSPHHydroBaseBindings(self.SPHHydroBase3d, 3)
        
        self.generateDimBindings(1)
        self.generateDimBindings(2)
        self.generateDimBindings(3)

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
        scalarfieldlist = "Spheral::FieldSpace::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::FieldSpace::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::FieldSpace::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::FieldSpace::SymTensorFieldList%id" % ndim
        connectivitymap = "Spheral::NeighborSpace::ConnectivityMap%id" % ndim
        tablekernel = "Spheral::KernelSpace::TableKernel%id" % ndim

        # SPH sum density.
        self.space.add_function("computeSPHSumMassDensity", None,
                                [constrefparam(connectivitymap, "connectivityMap"),
                                 constrefparam(tablekernel, "W"),
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
        me = "Spheral::SPHSpace::SPHHydroBase%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        fieldbase = "Spheral::FieldSpace::FieldBase%id" % ndim
        intfield = "Spheral::FieldSpace::IntField%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        vectorfield = "Spheral::FieldSpace::VectorField%id" % ndim
        vector3dfield = "Spheral::FieldSpace::Vector3dField%id" % ndim
        tensorfield = "Spheral::FieldSpace::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::FieldSpace::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::FieldSpace::VectorDoubleField%id" % ndim
        vectorvectorfield = "Spheral::FieldSpace::VectorVectorField%id" % ndim
        vectorsymtensorfield = "Spheral::FieldSpace::VectorSymTensorField%id" % ndim
        symtensorfield = "Spheral::FieldSpace::SymTensorField%id" % ndim
        intfieldlist = "Spheral::FieldSpace::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::FieldSpace::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::FieldSpace::VectorFieldList%id" % ndim
        vector3dfieldlist = "Spheral::FieldSpace::Vector3dFieldList%id" % ndim
        tensorfieldlist = "Spheral::FieldSpace::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::FieldSpace::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::FieldSpace::ThirdRankTensorFieldList%id" % ndim
        vectordoublefieldlist = "Spheral::FieldSpace::VectorDoubleFieldList%id" % ndim
        vectorvectorfieldlist = "Spheral::FieldSpace::VectorVectorFieldList%id" % ndim
        vectorsymtensorfieldlist = "Spheral::FieldSpace::VectorSymTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeSpace::NodeList%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        database = "Spheral::DataBaseSpace::DataBase%id" % ndim
        connectivitymap = "Spheral::NeighborSpace::ConnectivityMap%id" % ndim
        key = "pair_NodeList%id_string" % ndim
        vectorkeys = "vector_of_pair_NodeList%id_string" % ndim
        tablekernel = "Spheral::KernelSpace::TableKernel%id" % ndim
        artificialviscosity = "Spheral::ArtificialViscositySpace::ArtificialViscosity%id" % ndim
        fileio = "Spheral::FileIOSpace::FileIO"
        smoothingscalebase = "Spheral::NodeSpace::SmoothingScaleBase%id" % ndim

        # Constructors.
        x.add_constructor([constrefparam(smoothingscalebase, "smoothingScaleMethod"),
                           constrefparam(tablekernel, "W"),
                           constrefparam(tablekernel, "WPi"),
                           refparam(artificialviscosity, "Q"),
                           param("double", "cfl", default_value="0.5"),
                           param("int", "useVelocityMagnitudeForDt", default_value="false"),
                           param("int", "compatibleEnergyEvolution", default_value="true"),
                           param("int", "gradhCorrection", default_value="false"),
                           param("int", "XSPH", default_value="true"),
                           param("MassDensityType", "densityUpdate", default_value="Spheral::PhysicsSpace::RigorousSumDensity"),
                           param("HEvolutionType", "HUpdate", default_value="Spheral::PhysicsSpace::IdealH"),
                           param("double", "epsTensile", default_value="0.3"),
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
                                         refparam("std::string", "pathName")],
                     is_const = True,
                     is_virtual = True)
        x.add_method("restoreState", None, [refparam(fileio, "fileIO"),
                                            refparam("std::string", "pathName")],
                     is_virtual = True)
        
        # Attributes.
        x.add_instance_attribute("densityUpdate", "MassDensityType", getter="densityUpdate", setter="densityUpdate")
        x.add_instance_attribute("sumForMassDensity", "MassDensityType", getter="densityUpdate", setter="densityUpdate")  # Provided for backwards compatibility
        x.add_instance_attribute("HEvolution", "HEvolutionType", getter="HEvolution", setter="HEvolution")
        x.add_instance_attribute("compatibleEnergyEvolution", "bool", getter="compatibleEnergyEvolution", setter="compatibleEnergyEvolution")
        x.add_instance_attribute("gradhCorrection", "bool", getter="gradhCorrection", setter="gradhCorrection")
        x.add_instance_attribute("XSPH", "bool", getter="XSPH", setter="XSPH")
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
        const_ref_return_value(x, me, "%s::Hideal" % me, symtensorfieldlist, [], "Hideal")
        const_ref_return_value(x, me, "%s::maxViscousPressure" % me, scalarfieldlist, [], "maxViscousPressure")
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
        const_ref_return_value(x, me, "%s::pairAccelerations" % me, vectorvectorfieldlist, [], "pairAccelerations")

        return

    # #---------------------------------------------------------------------------
    # # Bindings (TotalHydro).
    # #---------------------------------------------------------------------------
    # def generateTotalHydroBindings(self, x, ndim):

    #     # Object names.
    #     me = "Spheral::PhysicsSpace::TotalHydro%id" % ndim
    #     dim = "Spheral::Dim<%i>" % ndim
    #     vector = "Vector%id" % ndim
    #     tensor = "Tensor%id" % ndim
    #     symtensor = "SymTensor%id" % ndim
    #     fieldbase = "Spheral::FieldSpace::FieldBase%id" % ndim
    #     intfield = "Spheral::FieldSpace::IntField%id" % ndim
    #     scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
    #     vectorfield = "Spheral::FieldSpace::VectorField%id" % ndim
    #     vector3dfield = "Spheral::FieldSpace::Vector3dField%id" % ndim
    #     tensorfield = "Spheral::FieldSpace::TensorField%id" % ndim
    #     thirdranktensorfield = "Spheral::FieldSpace::ThirdRankTensorField%id" % ndim
    #     vectordoublefield = "Spheral::FieldSpace::VectorDoubleField%id" % ndim
    #     vectorvectorfield = "Spheral::FieldSpace::VectorVectorField%id" % ndim
    #     vectorsymtensorfield = "Spheral::FieldSpace::VectorSymTensorField%id" % ndim
    #     symtensorfield = "Spheral::FieldSpace::SymTensorField%id" % ndim
    #     intfieldlist = "Spheral::FieldSpace::IntFieldList%id" % ndim
    #     scalarfieldlist = "Spheral::FieldSpace::ScalarFieldList%id" % ndim
    #     vectorfieldlist = "Spheral::FieldSpace::VectorFieldList%id" % ndim
    #     vector3dfieldlist = "Spheral::FieldSpace::Vector3dFieldList%id" % ndim
    #     tensorfieldlist = "Spheral::FieldSpace::TensorFieldList%id" % ndim
    #     symtensorfieldlist = "Spheral::FieldSpace::SymTensorFieldList%id" % ndim
    #     thirdranktensorfieldlist = "Spheral::FieldSpace::ThirdRankTensorFieldList%id" % ndim
    #     vectordoublefieldlist = "Spheral::FieldSpace::VectorDoubleFieldList%id" % ndim
    #     vectorvectorfieldlist = "Spheral::FieldSpace::VectorVectorFieldList%id" % ndim
    #     vectorsymtensorfieldlist = "Spheral::FieldSpace::VectorSymTensorFieldList%id" % ndim
    #     nodelist = "Spheral::NodeSpace::NodeList%id" % ndim
    #     state = "Spheral::State%id" % ndim
    #     derivatives = "Spheral::StateDerivatives%id" % ndim
    #     database = "Spheral::DataBaseSpace::DataBase%id" % ndim
    #     connectivitymap = "Spheral::NeighborSpace::ConnectivityMap%id" % ndim
    #     key = "pair_NodeList%id_string" % ndim
    #     vectorkeys = "vector_of_pair_NodeList%id_string" % ndim
    #     tablekernel = "Spheral::KernelSpace::TableKernel%id" % ndim
    #     artificialviscosity = "Spheral::ArtificialViscositySpace::ArtificialViscosity%id" % ndim
    #     fileio = "Spheral::FileIOSpace::FileIO"

    #     # Constructors.
    #     x.add_constructor([constrefparam(tablekernel, "W"),
    #                        constrefparam(tablekernel, "WPi"),
    #                        refparam(artificialviscosity, "Q"),
    #                        param("HEvolutionType", "HUpdate", default_value="Spheral::PhysicsSpace::IdealH"),
    #                        param("double", "hmin", default_value="1.0e-100"),
    #                        param("double", "hmax", default_value="1.0e100"),
    #                        param("double", "hratiomin", default_value="0.1")])

    #     # Methods.
    #     x.add_method("evaluateDerivatives", None, [param("const double", "time"),
    #                                                param("const double", "dt"),
    #                                                constrefparam(database, "dataBase"),
    #                                                constrefparam(state, "state"),
    #                                                refparam(derivatives, "derivatives")],
    #                  is_const=True, is_virtual=True)
    #     x.add_method("postStateUpdate", None, [constrefparam(database, "dataBase"),
    #                                            refparam(state, "state"),
    #                                            constrefparam(derivatives, "derivatives")], is_const=True, is_virtual=True)
    #     x.add_method("registerState", None, [refparam(database, "dataBase"),
    #                                          refparam(state, "state")],
    #                  is_virtual=True)
    #     x.add_method("registerDerivatives", None, [refparam(database, "dataBase"),
    #                                                refparam(derivatives, "derivatives")],
    #                  is_virtual=True)
    #     x.add_method("applyGhostBoundaries", None, [refparam(state, "state"), refparam(derivatives, "derivatives")], is_virtual=True)
    #     x.add_method("enforceBoundaries", None, [refparam(state, "state"), refparam(derivatives, "derivatives")], is_virtual=True)
    #     x.add_method("label", "std::string", [], is_const=True, is_virtual=True)
    #     x.add_method("dumpState", None, [refparam(fileio, "fileIO"),
    #                                      refparam("std::string", "pathName")],
    #                  is_const = True,
    #                  is_virtual = True)
    #     x.add_method("restoreState", None, [refparam(fileio, "fileIO"),
    #                                         refparam("std::string", "pathName")],
    #                  is_virtual = True)
        
    #     # Attributes.
    #     x.add_instance_attribute("HEvolution", "HEvolutionType", getter="HEvolution", setter="HEvolution")
    #     x.add_instance_attribute("hmin", "double", getter="hmin", setter="hmin")
    #     x.add_instance_attribute("hmax", "double", getter="hmax", setter="hmax")
    #     x.add_instance_attribute("hratiomin", "double", getter="hratiomin", setter="hratiomin")
    #     x.add_instance_attribute("Hideal", symtensorfieldlist, getter="Hideal", is_const=True)
    #     x.add_instance_attribute("timeStepMask", intfieldlist, getter="timeStepMask", is_const=True)
    #     x.add_instance_attribute("pressure", scalarfieldlist, getter="pressure", is_const=True)
    #     x.add_instance_attribute("soundSpeed", scalarfieldlist, getter="soundSpeed", is_const=True)
    #     x.add_instance_attribute("weightedNeighborSum", scalarfieldlist, getter="weightedNeighborSum", is_const=True)
    #     x.add_instance_attribute("totalEnergy", scalarfieldlist, getter="totalEnergy", is_const=True)
    #     x.add_instance_attribute("volume", scalarfieldlist, getter="volume", is_const=True)
    #     x.add_instance_attribute("linearMomentum", vectorfieldlist, getter="linearMomentum", is_const=True)
    #     x.add_instance_attribute("DEDt", scalarfieldlist, getter="DEDt", is_const=True)
    #     x.add_instance_attribute("massSecondMoment", symtensorfieldlist, getter="massSecondMoment", is_const=True)

    #     return
