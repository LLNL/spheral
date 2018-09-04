from pybindgen import *

from PBGutils import *
from ref_return_value import *

from PhysicsModule import generatePhysicsVirtualBindings

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class TaylorSPH:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir, dims):

        self.dims = dims

        # Includes.
        mod.add_include('"%s/TaylorSPHTypes.hh"' % srcdir)
    
        # Namespace.
        self.space = mod.add_cpp_namespace("Spheral")
        for dim in self.dims:
            exec('''
generichydro%(dim)id = findObject(self.space, "GenericHydro%(dim)id")
self.TaylorSPHHydroBase%(dim)id = addObject(self.space, "TaylorSPHHydroBase%(dim)id", allow_subclassing=True, parent=generichydro%(dim)id)
''' % {"dim" : dim})

        return

    #---------------------------------------------------------------------------
    # Generate bindings.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        for dim in self.dims:
            exec('''
self.generateTaylorSPHHydroBaseBindings(self.TaylorSPHHydroBase%(dim)id, %(dim)i)
''' % {"dim" : dim})
            self.generateDimBindings(dim)

        return

    #---------------------------------------------------------------------------
    # Add the types per dimension.
    #---------------------------------------------------------------------------
    def generateDimBindings(self, ndim):

        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::SymTensorFieldList%id" % ndim
        connectivitymap = "Spheral::ConnectivityMap%id" % ndim
        tablekernel = "Spheral::TableKernel%id" % ndim

        # TaylorSPH corrections.
        self.space.add_function("computeTaylorSPHCorrections", None,
                                [constrefparam(connectivitymap, "connectivityMap"),
                                 constrefparam(tablekernel, "W"),
                                 constrefparam(scalarfieldlist, "weight"),
                                 constrefparam(vectorfieldlist, "position"),
                                 constrefparam(symtensorfieldlist, "H"),
                                 refparam(tensorfieldlist, "D")],
                                template_parameters = [dim],
                                custom_name = "computeTaylorSPHCorrections%id" % ndim)

        return

    #---------------------------------------------------------------------------
    # Bindings (TaylorSPHHydroBase).
    #---------------------------------------------------------------------------
    def generateTaylorSPHHydroBaseBindings(self, x, ndim):

        # Object names.
        me = "Spheral::TaylorSPHHydroBase%id" % ndim
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
                           constrefparam(tablekernel, "W"),
                           refparam(artificialviscosity, "Q"),
                           param("double", "cfl", default_value="0.5"),
                           param("int", "useVelocityMagnitudeForDt", default_value="false"),
                           param("int", "compatibleEnergyEvolution", default_value="true"),
                           param("int", "XSPH", default_value="true"),
                           param("HEvolutionType", "HUpdate", default_value="Spheral::HEvolutionType::IdealH")])

        # Methods.
        x.add_method("initializeProblemStartup", None, [refparam(database, "dataBase")], is_virtual=True)

        x.add_method("registerState", None, [refparam(database, "dataBase"),
                                             refparam(state, "state")],
                     is_virtual=True)
        x.add_method("registerDerivatives", None, [refparam(database, "dataBase"),
                                                   refparam(derivatives, "derivatives")],
                     is_virtual=True)
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
        x.add_method("label", "std::string", [], is_const=True, is_virtual=True)
        x.add_method("dumpState", None, [refparam(fileio, "fileIO"),
                                         refparam("std::string", "pathName")],
                     is_const = True,
                     is_virtual = True)
        x.add_method("restoreState", None, [refparam(fileio, "fileIO"),
                                            refparam("std::string", "pathName")],
                     is_virtual = True)
        
        # Attributes.
        x.add_instance_attribute("HEvolution", "HEvolutionType", getter="HEvolution", setter="HEvolution")
        x.add_instance_attribute("compatibleEnergyEvolution", "bool", getter="compatibleEnergyEvolution", setter="compatibleEnergyEvolution")
        x.add_instance_attribute("XSPH", "bool", getter="XSPH", setter="XSPH")

        const_ref_return_value(x, me, "%s::smoothingScaleMethod" % me, smoothingscalebase, [], "smoothingScaleMethod")
        const_ref_return_value(x, me, "%s::timeStepMask" % me, intfieldlist, [], "timeStepMask")
        const_ref_return_value(x, me, "%s::pressure" % me, scalarfieldlist, [], "pressure")
        const_ref_return_value(x, me, "%s::soundSpeed" % me, scalarfieldlist, [], "soundSpeed")
        const_ref_return_value(x, me, "%s::specificThermalEnergy0" % me, scalarfieldlist, [], "specificThermalEnergy0")
        const_ref_return_value(x, me, "%s::Hideal" % me, symtensorfieldlist, [], "Hideal")
        const_ref_return_value(x, me, "%s::maxViscousPressure" % me, scalarfieldlist, [], "maxViscousPressure")
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

        const_ref_return_value(x, me, "%s::D" % me, tensorfieldlist, [], "D")

        return
