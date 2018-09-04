from pybindgen import *

from PBGutils import *
from ref_return_value import *

from PhysicsModule import generatePhysicsVirtualBindings

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class SVPH:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir, dims):

        self.dims = dims

        # Includes.
        mod.add_include('"%s/SVPHTypes.hh"' % srcdir)
    
        # Namespace.
        self.space = mod.add_cpp_namespace("Spheral")
        
        for dim in self.dims:
            exec('''
generichydro%(dim)id = findObject(self.space, "GenericHydro%(dim)id")
self.SVPHFacetedHydroBase%(dim)id = addObject(self.space, "SVPHFacetedHydroBase%(dim)id", allow_subclassing=True, parent=generichydro%(dim)id)
''' % {"dim" : dim})

        return

    #---------------------------------------------------------------------------
    # Generate bindings.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        for dim in self.dims:
            exec('''
self.generateSVPHFacetedHydroBaseBindings(self.SVPHFacetedHydroBase%(dim)id, %(dim)i)
''' % {"dim" : dim})
            self.generateDimBindings(dim)

        return

    #---------------------------------------------------------------------------
    # Add the types per dimension.
    #---------------------------------------------------------------------------
    def generateDimBindings(self, ndim):

        dim = "Spheral::Dim<%i>" % ndim
        vector = "Spheral::Dim<%i>::Vector" % ndim
        tensor = "Spheral::Dim<%i>::Tensor" % ndim
        symtensor = "Spheral::Dim<%i>::SymTensor" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::SymTensorFieldList%id" % ndim
        connectivitymap = "Spheral::ConnectivityMap%id" % ndim
        tablekernel = "Spheral::TableKernel%id" % ndim
        mesh = {1 : "LineMesh", 2 : "PolygonalMesh", 3 : "PolyhedralMesh"}[ndim]

        # SVPH sampling.
        for fl in [scalarfieldlist,
                   vectorfieldlist,
                   tensorfieldlist,
                   symtensorfieldlist]:
            self.space.add_function("sampleFieldListSVPH", 
                                    fl,
                                    [constrefparam(fl, "fieldList"),
                                     constrefparam(vectorfieldlist, "position"),
                                     constrefparam(symtensorfieldlist, "H"),
                                     constrefparam(connectivitymap, "connectivityMap"),
                                     constrefparam(tablekernel, "W"),
                                     constrefparam(mesh, "mesh"),
                                     param("bool", "firstOrderConsistent", default_value="true")],
                                    template_parameters = [dim],
                                    custom_name = "sampleFieldListSVPH%id" % ndim)

        # SVPH gradient.
        for (fl, gfl) in [(scalarfieldlist, vectorfieldlist),
                          (vectorfieldlist, tensorfieldlist)]:
            self.space.add_function("gradientFieldListSVPH", 
                                    gfl,
                                    [constrefparam(fl, "fieldList"),
                                     constrefparam(vectorfieldlist, "position"),
                                     constrefparam(symtensorfieldlist, "H"),
                                     constrefparam(connectivitymap, "connectivityMap"),
                                     constrefparam(tablekernel, "W"),
                                     constrefparam(mesh, "mesh"),
                                     param("bool", "firstOrderConsistent", default_value="true")],
                                    template_parameters = [dim],
                                    custom_name = "gradientFieldListSVPH%id" % ndim)

        return

    # #---------------------------------------------------------------------------
    # # Bindings (SVPHHydroBase).
    # #---------------------------------------------------------------------------
    # def generateSVPHHydroBaseBindings(self, x, ndim):

    #     # Object names.
    #     me = "Spheral::SVPHHydroBase%id" % ndim
    #     dim = "Spheral::Dim<%i>" % ndim
    #     vector = "Vector%id" % ndim
    #     tensor = "Tensor%id" % ndim
    #     symtensor = "SymTensor%id" % ndim
    #     fieldbase = "Spheral::FieldBase%id" % ndim
    #     intfield = "Spheral::IntField%id" % ndim
    #     scalarfield = "Spheral::ScalarField%id" % ndim
    #     vectorfield = "Spheral::VectorField%id" % ndim
    #     vector3dfield = "Spheral::Vector3dField%id" % ndim
    #     tensorfield = "Spheral::TensorField%id" % ndim
    #     thirdranktensorfield = "Spheral::ThirdRankTensorField%id" % ndim
    #     vectordoublefield = "Spheral::VectorDoubleField%id" % ndim
    #     vectorvectorfield = "Spheral::VectorVectorField%id" % ndim
    #     vectorsymtensorfield = "Spheral::VectorSymTensorField%id" % ndim
    #     symtensorfield = "Spheral::SymTensorField%id" % ndim
    #     intfieldlist = "Spheral::IntFieldList%id" % ndim
    #     scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
    #     vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
    #     vector3dfieldlist = "Spheral::Vector3dFieldList%id" % ndim
    #     tensorfieldlist = "Spheral::TensorFieldList%id" % ndim
    #     symtensorfieldlist = "Spheral::SymTensorFieldList%id" % ndim
    #     thirdranktensorfieldlist = "Spheral::ThirdRankTensorFieldList%id" % ndim
    #     vectordoublefieldlist = "Spheral::VectorDoubleFieldList%id" % ndim
    #     vectorvectorfieldlist = "Spheral::VectorVectorFieldList%id" % ndim
    #     vectorsymtensorfieldlist = "Spheral::VectorSymTensorFieldList%id" % ndim
    #     nodelist = "Spheral::NodeList%id" % ndim
    #     state = "Spheral::State%id" % ndim
    #     derivatives = "Spheral::StateDerivatives%id" % ndim
    #     database = "Spheral::DataBase%id" % ndim
    #     connectivitymap = "Spheral::ConnectivityMap%id" % ndim
    #     key = "pair_NodeList%id_string" % ndim
    #     vectorkeys = "vector_of_pair_NodeList%id_string" % ndim
    #     tablekernel = "Spheral::TableKernel%id" % ndim
    #     artificialviscosity = "Spheral::ArtificialViscosity%id" % ndim
    #     fileio = "Spheral::FileIO"
    #     smoothingscalebase = "Spheral::SmoothingScaleBase%id" % ndim
    #     mesh = {1 : "LineMesh", 2 : "PolygonalMesh", 3 : "PolyhedralMesh"}[ndim]

    #     # Constructors.
    #     x.add_constructor([constrefparam(smoothingscalebase, "smoothingScaleMethod"),
    #                        constrefparam(tablekernel, "W"),
    #                        refparam(artificialviscosity, "Q"),
    #                        param("double", "cfl"),
    #                        param("int", "useVelocityMagnitudeForDt"),
    #                        param("int", "compatibleEnergyEvolution"),
    #                        param("int", "XSVPH"),
    #                        param("int", "linearConsistent"),
    #                        param("MassDensityType", "densityUpdate"),
    #                        param("HEvolutionType", "HUpdate"),
    #                        param("double", "fcentroidal"),
    #                        param(vector, "xmin"),
    #                        param(vector, "xmax")])

    #     # Methods.
    #     x.add_method("initializeProblemStartup", None, [refparam(database, "dataBase")], is_virtual=True)
    #     x.add_method("registerState", None, [refparam(database, "dataBase"),
    #                                          refparam(state, "state")],
    #                  is_virtual=True)
    #     x.add_method("registerDerivatives", None, [refparam(database, "dataBase"),
    #                                                refparam(derivatives, "derivatives")],
    #                  is_virtual=True)
    #     x.add_method("initialize", None, [param("const double", "time"),
    #                                       param("const double", "dt"),
    #                                       constrefparam(database, "dataBase"),
    #                                       refparam(state, "state"),
    #                                       refparam(derivatives, "derivatives")], is_virtual=True)
    #     x.add_method("evaluateDerivatives", None, [param("const double", "time"),
    #                                                param("const double", "dt"),
    #                                                constrefparam(database, "dataBase"),
    #                                                constrefparam(state, "state"),
    #                                                refparam(derivatives, "derivatives")],
    #                  is_const=True, is_virtual=True)
    #     x.add_method("finalizeDerivatives", None, [param("const double", "time"),
    #                                                param("const double", "dt"),
    #                                                constrefparam(database, "dataBase"),
    #                                                constrefparam(state, "state"),
    #                                                refparam(derivatives, "derivatives")], is_const=True, is_virtual=True)
    #     x.add_method("finalize", None, [param("const double", "time"),
    #                                     param("const double", "dt"),
    #                                     refparam(database, "dataBase"),
    #                                     refparam(state, "state"),
    #                                     refparam(derivatives, "derivatives")], is_virtual=True)
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
    #     x.add_instance_attribute("densityUpdate", "MassDensityType", getter="densityUpdate", setter="densityUpdate")
    #     x.add_instance_attribute("HEvolution", "HEvolutionType", getter="HEvolution", setter="HEvolution")
    #     x.add_instance_attribute("compatibleEnergyEvolution", "bool", getter="compatibleEnergyEvolution", setter="compatibleEnergyEvolution")
    #     x.add_instance_attribute("XSVPH", "bool", getter="XSVPH", setter="XSVPH")
    #     x.add_instance_attribute("linearConsistent", "bool", getter="linearConsistent", setter="linearConsistent")
    #     x.add_instance_attribute("xmin", vector, getter="xmin", setter="xmin")
    #     x.add_instance_attribute("xmax", vector, getter="xmax", setter="xmax")
    #     x.add_instance_attribute("mesh", mesh, getter="mesh", is_const=True)

    #     const_ref_return_value(x, me, "%s::smoothingScaleMethod" % me, smoothingscalebase, [], "smoothingScaleMethod")
    #     const_ref_return_value(x, me, "%s::A" % me, scalarfieldlist, [], "A")
    #     const_ref_return_value(x, me, "%s::B" % me, vectorfieldlist, [], "B")
    #     const_ref_return_value(x, me, "%s::gradB" % me, tensorfieldlist, [], "gradB")
    #     const_ref_return_value(x, me, "%s::timeStepMask" % me, intfieldlist, [], "timeStepMask")
    #     const_ref_return_value(x, me, "%s::pressure" % me, scalarfieldlist, [], "pressure")
    #     const_ref_return_value(x, me, "%s::soundSpeed" % me, scalarfieldlist, [], "soundSpeed")
    #     const_ref_return_value(x, me, "%s::volume" % me, scalarfieldlist, [], "volume")
    #     const_ref_return_value(x, me, "%s::specificThermalEnergy0" % me, scalarfieldlist, [], "specificThermalEnergy0")
    #     const_ref_return_value(x, me, "%s::Hideal" % me, symtensorfieldlist, [], "Hideal")
    #     const_ref_return_value(x, me, "%s::maxViscousPressure" % me, scalarfieldlist, [], "maxViscousPressure")
    #     const_ref_return_value(x, me, "%s::massDensitySum" % me, scalarfieldlist, [], "massDensitySum")
    #     const_ref_return_value(x, me, "%s::weightedNeighborSum" % me, scalarfieldlist, [], "weightedNeighborSum")
    #     const_ref_return_value(x, me, "%s::massSecondMoment" % me, symtensorfieldlist, [], "massSecondMoment")
    #     const_ref_return_value(x, me, "%s::XSVPHDeltaV" % me, vectorfieldlist, [], "XSVPHDeltaV")
    #     const_ref_return_value(x, me, "%s::DxDt" % me, vectorfieldlist, [], "DxDt")
    #     const_ref_return_value(x, me, "%s::DvDt" % me, vectorfieldlist, [], "DvDt")
    #     const_ref_return_value(x, me, "%s::DmassDensityDt" % me, scalarfieldlist, [], "DmassDensityDt")
    #     const_ref_return_value(x, me, "%s::DspecificThermalEnergyDt" % me, scalarfieldlist, [], "DspecificThermalEnergyDt")
    #     const_ref_return_value(x, me, "%s::DHDt" % me, symtensorfieldlist, [], "DHDt")
    #     const_ref_return_value(x, me, "%s::DvDx" % me, tensorfieldlist, [], "DvDx")
    #     const_ref_return_value(x, me, "%s::internalDvDx" % me, tensorfieldlist, [], "internalDvDx")
    #     const_ref_return_value(x, me, "%s::pairAccelerations" % me, vectorvectorfieldlist, [], "pairAccelerations")

    #     return

    #---------------------------------------------------------------------------
    # Bindings (SVPHFacetedHydroBase).
    #---------------------------------------------------------------------------
    def generateSVPHFacetedHydroBaseBindings(self, x, ndim):

        # Object names.
        me = "Spheral::SVPHFacetedHydroBase%id" % ndim
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
        mesh = {1 : "LineMesh", 2 : "PolygonalMesh", 3 : "PolyhedralMesh"}[ndim]

        # Constructors.
        x.add_constructor([constrefparam(smoothingscalebase, "smoothingScaleMethod"),
                           constrefparam(tablekernel, "W"),
                           refparam(artificialviscosity, "Q"),
                           param("double", "cfl"),
                           param("int", "useVelocityMagnitudeForDt"),
                           param("int", "compatibleEnergyEvolution"),
                           param("int", "XSVPH"),
                           param("int", "linearConsistent"),
                           param("int", "generateVoid"),
                           param("MassDensityType", "densityUpdate"),
                           param("HEvolutionType", "HUpdate"),
                           param("double", "fcentroidal"),
                           param("double", "fcellPressure"),
                           param(vector, "xmin"),
                           param(vector, "xmax")])

        # Methods.
        x.add_method("dt", "pair_double_string", [constrefparam(database, "dataBase"),
                                                  constrefparam(state, "state"),
                                                  constrefparam(derivatives, "derivatives"),
                                                  param("double", "time")],
                     is_const=True, is_virtual=True)
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
        x.add_instance_attribute("densityUpdate", "MassDensityType", getter="densityUpdate", setter="densityUpdate")
        x.add_instance_attribute("HEvolution", "HEvolutionType", getter="HEvolution", setter="HEvolution")
        x.add_instance_attribute("compatibleEnergyEvolution", "bool", getter="compatibleEnergyEvolution", setter="compatibleEnergyEvolution")
        x.add_instance_attribute("XSVPH", "bool", getter="XSVPH", setter="XSVPH")
        x.add_instance_attribute("linearConsistent", "bool", getter="linearConsistent", setter="linearConsistent")
        x.add_instance_attribute("generateVoid", "bool", getter="generateVoid", setter="generateVoid")
        x.add_instance_attribute("fcentroidal", "double", getter="fcentroidal", setter="fcentroidal")
        x.add_instance_attribute("fcellPressure", "double", getter="fcellPressure", setter="fcellPressure")
        x.add_instance_attribute("xmin", vector, getter="xmin", setter="xmin")
        x.add_instance_attribute("xmax", vector, getter="xmax", setter="xmax")
        x.add_instance_attribute("mesh", mesh, getter="mesh", is_const=True)

        const_ref_return_value(x, me, "%s::smoothingScaleMethod" % me, smoothingscalebase, [], "smoothingScaleMethod")
        #const_ref_return_value(x, me, "%s::A" % me, scalarfieldlist, [], "A")
        #const_ref_return_value(x, me, "%s::B" % me, vectorfieldlist, [], "B")
        #const_ref_return_value(x, me, "%s::gradB" % me, tensorfieldlist, [], "gradB")
        const_ref_return_value(x, me, "%s::timeStepMask" % me, intfieldlist, [], "timeStepMask")
        const_ref_return_value(x, me, "%s::pressure" % me, scalarfieldlist, [], "pressure")
        const_ref_return_value(x, me, "%s::cellPressure" % me, scalarfieldlist, [], "cellPressure")
        const_ref_return_value(x, me, "%s::soundSpeed" % me, scalarfieldlist, [], "soundSpeed")
        const_ref_return_value(x, me, "%s::volume" % me, scalarfieldlist, [], "volume")
        const_ref_return_value(x, me, "%s::specificThermalEnergy0" % me, scalarfieldlist, [], "specificThermalEnergy0")
        const_ref_return_value(x, me, "%s::Hideal" % me, symtensorfieldlist, [], "Hideal")
        const_ref_return_value(x, me, "%s::maxViscousPressure" % me, scalarfieldlist, [], "maxViscousPressure")
        const_ref_return_value(x, me, "%s::massDensitySum" % me, scalarfieldlist, [], "massDensitySum")
        const_ref_return_value(x, me, "%s::weightedNeighborSum" % me, scalarfieldlist, [], "weightedNeighborSum")
        const_ref_return_value(x, me, "%s::massSecondMoment" % me, symtensorfieldlist, [], "massSecondMoment")
        const_ref_return_value(x, me, "%s::XSVPHDeltaV" % me, vectorfieldlist, [], "XSVPHDeltaV")
        const_ref_return_value(x, me, "%s::DxDt" % me, vectorfieldlist, [], "DxDt")
        const_ref_return_value(x, me, "%s::DvDt" % me, vectorfieldlist, [], "DvDt")
        const_ref_return_value(x, me, "%s::DmassDensityDt" % me, scalarfieldlist, [], "DmassDensityDt")
        const_ref_return_value(x, me, "%s::DspecificThermalEnergyDt" % me, scalarfieldlist, [], "DspecificThermalEnergyDt")
        const_ref_return_value(x, me, "%s::DHDt" % me, symtensorfieldlist, [], "DHDt")
        const_ref_return_value(x, me, "%s::DvDx" % me, tensorfieldlist, [], "DvDx")
        const_ref_return_value(x, me, "%s::internalDvDx" % me, tensorfieldlist, [], "internalDvDx")

        return

